#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_heavy_dynamics_benchmark.py
===============================
Phase 2 — Dynamic Biophysical Profiling & Force-Field Thermodynamics of 5
unindexed frontier therapeutic modalities (zero overlap with the Phase-1 set).

Pipeline:
  STAGE 1  Medicinal-chemistry triage: MW / cLogP / TPSA / Fsp3 / RotB /
           HBD / HBA + SAScore + QED + PAINS alerts (RDKit).
  STAGE 2  36-point relaxed torsional PES scan (0-360 deg, 10 deg steps)
           on each molecule's domain-hinge bond, MMFF94 (UFF fallback)
           constrained minimization per point, rotational barrier
           delta-E_barrier in kcal/mol.
  STAGE 3  OpenMM MD on the most flexible target: OpenFF Sage 2.1.0
           parameterization (via the dedicated `phase2ff` conda env +
           export_openff_system.py; MMFF94 charges when AM1-BCC backends
           are unavailable), GBSA/OBC2 implicit solvent, Langevin middle
           integrator at 300 K, gamma = 1/ps, dt = 2 fs, >= 100k steps,
           DCD every 500 steps; per-frame RMSD (Kabsch), radius of
           gyration (ASE), and intramolecular H-bond persistence.
  STAGE 4  Publication figures (300 DPI) in ./figures_phase2/.
  STAGE 5  Atomic JSON serialization + optional --auto_shutdown.

Usage:
    python run_heavy_dynamics_benchmark.py
    python run_heavy_dynamics_benchmark.py --md_steps 100000 --auto_shutdown
    python run_heavy_dynamics_benchmark.py --skip_md --skip_scan   # stage 1 + figures only
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

if sys.platform.startswith("win"):
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

import matplotlib

matplotlib.use("Agg")

SEED = 0x2026
RESULTS: Dict[str, Any] = {
    "phase": "2-dynamic-biophysical-profiling",
    "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
    "stage1_descriptors": {},
    "stage2_torsion_scans": {},
    "stage3_md": {},
    "meta": {"python": sys.version.split()[0]},
    "fatal_error": None,
    "all_stages_ok": False,
}

REGISTRY: List[Dict[str, str]] = [
    {"id": "T01", "name": "Molecular glue degrader (CRBN imide mimic)",
     "smiles": "O=C1NC(=O)C(N2C(=O)c3c(N4CCN(c5cccc(Cl)c5)CC4)c(F)ccc3C2=O)CC1"},
    {"id": "T02", "name": "ADC cathepsin-B linker (Val-Cit-PAB mimic)",
     "smiles": "CC(C)[C@@H](NC(=O)[C@H](CCCNC(=O)N)NC(=O)OCc1ccc(NC(=O)OCC2=CCN(CC2)C(=O)C)cc1)C(=O)O"},
    {"id": "T03", "name": "Bicyclic disulfide-constrained peptidomimetic",
     "smiles": "O=C1N[C@H]2CSSC[C@@H]3NC(=O)[C@H](CC(=O)N3)NC(=O)[C@H](CSSC2)NC1=O"},
    {"id": "T04", "name": "Allosteric covalent inhibitor core (switch-II pocket)",
     "smiles": "C=CC(=O)N1CCN(CC1)c2nc(Nc3ccc(F)c(Cl)c3)nc4c2c(C#C)c(c5ccccc5)n4C"},
    {"id": "T05", "name": "Non-thalidomide E3 binder (VHL-mimetic proline hybrid)",
     "smiles": "CC(C)(C)[C@H](NC(=O)[C@@H]1C[C@@H](O)CN1C(=O)c2ccc(c3nccs3)cc2)C(=O)NCc4ccc(C#N)cc4"},
]

BONDI_NM = {1: 0.120, 6: 0.170, 7: 0.155, 8: 0.152, 9: 0.147, 16: 0.180, 17: 0.175}
OBC_SCALE = {1: 0.85, 6: 0.72, 7: 0.79, 8: 0.85, 9: 0.88, 16: 0.96, 17: 0.80}


def _log(stage: str, msg: str) -> None:
    stamp = _dt.datetime.now().strftime("%H:%M:%S")
    print(f"[{stamp}] [{stage}] {msg}", flush=True)


def _hr(title: str) -> None:
    print(f"\n{'=' * 72}\n{title}\n{'=' * 72}", flush=True)


def write_json_atomic(path: Path, payload: Dict[str, Any]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str),
                   encoding="utf-8")
    os.replace(tmp, path)


# --------------------------------------------------------------------------- #
# STAGE 1 — medicinal chemistry triage
# --------------------------------------------------------------------------- #
def _load_sascorer():
    try:
        from rdkit.Contrib.SA_Score import sascorer
        return sascorer
    except ImportError:
        from rdkit import RDConfig
        sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
        import sascorer
        return sascorer


def stage1_descriptors() -> None:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, FilterCatalog, Lipinski, QED, rdMolDescriptors

    sascorer = _load_sascorer()
    fc_params = FilterCatalog.FilterCatalogParams()
    for catalog in (FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A,
                    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B,
                    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C):
        fc_params.AddCatalog(catalog)
    catalog = FilterCatalog.FilterCatalog(fc_params)

    _hr("STAGE 1 — medicinal chemistry & structural triage")
    for entry in REGISTRY:
        mol = Chem.MolFromSmiles(entry["smiles"])
        if mol is None:
            RESULTS["stage1_descriptors"][entry["id"]] = {"status": "parse_failed"}
            _log("stage1", f"{entry['id']}: PARSE FAILED")
            continue
        matches = catalog.GetMatches(mol)
        rec = {
            "status": "ok",
            "name": entry["name"],
            "mw": round(Descriptors.MolWt(mol), 2),
            "clogp": round(Crippen.MolLogP(mol), 2),
            "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "fsp3": round(Lipinski.FractionCSP3(mol), 3),
            "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
            "hbd": int(Lipinski.NumHDonors(mol)),
            "hba": int(Lipinski.NumHAcceptors(mol)),
            "sascore": round(float(sascorer.calculateScore(mol)), 2),
            "qed": round(float(QED.qed(mol)), 3),
            "pains_alerts": len(matches),
            "pains_families": [m.GetDescription() for m in matches][:3],
        }
        RESULTS["stage1_descriptors"][entry["id"]] = rec
        _log("stage1", f"{entry['id']}: MW={rec['mw']} cLogP={rec['clogp']} "
                       f"TPSA={rec['tpsa']} Fsp3={rec['fsp3']} RotB={rec['rotatable_bonds']} "
                       f"SAS={rec['sascore']} QED={rec['qed']} PAINS={rec['pains_alerts']}")


# --------------------------------------------------------------------------- #
# STAGE 2 — relaxed torsional PES scans
# --------------------------------------------------------------------------- #
def _pick_scan_torsion(mol: Any) -> Dict[str, Any]:
    """Pick the dihedral (i,j,k,l) spanning the molecule's 'domain hinge'.

    Heuristic: among heavy-atom single bonds that admit a dihedral, rank by
    the size of the smaller fragment produced on bond deletion (the bond that
    splits the molecule into its two largest halves is the linker hinge).
    Ring bonds score 0; fall back to an S-S bond (disulfide strain probe) or
    the most central single bond.
    """
    from rdkit import Chem

    best: Optional[Tuple[Tuple[int, int, int, int], float, int, Any]] = None
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        n1 = [nb for nb in a1.GetNeighbors()
              if nb.GetIdx() != a2.GetIdx() and nb.GetAtomicNum() > 1]
        n2 = [nb for nb in a2.GetNeighbors()
              if nb.GetIdx() != a1.GetIdx() and nb.GetAtomicNum() > 1]
        if not n1 or not n2:
            continue
        if bond.IsInRing():
            score, split = 0.0, 0
        else:
            rw = Chem.RWMol(mol)
            rw.RemoveBond(a1.GetIdx(), a2.GetIdx())
            frags = [len(f) for f in Chem.GetMolFrags(rw.GetMol())]
            split = min(frags) if len(frags) >= 2 else 0
            score = float(split) - 0.001 * abs(frags[0] - frags[1]) if len(frags) >= 2 else 0.0
        i = n1[0].GetIdx()
        k = a1.GetIdx()
        j = a2.GetIdx()
        l = n2[0].GetIdx()
        if score > 0 and (best is None or score > best[1]):
            best = ((i, k, j, l), score, split, bond)

    if best is not None and best[1] > 0:
        return {"torsion": best[0], "mode": "hinge (rotatable, max-split)",
                "split_sizes": best[2]}

    # fallback A: disulfide S-S bond (strain probe), fallback B: any candidate
    ss = None
    generic = None
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        n1 = [nb for nb in a1.GetNeighbors()
              if nb.GetIdx() != a2.GetIdx() and nb.GetAtomicNum() > 1]
        n2 = [nb for nb in a2.GetNeighbors()
              if nb.GetIdx() != a1.GetIdx() and nb.GetAtomicNum() > 1]
        if not n1 or not n2:
            continue
        cand = ((n1[0].GetIdx(), a1.GetIdx(), a2.GetIdx(), n2[0].GetIdx()), bond)
        if a1.GetSymbol() == "S" and a2.GetSymbol() == "S":
            ss = cand
        generic = generic or cand
    chosen = ss or generic
    if chosen is None:
        raise RuntimeError("no eligible dihedral for torsional scan")
    return {"torsion": chosen[0],
            "mode": "disulfide S-S strain probe" if ss else "central single bond",
            "split_sizes": None}


def stage2_torsion_scans(out_dir: Path, num_points: int) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    _hr("STAGE 2 — 36-point relaxed torsional PES scans (MMFF94)")
    for entry in REGISTRY:
        tid = entry["id"]
        try:
            mol = Chem.MolFromSmiles(entry["smiles"])
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            mh = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = SEED
            params.useSmallRingTorsions = True
            params.useMacrocycleTorsions = True
            if AllChem.EmbedMolecule(mh, params) < 0:
                retry = AllChem.ETKDGv3()
                retry.useRandomCoords = True
                retry.randomSeed = SEED + 1
                if AllChem.EmbedMolecule(mh, retry) < 0:
                    raise RuntimeError("ETKDGv3 embedding failed")

            props = AllChem.MMFFGetMoleculeProperties(mh)
            mmff_ok = props is not None

            def make_ff(m):
                if mmff_ok:
                    return AllChem.MMFFGetMoleculeForceField(m, props)
                return AllChem.UFFGetMoleculeForceField(m)

            def add_torsion_constraint(ff, i, j, k, l, lo, hi):
                if mmff_ok:
                    ff.MMFFAddTorsionConstraint(i, j, k, l, False, lo, hi, 1000.0)
                else:
                    ff.UFFAddTorsionConstraint(i, j, k, l, False, lo, hi, 1000.0)

            scan = _pick_scan_torsion(mh)
            i, j, k, l = scan["torsion"]
            conf = mh.GetConformer()

            angles: List[float] = []
            energies: List[float] = []
            achieved: List[float] = []
            for n in range(num_points):
                theta = 360.0 * n / num_points
                try:
                    AllChem.SetDihedralDeg(conf, i, j, k, l, theta)
                except ValueError:
                    pass  # in-ring torsion: flat-bottom constraint alone drives it
                ff = make_ff(mh)
                add_torsion_constraint(ff, i, j, k, l, theta - 5.0, theta + 5.0)
                ff.Minimize(800)
                energies.append(float(ff.CalcEnergy()))
                angles.append(theta)
                achieved.append(AllChem.GetDihedralDeg(conf, i, j, k, l))

            e_min = min(energies)
            rel = [e - e_min for e in energies]
            barrier = max(rel)
            i_min = rel.index(0.0)
            i_max = rel.index(barrier)

            # persist scan table + minimum conformer
            csv_path = out_dir / f"{tid}_torsion_scan.csv"
            lines = ["angle_deg,energy_kcal_mol,rel_energy_kcal_mol,achieved_deg"]
            lines += [f"{a:.1f},{e:.4f},{r:.4f},{g:.2f}"
                      for a, e, r, g in zip(angles, energies, rel, achieved)]
            csv_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

            writer = Chem.SDWriter(str(out_dir / f"{tid}_scan_min.sdf"))
            try:
                mh.SetProp("_Name", f"{tid} scan-minimum")
                mh.SetProp("dihedral_deg_min", f"{angles[i_min]:.1f}")
                writer.write(mh, confId=conf.GetId())
            finally:
                writer.close()

            RESULTS["stage2_torsion_scans"][tid] = {
                "status": "ok",
                "force_field": "MMFF94" if mmff_ok else "UFF",
                "scan_mode": scan["mode"],
                "torsion_atoms": [i, j, k, l],
                "n_points": num_points,
                "angle_of_minimum_deg": angles[i_min],
                "angle_of_barrier_deg": angles[i_max],
                "e_min_kcal_mol": round(e_min, 2),
                "delta_e_barrier_kcal_mol": round(barrier, 2),
                "achieved_max_dev_deg": round(max(abs(a - g) for a, g in zip(angles, achieved)), 2),
                "angles": angles,
                "rel_energies": [round(r, 3) for r in rel],
                "csv": str(csv_path),
            }
            _log("stage2", f"{tid}: barrier={barrier:.2f} kcal/mol "
                           f"(min@{angles[i_min]:.0f}°, max@{angles[i_max]:.0f}°, "
                           f"{scan['mode']})")
        except Exception as exc:
            RESULTS["stage2_torsion_scans"][tid] = {"status": "failed",
                                                    "error": f"{type(exc).__name__}: {exc}"}
            _log("stage2", f"{tid}: FAILED — {type(exc).__name__}: {exc}")


# --------------------------------------------------------------------------- #
# STAGE 3 — OpenMM MD (Sage 2.1.0 + GBSA/OBC2)
# --------------------------------------------------------------------------- #
def _kabsch_rmsd(p: Any, q: Any) -> float:
    """RMSD (same atom order) after optimal rigid-body superposition."""
    import numpy as np

    pc = p - p.mean(axis=0)
    qc = q - q.mean(axis=0)
    v, _s, wt = np.linalg.svd(pc.T @ qc)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    diff = pc @ rot - qc
    return float(np.sqrt((diff ** 2).sum(axis=1).mean()))


def _topology_from_rdkit(mol: Any):
    """Deterministic openmm Topology matching RDKit atom order."""
    from openmm import app

    topo = app.Topology()
    chain = topo.addChain()
    residue = topo.addResidue("LIG", chain)
    atoms = []
    for atom in mol.GetAtoms():
        element = app.Element.getByAtomicNumber(atom.GetAtomicNum())
        atoms.append(topo.addAtom(atom.GetSymbol(), element, residue))
    for bond in mol.GetBonds():
        topo.addBond(atoms[bond.GetBeginAtomIdx()], atoms[bond.GetEndAtomIdx()])
    return topo


class _AnalysisReporter:
    """OpenMM reporter: per-frame RMSD / Rg (ASE) / energies / IMHB census."""

    def __init__(self, interval: int, heavy_idx, heavy_masses, heavy_numbers,
                 ref_positions_nm, donors_h, acceptors, dof: int):
        import numpy as np
        from ase import Atoms
        from ase.data import atomic_masses

        self.interval = interval
        self._np = np
        self.heavy_idx = self._np.asarray(heavy_idx)
        self.ref = ref_positions_nm[self.heavy_idx]
        mass_weights = self._np.asarray([atomic_masses[z] for z in heavy_numbers])
        self.atoms = Atoms(numbers=list(heavy_numbers), masses=list(mass_weights),
                           positions=self._np.zeros((len(heavy_idx), 3)))
        self.donors_h = donors_h           # list[(d_idx, [h_idx, ...])]
        self.acceptors = self._np.asarray(acceptors)
        self.dof = dof
        self.rows: List[Dict[str, float]] = []
        self.hbond_counts: Dict[Tuple[int, int], int] = {}
        self._kB = 0.008314462618  # kJ/mol/K

    def describeNextReport(self, simulation):
        steps = self.interval - simulation.currentStep % self.interval
        # (steps, positions, velocities, forces, energy, enforcePeriodicBox)
        return (steps, True, False, False, True, False)

    def _imhb(self, pos_nm) -> List[Tuple[int, int]]:
        hits = []
        np = self._np
        for d_idx, h_list in self.donors_h:
            for a_idx in self.acceptors:
                if a_idx == d_idx:
                    continue
                dist_da = np.linalg.norm(pos_nm[d_idx] - pos_nm[a_idx])
                if dist_da > 0.35:  # 3.5 A
                    continue
                for h_idx in h_list:
                    v1 = pos_nm[d_idx] - pos_nm[h_idx]
                    v2 = pos_nm[a_idx] - pos_nm[h_idx]
                    cosang = float(np.dot(v1, v2) /
                                   (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12))
                    if math.degrees(math.acos(max(-1.0, min(1.0, cosang)))) >= 120.0:
                        hits.append((int(d_idx), int(a_idx)))
                        break
        return hits

    def report(self, simulation, state):
        import openmm.unit as unit

        np = self._np
        pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        self.atoms.set_positions(pos_nm[self.heavy_idx] * 10.0)  # Angstrom
        com = self.atoms.get_center_of_mass()
        rg = float(np.sqrt((self.atoms.get_masses() *
                            ((self.atoms.get_positions() - com) ** 2).sum(axis=1)
                            ).sum() / self.atoms.get_masses().sum()))
        rmsd = _kabsch_rmsd(pos_nm[self.heavy_idx], self.ref) * 10.0  # A
        hbonds = self._imhb(pos_nm)
        for pair in hbonds:
            self.hbond_counts[pair] = self.hbond_counts.get(pair, 0) + 1
        temp = 2.0 * ke / (self.dof * self._kB)
        self.rows.append({
            "step": int(simulation.currentStep),
            "time_ps": simulation.currentStep * 0.002,
            "pe_kj_mol": round(pe, 2),
            "ke_kj_mol": round(ke, 2),
            "te_kj_mol": round(pe + ke, 2),
            "temp_k": round(temp, 1),
            "rmsd_A": round(rmsd, 3),
            "rg_A": round(rg, 3),
            "n_imhb": len(hbonds),
        })


def stage3_openmm_md(out_dir: Path, args: argparse.Namespace) -> None:
    import openmm
    import openmm.unit as unit
    from openmm import app
    from rdkit import Chem

    _hr("STAGE 3 — OpenMM MD: OpenFF Sage 2.1.0 + GBSA/OBC2, 300 K")

    # --- target selection: most conformationally flexible entry ----------
    flex = {tid: rec.get("rotatable_bonds", -1)
            for tid, rec in RESULTS["stage1_descriptors"].items()}
    target_id = args.md_target if args.md_target != "auto" else max(flex, key=flex.get)
    entry = next(e for e in REGISTRY if e["id"] == target_id)
    _log("stage3", f"MD target = {target_id} ({entry['name']}; RotB={flex.get(target_id)})")

    prefix = str(out_dir / target_id)
    export_script = Path(__file__).resolve().parent / "export_openff_system.py"
    if not export_script.exists():
        raise RuntimeError(f"export script missing: {export_script}")

    ff_env = args.ff_env
    if not Path(ff_env).exists():
        candidates = [Path(sys.prefix) / "envs" / "phase2ff" / "python.exe",
                      Path(sys.prefix).parent / "envs" / "phase2ff" / "python.exe",
                      Path.home() / "miniconda3" / "envs" / "phase2ff" / "python.exe"]
        ff_env = next((str(c) for c in candidates if c.exists()), ff_env)
    if not Path(ff_env).exists():
        raise RuntimeError(f"phase2ff parameterization env not found (tried {ff_env}); "
                           f"pass --ff_env /path/to/python.exe")

    _log("stage3", f"parameterizing via {ff_env} ...")
    proc = subprocess.run([ff_env, str(export_script), entry["smiles"], prefix],
                          capture_output=True, text=True, errors="replace", timeout=1800)
    if proc.returncode != 0:
        raise RuntimeError(f"export failed: {(proc.stderr or proc.stdout)[-800:]}")
    meta = json.loads(Path(f"{prefix}_meta.json").read_text(encoding="utf-8"))
    _log("stage3", f"system exported: {meta['n_particles']} particles, "
                   f"{meta['n_constraints']} constraints, charges={meta['charge_source']}")

    # --- assemble in the main runtime ------------------------------------
    system = openmm.XmlSerializer.deserialize(Path(f"{prefix}_system.xml").read_text())
    pdb = app.PDBFile(f"{prefix}_start.pdb")

    mol = Chem.MolFromSmiles(entry["smiles"])
    mh = Chem.AddHs(mol)
    if mh.GetNumAtoms() != system.getNumParticles():
        raise RuntimeError("atom-count mismatch between RDKit and exported system")
    topology = _topology_from_rdkit(mh)
    positions = pdb.positions

    nonbonded = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    if nonbonded.getNonbondedMethod() == openmm.NonbondedForce.CutoffPeriodic:
        nonbonded.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    def _charge_value(q):
        try:
            return float(q.value_in_unit(unit.elementary_charge))
        except AttributeError:
            return float(q)

    charges = [_charge_value(nonbonded.getParticleParameters(i)[0])
               for i in range(system.getNumParticles())]

    gb = openmm.GBSAOBCForce()
    gb.setSolventDielectric(78.5)
    gb.setSoluteDielectric(1.0)
    for idx, atom in enumerate(topology.atoms()):
        z = atom.element.atomic_number
        gb.addParticle(charges[idx], BONDI_NM.get(z, 0.15), OBC_SCALE.get(z, 0.8))
    system.addForce(gb)
    _log("stage3", f"GBSA/OBC2 added (solvent eps 78.5, solute eps 1.0, "
                   f"{gb.getNumParticles()} particles)")

    integrator = openmm.LangevinMiddleIntegrator(300.0 * unit.kelvin,
                                                 1.0 / unit.picosecond,
                                                 2.0 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)

    t0 = time.time()
    simulation.minimizeEnergy(maxIterations=1000)
    _log("stage3", f"energy minimized in {time.time() - t0:.1f} s")

    state0 = simulation.context.getState(getPositions=True)
    ref_nm = state0.getPositions(asNumpy=True).value_in_unit(unit.nanometer)

    heavy_idx = [a.GetIdx() for a in mh.GetAtoms() if a.GetAtomicNum() > 1]
    heavy_numbers = [mh.GetAtomWithIdx(ix).GetAtomicNum() for ix in heavy_idx]

    donors_h: List[Tuple[int, List[int]]] = []
    acceptors: List[int] = []
    for atom in mh.GetAtoms():
        z = atom.GetAtomicNum()
        if z not in (7, 8):
            continue
        h_ids = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 1]
        if h_ids:
            donors_h.append((atom.GetIdx(), h_ids))
        q = charges[atom.GetIdx()]
        if q < 0.30:  # neutral/anionic heteroatoms can accept
            acceptors.append(atom.GetIdx())
    dof = 3 * system.getNumParticles() - system.getNumConstraints() - 3

    simulation.context.setVelocitiesToTemperature(300.0 * unit.kelvin)
    _log("stage3", f"equilibration: {args.equil_steps} steps ...")
    simulation.step(args.equil_steps)

    dcd_path = f"{prefix}_trajectory.dcd"
    analysis = _AnalysisReporter(args.report_interval, heavy_idx, None, heavy_numbers,
                                 ref_nm, donors_h, acceptors, dof)
    simulation.reporters.append(app.DCDReporter(dcd_path, args.report_interval))
    simulation.reporters.append(analysis)
    simulation.reporters.append(app.StateDataReporter(
        sys.stdout, 5000, step=True, time=True, speed=True, temperature=True,
        potentialEnergy=True, remainingTime=True, totalSteps=args.md_steps))

    _log("stage3", f"production: {args.md_steps} steps "
                   f"({args.md_steps * 0.002:.0f} ps), recording every {args.report_interval}")
    t0 = time.time()
    simulation.step(args.md_steps)
    elapsed = time.time() - t0
    _log("stage3", f"production finished in {elapsed:.1f} s "
                   f"({args.md_steps / elapsed:.0f} steps/s)")

    final_state = simulation.context.getState(getPositions=True)
    with open(f"{prefix}_final.pdb", "w") as fh:
        app.PDBFile.writeFile(topology, final_state.getPositions(), fh, keepIds=True)

    rows = analysis.rows
    n_frames = len(rows)
    csv_path = Path(f"{prefix}_md_metrics.csv")
    header = "step,time_ps,pe_kj_mol,ke_kj_mol,te_kj_mol,temp_k,rmsd_A,rg_A,n_imhb"
    csv_path.write_text(
        "\n".join([header] + [",".join(str(r[c]) for c in header.split(",")) for r in rows]) + "\n",
        encoding="utf-8")

    persistence = sorted(
        ({"donor": d, "acceptor": a, "persistence": c / n_frames}
         for (d, a), c in analysis.hbond_counts.items()),
        key=lambda h: -h["persistence"])

    import statistics as _st
    rmsd_tail = [r["rmsd_A"] for r in rows[len(rows) // 2:]] or [0.0]
    rg_all = [r["rg_A"] for r in rows] or [0.0]
    RESULTS["stage3_md"] = {
        "target": target_id,
        "name": entry["name"],
        "force_field": "OpenFF Sage 2.1.0 (valence+vdW) + MMFF94 charges + GBSA/OBC2",
        "charge_source": meta["charge_source"],
        "integrator": "LangevinMiddleIntegrator 300 K, gamma=1/ps, dt=2 fs",
        "n_particles": system.getNumParticles(),
        "n_constraints": system.getNumConstraints(),
        "equil_steps": args.equil_steps,
        "production_steps": args.md_steps,
        "production_ps": args.md_steps * 0.002,
        "steps_per_second": round(args.md_steps / elapsed, 1),
        "n_frames": n_frames,
        "rmsd_mean_A": round(_st.mean(rmsd_tail), 3),
        "rmsd_max_A": round(max(r["rmsd_A"] for r in rows), 3),
        "rg_mean_A": round(_st.mean(rg_all), 3),
        "rg_std_A": round(_st.pstdev(rg_all), 3),
        "temp_mean_k": round(_st.mean([r["temp_k"] for r in rows]), 1),
        "pe_mean_kj_mol": round(_st.mean([r["pe_kj_mol"] for r in rows]), 1),
        "imhb_mean_per_frame": round(_st.mean([r["n_imhb"] for r in rows]), 2),
        "imhb_persistent_pairs": [p for p in persistence if p["persistence"] >= 0.10][:8],
        "dcd": dcd_path,
        "metrics_csv": str(csv_path),
        "final_pdb": f"{prefix}_final.pdb",
        "rows": rows,
    }
    _log("stage3", f"MD metrics: <RMSD>={RESULTS['stage3_md']['rmsd_mean_A']} A "
                   f"(tail), <Rg>={RESULTS['stage3_md']['rg_mean_A']} +/- "
                   f"{RESULTS['stage3_md']['rg_std_A']} A, "
                   f"T={RESULTS['stage3_md']['temp_mean_k']} K, "
                   f"IMHB/frame={RESULTS['stage3_md']['imhb_mean_per_frame']}")


# --------------------------------------------------------------------------- #
# STAGE 4 — publication figures
# --------------------------------------------------------------------------- #
def stage4_figures(fig_dir: Path) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    fig_dir.mkdir(parents=True, exist_ok=True)

    # ---- fig 1: torsional scans -----------------------------------------
    scans = RESULTS["stage2_torsion_scans"]
    fig, axes = plt.subplots(3, 2, figsize=(12.5, 13.5))
    slots = [ax for row in axes for ax in row]
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(REGISTRY)))
    barriers = []
    for slot, entry, color in zip(slots, REGISTRY, colors):
        rec = scans.get(entry["id"], {})
        if rec.get("status") != "ok":
            slot.set_axis_off()
            slot.text(0.5, 0.5, f"{entry['id']}: scan failed", ha="center", va="center",
                      transform=slot.transAxes, color="firebrick")
            continue
        angles = np.asarray(rec["angles"])
        rel = np.asarray(rec["rel_energies"])
        order = np.argsort(angles % 360)
        slot.plot(np.concatenate([angles[order], angles[order] + 360]),
                  np.concatenate([rel[order], rel[order]]), "-o", color=color,
                  markersize=3.5, linewidth=1.4)
        i_min = int(np.argmin(rel))
        i_max = int(np.argmax(rel))
        slot.plot(angles[i_min], rel[i_min], "*", color="seagreen", markersize=15,
                  zorder=5, label=f"min @ {angles[i_min]:.0f}°")
        slot.plot(angles[i_max], rel[i_max], "*", color="firebrick", markersize=15,
                  zorder=5, label=f"barrier @ {angles[i_max]:.0f}°")
        slot.annotate(f"ΔE‡ = {rec['delta_e_barrier_kcal_mol']:.1f} kcal/mol",
                      xy=(angles[i_max], rel[i_max]),
                      xytext=(0.03, 0.90), textcoords="axes fraction",
                      fontsize=10, fontweight="bold", color="firebrick")
        slot.set_title(f"{entry['id']} — {entry['name']}", fontsize=10.5)
        slot.set_xlabel("dihedral angle (°)")
        slot.set_ylabel("ΔE (kcal/mol)")
        slot.set_xticks(range(0, 361, 60))
        slot.set_ylim(top=max(rel) * 1.12)  # headroom for the barrier star
        slot.legend(fontsize=8, loc="upper right")
        slot.grid(alpha=0.35)
        barriers.append((entry["id"], rec["delta_e_barrier_kcal_mol"]))
    summary = slots[len(REGISTRY)]
    if barriers:
        ids = [b[0] for b in barriers]
        vals = [b[1] for b in barriers]
        bars = summary.bar(ids, vals, color=colors[:len(ids)], edgecolor="black", linewidth=0.6)
        for bar, val in zip(bars, vals):
            summary.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                         f"{val:.1f}", ha="center", fontsize=9, fontweight="bold")
        summary.set_ylabel("ΔE‡ barrier (kcal/mol)")
        summary.set_title("Rotational barriers across modalities", fontsize=10.5)
        summary.grid(alpha=0.35, axis="y")
    else:
        summary.set_axis_off()
    fig.suptitle("Relaxed Torsional Energy Landscapes (36-point MMFF94 scans)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out1 = fig_dir / "fig1_torsion_scans.png"
    fig.savefig(out1, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log("stage4", f"wrote {out1}")

    # ---- fig 2: MD trajectories ------------------------------------------
    md = RESULTS["stage3_md"]
    if md.get("rows"):
        rows = md["rows"]
        t = [r["time_ps"] for r in rows]
        fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
        axes[0].plot(t, [r["rmsd_A"] for r in rows], color="#0173B2", linewidth=0.9)
        tail = float(md["rmsd_mean_A"])
        axes[0].axhline(tail, color="firebrick", linestyle="--", linewidth=1,
                        label=f"tail mean = {tail:.2f} Å")
        axes[0].set_ylabel("RMSD (Å)")
        axes[0].set_title(f"A — RMSD vs initial relaxed state — {md['target']} "
                          f"({md['production_ps']:.0f} ps, OBC2, 300 K)", fontsize=11)
        axes[0].legend(fontsize=9)

        axes[1].plot(t, [r["rg_A"] for r in rows], color="#D55E00", linewidth=0.9)
        axes[1].set_ylabel("R$_g$ (Å)")
        axes[1].set_title(f"B — Radius of gyration (mean {md['rg_mean_A']:.2f} "
                          f"± {md['rg_std_A']:.2f} Å)", fontsize=11)

        axes[2].plot(t, [r["pe_kj_mol"] for r in rows], color="#0173B2",
                     linewidth=0.8, label="Potential")
        axes[2].plot(t, [r["ke_kj_mol"] for r in rows], color="#009E73",
                     linewidth=0.8, label="Kinetic")
        axes[2].plot(t, [r["te_kj_mol"] for r in rows], color="black",
                     linewidth=1.0, label="Total")
        axes[2].set_ylabel("Energy (kJ/mol)")
        axes[2].set_xlabel("time (ps)")
        axes[2].set_title("C — Energy components", fontsize=11)
        axes[2].legend(fontsize=9, ncol=3)
        for ax in axes:
            ax.grid(alpha=0.35)
        fig.tight_layout()
        out2 = fig_dir / "fig2_md_trajectories.png"
        fig.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        _log("stage4", f"wrote {out2}")

    # ---- fig 3: medchem radar --------------------------------------------
    descs = RESULTS["stage1_descriptors"]
    radar_axes = ["SAScore", "QED", "Fsp³", "TPSA (Å²)", "MW (Da)"]
    keys = ["sascore", "qed", "fsp3", "tpsa", "mw"]
    valid = [e for e in REGISTRY if descs.get(e["id"], {}).get("status") == "ok"]
    matrix = np.array([[descs[e["id"]][k] for k in keys] for e in valid], dtype=float)
    maxima = matrix.max(axis=0)
    maxima[maxima == 0] = 1.0

    angles = np.linspace(0, 2 * np.pi, len(radar_axes), endpoint=False).tolist()
    angles += angles[:1]
    fig, ax = plt.subplots(figsize=(9.2, 8.0), subplot_kw=dict(polar=True))
    palette = ["#0173B2", "#D55E00", "#009E73", "#CC78BC", "#E63312"]
    short_names = {"T01": "CRBN molecular glue", "T02": "ADC Val-Cit-PAB linker",
                   "T03": "bicyclic disulfide peptide", "T04": "covalent inhibitor core",
                   "T05": "VHL-mimetic hybrid"}
    for entry, color in zip(valid, palette):
        vals = (np.array([descs[entry["id"]][k] for k in keys]) / maxima).tolist()
        vals += vals[:1]
        ax.plot(angles, vals, linewidth=2.0,
                label=f"{entry['id']} · {short_names[entry['id']]}",
                color=color)
        ax.fill(angles, vals, alpha=0.10, color=color)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(radar_axes, fontsize=11)
    ax.set_rgrids([0.25, 0.50, 0.75, 1.00], angle=90, fontsize=8, color="0.45")
    ax.set_ylim(0, 1.02)
    ax.set_title("Drug-likeness & Synthesizability Fingerprint\n"
                 "(axes normalized to set maximum; SAScore: larger = harder to synthesize)",
                 fontsize=11.5, fontweight="bold", pad=26)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=2,
              fontsize=9.5, frameon=False)
    fig.tight_layout()
    out3 = fig_dir / "fig3_medchem_radar.png"
    fig.savefig(out3, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log("stage4", f"wrote {out3}")


# --------------------------------------------------------------------------- #
# STAGE 5 — shutdown protocol (opt-in)
# --------------------------------------------------------------------------- #
def schedule_shutdown(delay_s: int) -> Dict[str, Any]:
    info: Dict[str, Any] = {"scheduled": False, "delay_seconds": int(delay_s)}
    try:
        if sys.platform.startswith("win"):
            cmd = ["shutdown", "/s", "/t", str(int(delay_s)),
                   "/c", "run_heavy_dynamics_benchmark.py finished."]
            info["cancel_command"] = "shutdown /a"
        else:
            cmd = ["shutdown", "-h", f"+{max(1, int(round(delay_s / 60)))}"]
            info["cancel_command"] = "sudo shutdown -c"
        proc = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=30)
        info["command"] = " ".join(cmd)
        info["scheduled"] = proc.returncode == 0
        if not info["scheduled"]:
            info["error"] = (proc.stderr or proc.stdout or "").strip()
    except Exception as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"
    if info["scheduled"]:
        print("\n" + "!" * 62)
        print(f"  SYSTEM SHUTDOWN SCHEDULED -- power-off in ~{delay_s} s")
        print(f"  cancel with : {info.get('cancel_command')}")
        print("!" * 62)
    return info


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Phase-2 heavy dynamics benchmark "
                                            "(torsion scans + OpenMM MD + figures)")
    p.add_argument("--out_dir", default="results_phase2")
    p.add_argument("--fig_dir", default="figures_phase2")
    p.add_argument("--ff_env", default="phase2ff-missing",
                   help="python of the OpenFF parameterization env")
    p.add_argument("--scan_points", type=int, default=36)
    p.add_argument("--md_target", default="auto",
                   help="molecule id for MD (default: highest RotB)")
    p.add_argument("--equil_steps", type=int, default=5000)
    p.add_argument("--md_steps", type=int, default=100000)
    p.add_argument("--report_interval", type=int, default=500)
    p.add_argument("--skip_scan", action="store_true")
    p.add_argument("--skip_md", action="store_true")
    p.add_argument("--fig_only", action="store_true",
                   help="regenerate figures from an existing phase2_results.json")
    p.add_argument("--auto_shutdown", action="store_true",
                   help="clean power-off --shutdown_delay s after a completed run "
                        "(default False)")
    p.add_argument("--shutdown_delay", type=int, default=60)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exit_code = 0
    try:
        _hr(f"PHASE 2 HEAVY DYNAMICS BENCHMARK — {RESULTS['timestamp']}")
        _log("main", f"python {sys.version.split()[0]} | out_dir={out_dir} | "
                     f"fig_dir={args.fig_dir}")

        if args.fig_only:
            payload = json.loads((out_dir / "phase2_results.json").read_text(encoding="utf-8"))
            RESULTS["stage1_descriptors"] = payload.get("stage1_descriptors", {})
            RESULTS["stage2_torsion_scans"] = payload.get("stage2_torsion_scans", {})
            RESULTS["stage3_md"] = payload.get("stage3_md", {})
            _log("main", "--fig_only: loaded previous results, regenerating figures")
            stage4_figures(Path(args.fig_dir))
            return 0

        stage1_descriptors()
        RESULTS["all_stages_ok"] = True

        if not args.skip_scan:
            stage2_torsion_scans(out_dir, args.scan_points)
            if any(r.get("status") != "ok" for r in RESULTS["stage2_torsion_scans"].values()):
                RESULTS["all_stages_ok"] = False

        if not args.skip_md:
            stage3_openmm_md(out_dir, args)
            if RESULTS["stage3_md"].get("n_frames", 0) < args.md_steps / args.report_interval - 2:
                RESULTS["all_stages_ok"] = False

        stage4_figures(Path(args.fig_dir))
    except Exception as exc:
        RESULTS["fatal_error"] = f"{type(exc).__name__}: {exc}"
        RESULTS["all_stages_ok"] = False
        print(f"\n[FATAL] {RESULTS['fatal_error']}")
        traceback.print_exc()
        exit_code = 1
    finally:
        _hr("STAGE 5 — serializing results")
        results_path = out_dir / "phase2_results.json"
        write_json_atomic(results_path, RESULTS)
        _log("stage5", f"results serialized to {results_path}")
        if args.auto_shutdown:
            if RESULTS.get("all_stages_ok"):
                RESULTS["shutdown"] = schedule_shutdown(args.shutdown_delay)
                write_json_atomic(results_path, RESULTS)
            else:
                RESULTS["shutdown"] = {"scheduled": False,
                                       "reason": "run incomplete — power-off suppressed"}
                write_json_atomic(results_path, RESULTS)
                print("[shutdown] --auto_shutdown set but run incomplete — suppressed.")
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
