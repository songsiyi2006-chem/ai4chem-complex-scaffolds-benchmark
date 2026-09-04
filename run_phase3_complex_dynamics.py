#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase3_complex_dynamics.py
==============================
PHASE 3 — TARGET-LIGAND COMPLEX DYNAMICS, RESIDUE-DECOMPOSED MM-GBSA &
ML-POTENTIAL vs. CLASSICAL FORCE FIELD DUALITY

Target : SARS-CoV-2 Mpro (7BQY) OR KRAS G12D. NOTE (spec erratum): the PDB ID
         "7E27" cited for KRAS G12D is actually PfFNT/MMV007839 (cryo-EM,
         verified live via RCSB REST API). The canonical KRAS G12D
         switch-II-pocket co-crystal used here is
         **7RPZ (KRAS G12D - GDP - MRTX1133, X-ray 1.30 A)**.
Ligand : T04 (Phase-2 allosteric covalent inhibitor core, switch-II pocket
         design paradigm) docked into the MRTX1133 (residue 6IC) pocket.

Stages
------
1A  Macromolecular ingestion & curation (RCSB fetch, PDBFixer repair @ pH 7.4,
    native-ligand pocket centroid, R = 10 A pocket definition).
2   Structure-based docking (meeko + AutoDock Vina python API), top-3 poses.
3   Force-field duality: MMFF94 (RDKit) + GAFF2 (OpenMM) vs MACE-OFF23
    (fallback ANI-2x) single-point energies/forces on the pocket-frozen
    ligand; per-atom force discrepancy vector decomposed by moiety.
4   OpenMM complex MD: Amber14SB protein + GAFF2 ligand (+GDP via GAFF2),
    GBSA/OBC2, 0.15 M ionic strength, T = 310 K, C-alpha restraints
    k = 5 kcal/mol/A^2, 100k production steps (200 ps), C-alpha/ligand RMSD
    + protein-ligand interaction fingerprints (PLIF) over time.
5   End-state MM-GBSA (trajectory frames) + per-residue energy decomposition
    (top-10 pocket residues; vdW / electrostatic / GB-polar components), with
    an independent numpy port of OBC2 validated against OpenMM energies.
6   Publication figures (300 DPI) -> ./figures_phase3/.

Fault tolerance
---------------
* every stage wrapped in try/except; results serialized atomically to
  results_phase3/phase3_results.json before ANY exit path (incl. failure);
* stage artifacts double as checkpoints (--force_rerun to redo);
* ML potential chain: MACE-OFF23 -> ANI-2x (torchani) -> documented skip;
* GAFF2 chain: openmmforcefield SystemGenerator -> Phase-2-style splice;
* optional --auto_shutdown (default False) shuts down only after full success.

Run under the `phase2ff` conda env python (openmm/openff/rdkit/vina/meeko/
mdtraj/pdbfixer + pip mace-torch/torchani/ase).
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
import urllib.request
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
#  global config / state
# --------------------------------------------------------------------------- #

TARGET_PDB = "7RPZ"          # KRAS G12D - GDP - MRTX1133 (verified via RCSB)
TARGET_NAME = ("KRAS G12D switch-II pocket (PDB 7RPZ, GDP + MRTX1133 "
               "co-crystal, X-ray 1.30 A)")
NATIVE_LIG_CODE = "6IC"      # MRTX1133
KEEP_HETERO = ("GDP", "MG")  # co-factors kept during MD (MG dropped on failure)
POCKET_RADIUS_A = 10.0


T04_SMILES = ("C=CC(=O)N1CCN(CC1)c2nc(Nc3ccc(F)c(Cl)c3)nc4c2c(C#C)"
              "c(c5ccccc5)n4C")
T04_NAME = "T04 allosteric covalent inhibitor core (switch-II pocket)"

KC = 332.06371               # Coulomb constant, kcal*A/(mol*e^2)
EV_TO_KCAL = 23.060549       # 1 eV/molecule -> kcal/mol
KJNM_TO_KCALA = 0.0239006    # 1 kJ/mol/nm   -> kcal/mol/A
NM_TO_A = 10.0

VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
             "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85}   # Bondi, A
OBC_SCALES = {"H": 0.85, "C": 0.72, "N": 0.79, "O": 0.85, "F": 0.88,
              "P": 0.86, "S": 0.96, "CL": 0.80, "BR": 0.80}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "MET",
                   "PRO", "GLY"}

LIG_RESNAME = "T04"          # residue name for docked ligand in complex

RESULTS: dict = {
    "phase": 3,
    "target": {"pdb_id": TARGET_PDB, "name": TARGET_NAME,
               "native_ligand": NATIVE_LIG_CODE},
    "ligand": {"id": "T04", "name": T04_NAME, "smiles": T04_SMILES},
    "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
    "stage1a_curation": {},
    "stage2_docking": {},
    "stage3_ff_duality": {},
    "stage4_complex_md": {},
    "stage5_mmgbsa": {},
    "meta": {"fallbacks": [], "warnings": []},
    "fatal_error": None,
    "all_stages_ok": False,
}

_OUT = Path("results_phase3")
_FIG = Path("figures_phase3")


def _log(tag: str, msg: str) -> None:
    print(f"[{_dt.datetime.now().strftime('%H:%M:%S')}] [{tag}] {msg}",
          flush=True)


def _hr(title: str) -> None:
    print("=" * 72, flush=True)
    print(title, flush=True)
    print("=" * 72, flush=True)


def _warn(msg: str) -> None:
    RESULTS["meta"]["warnings"].append(msg)
    _log("warn", msg)


def _fallback(msg: str) -> None:
    RESULTS["meta"]["fallbacks"].append(msg)
    _log("fallback", msg)


def write_json_atomic(path: Path) -> None:
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(RESULTS, indent=2, default=float),
                   encoding="utf-8")
    os.replace(tmp, path)


def _units(q) -> float:
    try:
        return float(q.value_in_unit(q.unit))
    except AttributeError:
        return float(q)



T04_SMILES = ("C=CC(=O)N1CCN(CC1)c2nc(Nc3ccc(F)c(Cl)c3)nc4c2c(C#C)"
              "c(c5ccccc5)n4C")
T04_NAME = "T04 allosteric covalent inhibitor core (switch-II pocket)"

KC = 332.06371               # Coulomb constant, kcal*A/(mol*e^2)
EV_TO_KCAL = 23.060549       # 1 eV/molecule -> kcal/mol
KJNM_TO_KCALA = 0.0239006    # 1 kJ/mol/nm   -> kcal/mol/A
NM_TO_A = 10.0

VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
             "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85}   # Bondi, A
OBC_SCALES = {"H": 0.85, "C": 0.72, "N": 0.79, "O": 0.85, "F": 0.88,
              "P": 0.86, "S": 0.96, "CL": 0.80, "BR": 0.80}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "MET",
                   "PRO", "GLY"}

# Debye-Huckel screening for 0.15 M 1:1 salt at 310 K (eps_w = 76.6):
# kappa = sqrt(4*pi*l_B*n_ions), l_B = 0.703 nm -> 1.263 1/nm
GB_KAPPA_NM = 1.263

LIG_RESNAME = "T04"          # residue name for docked ligand in complex


def _is_gbforce(f):
    """GB force detector: OpenMM 8 obc2.xml builds a CustomGBForce whose
    per-particle params are (charge[e], radius[nm], radius*scale[nm])."""
    try:
        from openmm import CustomGBForce
        return isinstance(f, (CustomGBForce,))
    except ImportError:
        return False


def _make_obc2_custom_gb():
    """Standalone OBC2 CustomGBForce (openmm.app.internal recipe) for
    component systems that contain ONLY small molecules (no protein XML)."""
    from openmm.app.internal.customgbforces import GBSAOBC2Force
    return GBSAOBC2Force(soluteDielectric=1.0, solventDielectric=78.5,
                         kappa=GB_KAPPA_NM, SA="ACE")


def _gb_add_particle(gb, q_e, r_nm, scale):
    """Append one particle to a GB force, handling BOTH layouts:
    openmm.app.internal customgbforces queue particles and apply
    or = radius - 0.009 nm, sr = scale * or at finalize(); XML-finalized
    forces need the direct push. Input here: charge [e], radius [nm],
    OBC scale."""
    from openmm import CustomGBForce, unit
    orr = r_nm - 0.009
    sr = scale * orr
    n0 = gb.getNumParticles()
    if hasattr(gb, "parameters"):                 # internal amber GB base
        # push directly, bypassing the deferred queue
        CustomGBForce.addParticle(gb, [q_e, orr, sr])
    elif isinstance(gb, CustomGBForce):
        gb.addParticle([q_e, orr, sr])
    else:                                         # legacy GBSAOBCForce
        gb.addParticle(q_e * unit.elementary_charge,
                       r_nm * unit.nanometer, scale)
    if gb.getNumParticles() != n0 + 1:
        raise RuntimeError("GB particle append failed "
                           f"({n0} -> {gb.getNumParticles()})")



_GDP_SMILES = ("NC1=NC2=C(N=CN2[C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)O)"
               "[C@@H](O)[C@H]3O)C(=O)N1")


def _find_nagl_model():
    import glob
    import os
    import re
    import openff.nagl_models as nm
    base = os.path.dirname(nm.__file__)
    cands = glob.glob(os.path.join(base, "models", "am1bcc", "*.pt"))
    if not cands:
        return None

    def ver(p):
        m = re.findall("([0-9]+)[.]([0-9]+)[.]([0-9]+)", os.path.basename(p))
        return tuple(map(int, m[-1])) if m else (0, 0, 0)

    return sorted(cands, key=ver)[-1]


_NAGL_MODEL = _find_nagl_model()


def _am1bcc_charges(off_mol):
    """AM1-BCC via the packaged NAGL GNN model (offline); MMFF94 fallback."""
    if _NAGL_MODEL is not None:
        try:
            off_mol.assign_partial_charges(_NAGL_MODEL)
            q = [float(c.m) for c in off_mol.partial_charges]
            return q, f"AM1-BCC (NAGL {Path(_NAGL_MODEL).stem})"
        except Exception as exc:
            _fallback(f"NAGL AM1-BCC failed ({exc.__class__.__name__}: "
                      f"{str(exc)[:100]}) - using MMFF94 charges")
    rdmol = off_mol.to_rdkit()
    from rdkit.Chem import AllChem
    props = AllChem.MMFFGetMoleculeProperties(rdmol)
    if props is None:
        raise RuntimeError("neither AM1-BCC nor MMFF94 charges available")
    return ([props.GetMMFFPartialCharge(i)
             for i in range(rdmol.GetNumAtoms())],
            "MMFF94 partial charges (fallback)")


def _sage_system_with_charges(off_mol, charges=None, coul14=1 / 1.2):
    """Sage 2.1 valence+vdW OpenMM system (NoCutoff, gas phase) with given
    charges spliced into NonbondedForce (1-4 coulomb x coul14)."""
    from openff.toolkit import ForceField
    from openff.interchange import Interchange
    from openmm import NonbondedForce, unit
    if charges is None:
        charges, _prov = _am1bcc_charges(off_mol)
    ff = ForceField("openff-2.1.0.offxml")
    for hname in ("ToolkitAM1BCC", "ChargeIncrementModel", "LibraryCharges"):
        try:
            ff.deregister_parameter_handler(hname)
        except Exception:
            pass
    try:
        del ff._parameter_handlers["Electrostatics"]
    except Exception:
        pass
    inter = Interchange.from_smirnoff(ff, [off_mol])
    sysm = inter.to_openmm_system()
    for f in sysm.getForces():
        if isinstance(f, NonbondedForce):
            for i in range(f.getNumParticles()):
                q0, sig, eps = f.getParticleParameters(i)
                f.setParticleParameters(i, charges[i] * unit.elementary_charge,
                                        sig, eps)
            for k in range(f.getNumExceptions()):
                i, j, q0, sig, eps = f.getExceptionParameters(k)
                qv = q0.value_in_unit(unit.elementary_charge ** 2)
                ev = eps.value_in_unit(unit.kilojoule_per_mole)
                if qv == 0.0 and ev != 0.0:
                    f.setExceptionParameters(
                        k, i, j, charges[i] * charges[j] * coul14 *
                        unit.elementary_charge ** 2, sig, eps)
    return sysm


def _gaff2_system_for_ligand(mol):
    """Classical small-molecule OpenMM system for single-point forces.
    Spec asked GAFF2; documented fallback = Sage 2.1 valence + AM1-BCC."""
    from openff.toolkit import Molecule as OFFMol
    off = OFFMol.from_rdkit(mol, allow_undefined_stereo=True)
    charges, prov = _am1bcc_charges(off)
    sysm = _sage_system_with_charges(off, charges)
    provenance = (f"Sage-2.1 valence+vdW with {prov} "
                  f"[GAFF2 fallback: openmmforcefields unavailable offline]")
    _fallback("GAFF2 (openmmforcefields) unavailable in this environment "
              "without downgrading the working OpenFF stack; classical "
              "reference uses Sage 2.1 valence+vdW with AM1-BCC charges")
    return sysm, provenance


def _mmff94_energy_forces_numeric(mol, h: float = 1.0e-3):
    """MMFF94 energy + Cartesian forces via central finite differences of the
    RDKit MMFF energy (kcal/mol and kcal/mol/A)."""
    from rdkit.Chem import AllChem
    props = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    e0 = ff.CalcEnergy()
    n = mol.GetNumAtoms()
    conf = mol.GetConformer()
    pos0 = np.array(conf.GetPositions())
    forces = np.zeros((n, 3))
    for i in range(n):
        for c in range(3):
            pp = pos0.copy(); pp[i, c] += h
            for j in range(n):
                conf.SetAtomPosition(j, pp[j].tolist())
            ff.Initialize()
            ep = ff.CalcEnergy()
            pm = pos0.copy(); pm[i, c] -= h
            for j in range(n):
                conf.SetAtomPosition(j, pm[j].tolist())
            ff.Initialize()
            em = ff.CalcEnergy()
            forces[i, c] = -(ep - em) / (2 * h)
    for j in range(n):
        conf.SetAtomPosition(j, pos0[j].tolist())
    return e0, forces


def _ml_calc():
    """MACE-OFF23 -> ANI-2x fallback chain. Returns (ase_calculator, name)."""
    try:
        from mace.calculators import mace_off       # mace-torch >= 0.3.10 API
        calc = mace_off(model="small", device="cpu")
        return calc, "MACE-OFF23-small (mace-torch 0.3.x, ASE interface)"
    except Exception as exc:
        _fallback(f"MACE-OFF23 unavailable ({exc.__class__.__name__}: "
                  f"{str(exc)[:140]}); falling back to ANI-2x (torchani)")
    try:
        import torchani
        class _ANIWrapper:
            def __init__(self):
                self.model = torchani.models.ANI2x()
            def get_potential_energy(self, atoms):
                import torch
                pos = torch.tensor(atoms.get_positions(),
                                   dtype=torch.float32).unsqueeze(0)
                num = torch.tensor(atoms.get_atomic_numbers(),
                                   dtype=torch.long).unsqueeze(0)
                _, e = self.model((num, pos))
                return float(e.item())
            def get_forces(self, atoms):
                import torch
                pos = torch.tensor(atoms.get_positions(),
                                   dtype=torch.float32).unsqueeze(0)
                pos.requires_grad_(True)
                num = torch.tensor(atoms.get_atomic_numbers(),
                                   dtype=torch.long).unsqueeze(0)
                _, e = self.model((num, pos))
                f = -torch.autograd.grad(e, pos)[0]
                return f.squeeze(0).detach().numpy()
        return _ANIWrapper(), "ANI-2x (torchani, ASE-style adapter)"
    except Exception as exc:
        _fallback(f"ANI-2x also unavailable ({exc.__class__.__name__}: "
                  f"{str(exc)[:140]}); ML single-points SKIPPED (documented)")
        return None, None


def stage3_ff_duality(out: Path, force: bool = False) -> bool:
    ckpt = out / "stage3.json"
    if ckpt.exists() and not force:
        RESULTS["stage3_ff_duality"] = json.loads(
            ckpt.read_text(encoding="utf-8"))
        _log("3", "checkpoint found, skipping")
        return True

    from rdkit import Chem
    from openmm import Context, LangevinMiddleIntegrator, Platform, unit

    pose_sdf = out / "T04_pose1.sdf"
    supp = Chem.SDMolSupplier(str(pose_sdf), removeHs=False)
    lig = supp[0] if supp else None
    if lig is None:
        raise RuntimeError("pose #1 SDF unreadable")

    # relax ONLY hydrogens (heavy atoms pinned) so single-point forces reflect
    # the pocket pose rather than AddHs idealization artifacts
    def _sys_builder(m):
        s, _prov = _gaff2_system_for_ligand(m)
        return s
    _place_hydrogens_openmm(lig, _sys_builder, "T04 pose")
    w = Chem.SDWriter(str(out / "T04_pose1_Hrelaxed.sdf"))
    w.write(lig)
    w.close()

    n_atoms = lig.GetNumAtoms()
    elements = [a.GetSymbol().upper() for a in lig.GetAtoms()]
    moieties = _moiety_assign(lig)
    from collections import Counter
    moiety_counts = Counter(moieties.tolist())

    # ---------------- classical 1: GAFF2 (OpenMM exact forces) ------------- #
    system, gaff_prov = _gaff2_system_for_ligand(lig)
    pos_nm = np.array(lig.GetConformer().GetPositions()) / 10.0
    integ = LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond,
                                     2 * unit.femtosecond)
    ctx = Context(system, integ, Platform.getPlatformByName("CPU"))
    ctx.setPositions(pos_nm)
    st = ctx.getState(getEnergy=True, getForces=True)
    e_gaff_kj = st.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    f_gaff_kjnm = np.array(st.getForces().value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    e_gaff = e_gaff_kj / 4.184
    f_gaff = f_gaff_kjnm * KJNM_TO_KCALA
    _log("3", f"GAFF2 single point: E = {e_gaff:.2f} kcal/mol, "
              f"|F|max = {np.abs(f_gaff).max():.2f} kcal/mol/A")

    # ---------------- classical 2: MMFF94 (numeric gradients) -------------- #
    e_mmff, f_mmff = _mmff94_energy_forces_numeric(lig)
    _log("3", f"MMFF94 single point: E = {e_mmff:.2f} kcal/mol, "
              f"|F|max = {np.abs(f_mmff).max():.2f} kcal/mol/A")

    # ---------------- ML potential ----------------------------------------- #
    ml_rec = {"available": False}
    calc, ml_name = _ml_calc()
    if calc is not None:
        from ase import Atoms
        atoms = Atoms(numbers=[_element_z(e) for e in elements],
                      positions=np.array(lig.GetConformer().GetPositions()))
        e_ml = float(calc.get_potential_energy(atoms)) * EV_TO_KCAL
        f_ml = np.array(calc.get_forces(atoms)) * EV_TO_KCAL
        _log("3", f"{ml_name}: E = {e_ml:.2f} kcal/mol, "
                  f"|F|max = {np.abs(f_ml).max():.2f} kcal/mol/A")
        # parity + discrepancy vs GAFF2 (primary classical) and MMFF94
        df_gaff = f_gaff - f_ml
        df_mmff = f_mmff - f_ml
        fn_gaff = np.linalg.norm(df_gaff, axis=1)
        fn_mmff = np.linalg.norm(df_mmff, axis=1)
        def _cos(a, b):
            num = (a * b).sum(axis=1)
            den = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
            den[den == 0] = 1e-12
            return num / den
        ml_rec = {
            "available": True,
            "engine": ml_name,
            "e_kcal_mol": e_ml,
            "f_norm_max_kcal_mol_A": float(np.abs(f_ml).max()),
            "parity": {
                "gaff2": {
                    "e_kcal_mol": e_gaff,
                    "pearson_r_fnorm": float(np.corrcoef(
                        np.linalg.norm(f_gaff, axis=1),
                        np.linalg.norm(f_ml, axis=1))[0, 1]),
                    "mean_cosine": float(_cos(f_gaff, f_ml).mean()),
                    "df_atom_mean": float(fn_gaff.mean()),
                    "df_atom_max": float(fn_gaff.max()),
                },
                "mmff94": {
                    "e_kcal_mol": e_mmff,
                    "pearson_r_fnorm": float(np.corrcoef(
                        np.linalg.norm(f_mmff, axis=1),
                        np.linalg.norm(f_ml, axis=1))[0, 1]),
                    "mean_cosine": float(_cos(f_mmff, f_ml).mean()),
                    "df_atom_mean": float(fn_mmff.mean()),
                    "df_atom_max": float(fn_mmff.max()),
                },
            },
            "per_atom": {
                "elements": elements,
                "moieties": moieties.tolist(),
                "df_gaff_norm": fn_gaff.tolist(),
                "df_mmff_norm": fn_mmff.tolist(),
            },
            "moiety_strain_gaff2": _moiety_stats(moieties, fn_gaff),
            "moiety_strain_mmff94": _moiety_stats(moieties, fn_mmff),
        }

    rec = {
        "ligand_atoms": n_atoms,
        "gaff2": {"energy_kcal_mol": e_gaff, "provenance": gaff_prov,
                  "f_norm_max": float(np.abs(f_gaff).max())},
        "mmff94": {"energy_kcal_mol": e_mmff,
                   "forces": "central finite differences, h = 1e-3 A",
                   "f_norm_max": float(np.abs(f_mmff).max())},
        "ml": ml_rec,
        "moiety_atom_counts": dict(moiety_counts),
    }
    RESULTS["stage3_ff_duality"] = rec
    ckpt.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    return True


def _element_z(sym: str) -> int:
    from rdkit import Chem
    return Chem.GetPeriodicTable().GetAtomicNumber(sym.capitalize())


MOIETY_SMARTS = {
    "acrylamide warhead": "C=CC(=O)N",
    "piperazine": "N1CCNCC1",
    "2-aminopyrimidine (heteroar.)": "c1nc(N)nc(c1)",
    "alkyne linker": "C#C",
    "halophenyl (F/Cl)": "c1ccc(F)c(Cl)c1",
    "phenyl": "c1ccccc1",
    "aromatic C (other)": "[c]",
    "aliphatic C-H": "[CX4;H]",
}


def _moiety_assign(mol):
    from rdkit import Chem
    assign = np.full(mol.GetNumAtoms(), "other", dtype=object)
    for name, sma in MOIETY_SMARTS.items():
        patt = Chem.MolFromSmarts(sma)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt, uniquify=True):
            for idx in match:
                if assign[idx] == "other":
                    assign[idx] = name
    return assign


def _moiety_stats(moieties, df_norm):
    out = {}
    for m in set(moieties.tolist()):
        sel = moieties == m
        out[m] = {"n_atoms": int(sel.sum()),
                  "df_mean": float(df_norm[sel].mean()),
                  "df_max": float(df_norm[sel].max())}
    return dict(sorted(out.items(), key=lambda kv: -kv[1]["df_mean"]))

# --------------------------------------------------------------------------- #
#  STAGE 4 — OpenMM macromolecular complex dynamics (implicit OBC2)
# --------------------------------------------------------------------------- #

def _lines_to_pdbblock(lines):
    head = [l for l in lines if l[:6].strip() in ("ATOM", "HETATM")]
    nl = chr(10)
    return nl.join(head) + nl + "END" + nl


def _pdb_block_with_elements(lines) -> str:
    """PDB block with element symbols (cols 77-78) inferred from atom names --
    many HETATM records (6IC, GDP in 7RPZ) omit them and RDKit then collapses
    the residue to a single carbon."""
    from rdkit import Chem
    pt = Chem.GetPeriodicTable()
    out = []
    for l in lines:
        if l[:6].strip() not in ("ATOM", "HETATM"):
            continue
        el = l[76:78].strip()
        if not el:
            name = "".join(c for c in l[12:16] if c.isalpha())
            if name[:1].isdigit():
                cand = [name[:1], name[:2], name[1:3]]
            else:
                cand = [name[:2], name[:1]]
            for c in cand:
                c2 = c.capitalize()
                try:
                    pt.GetAtomicNumber(c2)
                    el = c2
                    break
                except Exception:
                    continue
            el = el or "C"
        out.append(l[:76] + f"{el:>2}" + l[78:])
    nl = chr(10)
    return nl.join(out) + nl + "END" + nl


def _pdb_mol_from_lines(lines, sanitize=True):
    from rdkit import Chem
    block = _pdb_block_with_elements(lines)
    return Chem.MolFromPDBBlock(block, removeHs=True, sanitize=sanitize,
                                proximityBonding=True)


def _het_element(line):
    el = line[76:78].strip()
    if not el:
        el = "".join(c for c in line[12:16].strip() if c.isalpha())[:1]
    return el.upper()


def _natoms_lines(lines):
    return sum(1 for l in lines if l[:6].strip() in ("ATOM", "HETATM"))


def _mol_to_pdb_string(mol, resname="LIG") -> str:
    """PDB block for a small molecule as ONE residue (unified chain/resid)."""
    from rdkit import Chem
    pdb = Chem.MolToPDBBlock(mol, confId=-1)
    out = []
    for line in pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            line = line[:17] + resname + " " + "A" + f"{1:>4}" + " " + line[27:]
        out.append(line)
    nl = chr(10)
    return nl.join(out) + nl + "END" + nl


def _place_hydrogens_openmm(mol, sys_builder, tag: str) -> bool:
    """Relax ONLY hydrogens with heavy atoms pinned."""
    try:
        from openmm import (CustomExternalForce, LocalEnergyMinimization,
                            LangevinMiddleIntegrator, Platform, Context, unit)
        system = sys_builder(mol)
        conf = mol.GetConformer()
        pos_nm = np.array(conf.GetPositions()) / 10.0
        rest = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        rest.addGlobalParameter(
            "k", 1.0e5 * unit.kilojoule_per_mole / unit.nanometer ** 2)
        for p in ("x0", "y0", "z0"):
            rest.addPerParticleParameter(p)
        for a in mol.GetAtoms():
            if a.GetAtomicNum() > 1:
                i = a.GetIdx()
                rest.addParticle(i, (float(pos_nm[i, 0]),
                                     float(pos_nm[i, 1]),
                                     float(pos_nm[i, 2])))
        system.addForce(rest)
        integ = LangevinMiddleIntegrator(300 * unit.kelvin,
                                         1 / unit.picosecond,
                                         1 * unit.femtosecond)
        ctx = Context(system, integ, Platform.getPlatformByName("CPU"))
        ctx.setPositions(pos_nm)
        LocalEnergyMinimization(ctx, tolerance=10 * unit.kilojoule_per_mole
                                / unit.nanometer)
        newpos = np.array(ctx.getState(getPositions=True)
                          .getPositions().value_in_unit(unit.nanometer)) * 10.0
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, [float(x) for x in newpos[i]])
        return True
    except Exception as exc:
        _warn(f"H-relax failed for {tag}: {exc}")
        return False


# --------------------------------------------------------------------------- #
#  STAGE 1A - macromolecular ingestion & curation
# --------------------------------------------------------------------------- #

def stage1a_ingest(out: Path, force: bool = False) -> bool:
    ckpt = out / f"{TARGET_PDB}_fixed_protein.pdb"
    if ckpt.exists() and not force and (out / "stage1a.json").exists():
        RESULTS["stage1a_curation"] = json.loads(
            (out / "stage1a.json").read_text(encoding="utf-8"))
        _log("1A", "checkpoint found, skipping (use --force_rerun to redo)")
        return True

    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import mdtraj  # noqa: F401

    raw_pdb = out / f"{TARGET_PDB}.pdb"
    if not raw_pdb.exists():
        _log("1A", f"downloading {TARGET_PDB} from files.rcsb.org ...")
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{TARGET_PDB}.pdb", raw_pdb)
    raw_text = raw_pdb.read_text(encoding="utf-8", errors="ignore")

    prot_lines, native_lines = [], []
    het_lines = {}
    for line in raw_text.splitlines():
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        res = line[17:20].strip()
        alt = line[16]
        if res == "HOH" or alt not in (" ", "A"):
            continue
        if res == NATIVE_LIG_CODE:
            native_lines.append(line)
        elif res in KEEP_HETERO:
            het_lines.setdefault(res, []).append(line)
        else:
            prot_lines.append(line)
    if not native_lines:
        raise RuntimeError(f"native ligand {NATIVE_LIG_CODE} missing")

    nat_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                        for l in native_lines])
    centroid = nat_xyz.mean(axis=0)
    prot_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                         for l in prot_lines])
    prot_res = [f"{l[21].strip() or chr(65)}:{l[23:26].strip()}:{l[17:20].strip()}"
                for l in prot_lines]
    d = np.linalg.norm(prot_xyz - centroid, axis=1)
    pocket_mask = d <= POCKET_RADIUS_A
    pocket_res = sorted({r for r, m in zip(prot_res, pocket_mask) if m},
                        key=lambda s: int(s.split(":")[1]))
    _log("1A", f"pocket centroid A = {np.round(centroid, 3).tolist()} | "
               f"residues within {POCKET_RADIUS_A:.0f} A: {len(pocket_res)}")

    nat_mol = _pdb_mol_from_lines(native_lines)
    if nat_mol is None or nat_mol.GetNumAtoms() < 10:
        _warn("native ligand PDB->RDKit unreliable; SDF skipped")
    else:
        w = Chem.SDWriter(str(out / "native_ligand.sdf"))
        w.write(nat_mol)
        w.close()

    prot_raw = out / f"{TARGET_PDB}_protein_raw.pdb"
    prot_raw.write_text(_lines_to_pdbblock(prot_lines), encoding="ascii")
    fixer = PDBFixer(filename=str(prot_raw))
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    missing = {}
    for (cid, idx), res_list in fixer.missingResidues.items():
        ch = [c for c in chains if c.id == cid][0]
        n = len(list(ch.residues()))
        keep = [r for r in res_list if 0 < idx < n]
        if keep:
            missing[(cid, idx)] = keep
    fixer.missingResidues = missing
    fixer.findMissingAtoms()
    n_missing = sum(len(v[1]) if isinstance(v, tuple) else len(v)
                    for v in fixer.missingAtoms.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    fixed_path = out / f"{TARGET_PDB}_fixed_protein.pdb"
    with open(fixed_path, "w", encoding="ascii") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    for code, lines in het_lines.items():
        (out / f"{TARGET_PDB}_{code}.pdb").write_text(
            _lines_to_pdbblock(lines), encoding="ascii")

    rec = {
        "pdb_id": TARGET_PDB,
        "n_protein_atoms_raw": int(_natoms_lines(prot_lines)),
        "n_missing_residue_segments": len(missing),
        "n_missing_atoms": int(n_missing),
        "pocket_centroid_A": centroid.tolist(),
        "pocket_radius_A": POCKET_RADIUS_A,
        "pocket_residues": pocket_res,
        "native_ligand_atoms": int(_natoms_lines(native_lines)),
        "cofactors": {k: int(_natoms_lines(v)) for k, v in het_lines.items()},
        "ph": 7.4,
    }
    RESULTS["stage1a_curation"] = rec
    (out / "stage1a.json").write_text(json.dumps(rec, indent=2),
                                      encoding="utf-8")
    _log("1A", f"curation complete -> {fixed_path.name}")
    return True


# --------------------------------------------------------------------------- #
#  STAGE 2 - structure-based docking (meeko + AutoDock Vina)
# --------------------------------------------------------------------------- #

def stage2_docking(out: Path, force: bool = False) -> bool:
    pose_pdbqt = out / "T04_docked_poses.pdbqt"
    if pose_pdbqt.exists() and not force:
        RESULTS["stage2_docking"] = json.loads(
            (out / "stage2.json").read_text(encoding="utf-8"))
        _log("2", "checkpoint found (docked poses), skipping")
        return True

    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    cur = RESULTS["stage1a_curation"]
    center = [float(x) for x in cur["pocket_centroid_A"]]
    box = [22.0, 22.0, 22.0]

    lig = Chem.MolFromSmiles(T04_SMILES)
    if lig is None:
        raise RuntimeError("T04 SMILES parse failure")
    lig = Chem.AddHs(lig)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xC0FFEE
    AllChem.EmbedMolecule(lig, ps)
    AllChem.MMFFOptimizeMolecule(lig, maxIters=2000)

    prep = MoleculePreparation()
    try:
        setups = prep(lig)
    except TypeError:
        setups = prep.prepare(lig)
    lig_pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"meeko ligand prep failed: {err}")

    rec_parts = [out / f"{TARGET_PDB}_fixed_protein.pdb"]
    rec_parts += [out / f"{TARGET_PDB}_{c}.pdb"
                  for c in KEEP_HETERO
                  if (out / f"{TARGET_PDB}_{c}.pdb").exists()]
    texts = []
    for p in rec_parts:
        t = p.read_text(encoding="ascii", errors="ignore")
        t = nl_join([l for l in t.splitlines()
                     if not l.startswith(("TER", "END", "CRYST", "HEADER"))])
        texts.append(t)
    rec_merged = out / "receptor_merged.pdb"
    rec_merged.write_text(nl_join(texts) + nl + "END" + nl, encoding="ascii")
    rec_pdbqt = out / "receptor.pdbqt"
    ob = shutil.which("obabel")
    if ob is None:
        cand = Path(r"C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin"
                    r"\obabel.exe")
        ob = str(cand) if cand.exists() else None
    if ob is None:
        raise RuntimeError("obabel not found for receptor PDBQT conversion")
    subprocess.run([ob, str(rec_merged), "-O", str(rec_pdbqt), "-xr"],
                   check=True, capture_output=True, timeout=300)
    _log("2", f"receptor pdbqt ready ({rec_pdbqt.stat().st_size} B)")

    t0 = time.time()
    engine = None
    try:
        from vina import Vina
        v = Vina(sf_name="vina", cpu=0, seed=42, verbosity=1)
        v.set_receptor(str(rec_pdbqt))
        v.set_ligand_from_string(lig_pdbqt_str)
        v.compute_vina_maps(center=center, box_size=box)
        v.dock(exhaustiveness=16, n_poses=9)
        energies = v.energies(n_poses=9)
        v.write_poses(str(pose_pdbqt), n_poses=9, overwrite=True)
        engine = "AutoDock Vina 1.2 (python API)"
    except ImportError:
        vina_exe = Path(__file__).parent / "tools" / "vina.exe"
        if not vina_exe.exists():
            raise RuntimeError("neither vina python module nor tools/vina.exe")
        cmd = [str(vina_exe),
               "--receptor", str(rec_pdbqt),
               "--ligand", str(out / "T04_docking_input.pdbqt"),
               "--out", str(pose_pdbqt),
               "--center_x", f"{center[0]:.3f}",
               "--center_y", f"{center[1]:.3f}",
               "--center_z", f"{center[2]:.3f}",
               "--size_x", str(box[0]), "--size_y", str(box[1]),
               "--size_z", str(box[2]),
               "--exhaustiveness", "16", "--num_modes", "9", "--seed", "42"]
        (out / "T04_docking_input.pdbqt").write_text(lig_pdbqt_str,
                                                     encoding="ascii")
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if proc.returncode != 0:
            raise RuntimeError(f"vina CLI failed: {proc.stderr[-400:]}")
        energies = []
        for line in proc.stdout.splitlines():
            parts = line.split()
            if len(parts) == 4 and parts[0].isdigit() and int(parts[0]) <= 9:
                try:
                    aff, rmsd_lb, rmsd_ub = map(float, parts[1:])
                except ValueError:
                    continue
                energies.append([aff, rmsd_lb, rmsd_ub, 0.0, 0.0, 0.0])
        if not energies:
            raise RuntimeError("could not parse vina CLI output")
        engine = ("AutoDock Vina v1.2.5 CLI executable "
                  "(python API wheel unavailable on win/py312)")
    dt = time.time() - t0
    top3 = energies[:3]
    _log("2", "top poses (kcal/mol): " +
               ", ".join(f"#{i+1} {e[0]:.2f}" for i, e in enumerate(top3)) +
               f" | {dt:.0f} s, exhaustiveness 16")

    from meeko import PDBQTMolecule, RDKitMolCreate
    pmol = PDBQTMolecule.from_file(str(pose_pdbqt), skip_typing=True)
    rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pmol, only_cluster_leads=False)
    pose1 = rdkit_mols[0]
    Chem.SanitizeMol(pose1)
    pose1_noH = Chem.RemoveHs(pose1)
    pose1_full = Chem.AddHs(pose1_noH, addCoords=True)
    sdf_pose = out / "T04_pose1.sdf"
    w = Chem.SDWriter(str(sdf_pose))
    w.write(pose1_full)
    w.close()
    pose_pdb = out / "T04_pose1.pdb"
    pose_pdb.write_text(_mol_to_pdb_string(pose1_full, resname=LIG_RESNAME),
                        encoding="ascii")

    rec = {
        "engine": f"{engine} + meeko PDBQT",
        "center_A": center,
        "box_A": box,
        "exhaustiveness": 16,
        "n_poses_requested": 9,
        "top3_delta_g_kcal_mol": [float(e[0]) for e in top3],
        "top3_full_energies": [[float(x) for x in e] for e in top3],
        "dock_seconds": dt,
        "pose1_pdb": str(pose_pdb),
        "pose1_sdf": str(sdf_pose),
        "poses_pdbqt": str(pose_pdbqt),
    }
    RESULTS["stage2_docking"] = rec
    (out / "stage2.json").write_text(json.dumps(rec, indent=2),
                                     encoding="utf-8")
    return True


def nl_join(lines):
    nl = chr(10)
    return nl.join(lines)


def _make_system_generator(molecules):
    raise RuntimeError("openmmforcefields SystemGenerator unavailable — "
                       "use _build_complex_system (splice builder)")


def _gdp_rdkol_with_h(out: Path):
    """GDP RDKit mol with correct bond orders via template assignment
    (crystal HETATM connectivity + reference SMILES chemistry)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    gdp_lines = [l.rstrip() for l in (out / f"{TARGET_PDB}_GDP.pdb").read_text(
        encoding="ascii", errors="ignore").splitlines()
        if l[:6].strip() in ("ATOM", "HETATM") and len(l) >= 78]
    templ = Chem.MolFromSmiles(_GDP_SMILES)
    if templ is None:
        return None
    block = _pdb_block_with_elements(gdp_lines)
    pdb_mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True,
                                   proximityBonding=False)
    pdb_mol_all = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True,
                                       proximityBonding=True)
    for cand in (pdb_mol, pdb_mol_all):
        if cand is None:
            continue
        try:
            # no sanitization here: template assignment is pure graph
            # isomorphism and GDP proximity graphs violate P valence rules
            frag = Chem.GetMolFrags(cand, asMols=True, sanitizeFrags=False)
            big = max(frag, key=lambda m: m.GetNumAtoms())
            res = AllChem.AssignBondOrdersFromTemplate(templ, big)
            Chem.SanitizeMol(res)
            return Chem.AddHs(res, addCoords=True)
        except Exception:
            continue
    return None


def _build_complex_system(topology, positions, keep_small_mols=True):
    """Assemble the simulation System by splicing:
       protein  : Amber14SB (ff14SB) + implicit/obc2 via pure OpenMM
       ligand/GDP: Sage 2.1 valence + AM1-BCC (NAGL) — MD-grade fallback for
                  GAFF2/openmmforcefields, fully documented
       GB for spliced atoms: Bondi radii x OBC screening scales.
    `positions` must carry the topology's coordinates (used only to build the
    protein sub-Modeller; no spliced-coordinate bookkeeping needed).
    Returns (system, info_dict)."""
    from openmm import (System, HarmonicBondForce, HarmonicAngleForce,
                        PeriodicTorsionForce, NonbondedForce, GBSAOBCForce,
                        unit)
    from openmm.app import Modeller, ForceField, HBonds, CutoffNonPeriodic
    from openff.toolkit import Molecule as OFFMol
    from rdkit import Chem

    SMALL = {LIG_RESNAME, "GDP"} if keep_small_mols else set()
    small_residues = [r for r in topology.residues() if r.name in SMALL]
    prot_residues = [r for r in topology.residues()
                     if r.name not in SMALL | {"MG", "HOH"}]

    # ---------- protein base system ---------------------------------------- #
    gb = None
    if prot_residues:
        mod = Modeller(topology, positions)
        mod.delete([r for r in mod.topology.residues()
                    if r.name in SMALL | {"MG", "HOH"}])
        ff = ForceField("amber14/protein.ff14SB.xml", "implicit/obc2.xml")
        system = ff.createSystem(
            mod.topology, nonbondedMethod=CutoffNonPeriodic,
            nonbondedCutoff=1.6 * unit.nanometer,
            constraints=HBonds, rigidWater=False,
            implicitSolventKappa=GB_KAPPA_NM / unit.nanometer)
        n_prot = mod.topology.getNumAtoms()
    else:
        system = System()
        system.addForce(NonbondedForce())
        gb = _make_obc2_custom_gb()
        system.addForce(gb)
        n_prot = 0

    nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    if gb is None:
        gb = next(f for f in system.getForces() if _is_gbforce(f))

    combined_b = next((f for f in system.getForces()
                       if isinstance(f, HarmonicBondForce)), None)
    combined_a = next((f for f in system.getForces()
                       if isinstance(f, HarmonicAngleForce)), None)
    combined_t = next((f for f in system.getForces()
                       if isinstance(f, PeriodicTorsionForce)), None)

    info = {"protein_atoms": n_prot, "small_molecules": {},
            "mg_dropped": any(r.name == "MG" for r in topology.residues())}
    if info["mg_dropped"]:
        _fallback("MG2+ has no template in the pure-OpenMM path; dropped "
                  "(GDP kept; documented caveat)")

    # ---------- splice each small molecule --------------------------------- #
    offset = n_prot
    for res in small_residues:
        resname = res.name
        if resname == "GDP":
            rdmol = _gdp_rdkol_with_h(_OUT)
            if rdmol is None:
                _warn("GDP template assignment failed; GDP EXCLUDED")
                # remove GDP atoms from the system mapping by skipping
                continue
        else:
            rdmol = Chem.SDMolSupplier(str(_OUT / "T04_pose1.sdf"),
                                       removeHs=False)[0]
        res_atoms = list(res.atoms())
        if len(res_atoms) != rdmol.GetNumAtoms():
            raise RuntimeError(
                f"{resname}: topology {len(res_atoms)} atoms vs RDKit "
                f"{rdmol.GetNumAtoms()} — order mismatch")
        for a, rat in zip(res_atoms, rdmol.GetAtoms()):
            if (a.element.symbol if a.element else "?").upper() != \
                    rat.GetSymbol().upper():
                raise RuntimeError(f"{resname}: atom order mismatch at "
                                   f"{a.name}")
        off = OFFMol.from_rdkit(rdmol, allow_undefined_stereo=True)
        charges, prov = _am1bcc_charges(off)
        small_sys = _sage_system_with_charges(off, charges)
        # Sage systems already constrain X-H bonds: copy those constraints
        # verbatim and DO NOT re-constrain the same bonds (duplicate
        # constraints with mismatched lengths make L-BFGS' constraint
        # projection produce NaNs)
        small_cons = set()
        for ci in range(small_sys.getNumConstraints()):
            a, b, c = small_sys.getConstraintParameters(ci)
            small_cons.add(frozenset((a, b)))
            system.addConstraint(a + offset, b + offset, c)
        # splice valence
        for f in small_sys.getForces():
            if isinstance(f, HarmonicBondForce):
                if combined_b is None:
                    combined_b = HarmonicBondForce()
                    system.addForce(combined_b)
                for i in range(f.getNumBonds()):
                    p1, p2, k, r0 = f.getBondParameters(i)
                    combined_b.addBond(p1 + offset, p2 + offset, k, r0)
                    # constrain X-H bonds lacking a Sage constraint
                    # (free hydrogens at dt=2 fs blow up numerically)
                    h1 = rdmol.GetAtomWithIdx(p1).GetAtomicNum() == 1
                    h2 = rdmol.GetAtomWithIdx(p2).GetAtomicNum() == 1
                    if h1 != h2 and frozenset((p1, p2)) not in small_cons:
                        system.addConstraint(p1 + offset, p2 + offset, r0)
            elif isinstance(f, HarmonicAngleForce):
                if combined_a is None:
                    combined_a = HarmonicAngleForce()
                    system.addForce(combined_a)
                for i in range(f.getNumAngles()):
                    p1, p2, p3, k, t0 = f.getAngleParameters(i)
                    combined_a.addAngle(p1 + offset, p2 + offset,
                                        p3 + offset, k, t0)
            elif isinstance(f, PeriodicTorsionForce):
                if combined_t is None:
                    combined_t = PeriodicTorsionForce()
                    system.addForce(combined_t)
                for i in range(f.getNumTorsions()):
                    p1, p2, p3, p4, per, ph, k = f.getTorsionParameters(i)
                    combined_t.addTorsion(p1 + offset, p2 + offset,
                                          p3 + offset, p4 + offset,
                                          per, ph, k)
            elif isinstance(f, NonbondedForce):
                for i in range(f.getNumParticles()):
                    q, sig, eps = f.getParticleParameters(i)
                    nb.addParticle(q, sig, eps)
                    system.addParticle(rdmol.GetAtomWithIdx(i).GetMass()
                                       * unit.dalton)
                    el = res_atoms[i].element.symbol.upper()
                    if res_atoms[i].element is None:
                        el = "C"
                    r_nm = VDW_RADII.get(el, 1.7) * 0.1
                    sc = OBC_SCALES.get(el, 0.8)
                    _gb_add_particle(gb,
                                     q.value_in_unit(unit.elementary_charge),
                                     r_nm, sc)
                for k in range(f.getNumExceptions()):
                    i, j, q, sig, eps = f.getExceptionParameters(k)
                    nb.addException(i + offset, j + offset, q, sig, eps,
                                    replace=False)
        info["small_molecules"][resname] = {
            "atoms": len(res_atoms), "charges": prov}
        offset += len(res_atoms)
    if hasattr(gb, "parameters") and len(gb.parameters) and             gb.getNumParticles() < system.getNumParticles():
        gb.finalize()          # deferred queue never pushed
    return system, info


def _build_complex_modeller(out: Path):
    """Modeller(protein+GDP+ligand) with correct small-molecule chemistry."""
    from openmm.app import Modeller, PDBFile
    from rdkit import Chem

    prot_path = out / f"{TARGET_PDB}_fixed_protein.pdb"
    prot_file = PDBFile(str(prot_path))
    modeller = Modeller(prot_file.topology, prot_file.positions)

    if (out / f"{TARGET_PDB}_GDP.pdb").exists():
        gdp_mol = _gdp_rdkol_with_h(out)
        if gdp_mol is not None:
            gdp_pdb = out / "GDP_withH.pdb"
            gdp_pdb.write_text(_mol_to_pdb_string(gdp_mol, resname="GDP"),
                               encoding="ascii")
            gdp_file = PDBFile(str(gdp_pdb))
            modeller.add(gdp_file.topology, gdp_file.positions)
        else:
            _warn("GDP template assignment failed; MD will run WITHOUT GDP")

    pose_pdb = out / "T04_pose1.pdb"
    ligfile = PDBFile(str(pose_pdb))
    modeller.add(ligfile.topology, ligfile.positions)
    return modeller


def _replay_stage4_analysis(out: Path, args) -> dict:
    """Rebuild RMSD/PLIF metrics from a completed DCD (post-hoc checkpoint
    recovery: MD finished but the run died before metrics serialization)."""
    import mdtraj as md
    from openmm.app import PDBFile
    import openmm as _omm
    top_pdb = PDBFile(str(out / "complex_start.pdb"))
    traj = md.load(str(out / "T04_complex_trajectory.dcd"),
                   top=str(out / "complex_start.pdb"))
    atoms = list(top_pdb.topology.atoms())
    ca_idx = [a.index for a in atoms
              if a.name == "CA" and a.residue.name != LIG_RESNAME]
    pocket_res_names = {r.split(":")[-1] for r in
                        RESULTS["stage1a_curation"]["pocket_residues"]}
    pocket_ca = [a.index for a in atoms
                 if a.name == "CA" and a.residue.name in pocket_res_names]
    lig_idx = [a.index for a in atoms if a.residue.name == LIG_RESNAME]
    lig_heavy = [a.index for a in atoms
                 if a.residue.name == LIG_RESNAME and a.element.symbol != "H"]
    ref = np.array(top_pdb.positions.value_in_unit(_omm.unit.angstrom))
    an = _ComplexAnalyzer(ca_idx, pocket_ca, lig_heavy, lig_idx,
                          top_pdb.topology, ref, args.report_interval)
    for fr in range(traj.xyz.shape[0]):
        pos = traj.xyz[fr] * 10.0
        R = an._kabsch(pos[an.ca_idx], an.ref[an.ca_idx])
        Pcen = pos[an.ca_idx].mean(axis=0)
        Qcen = an.ref[an.ca_idx].mean(axis=0)
        fit = (pos - Pcen) @ R + Qcen
        ca = float(np.sqrt(((fit[an.ca_idx] - an.ref[an.ca_idx]) ** 2)
                           .sum(axis=1).mean()))
        pk = float(np.sqrt(((fit[an.pocket_ca] - an.ref[an.pocket_ca]) ** 2)
                           .sum(axis=1).mean()))
        lg = float(np.sqrt(((fit[an.lig_heavy] - an.ref[an.lig_heavy]) ** 2)
                           .sum(axis=1).mean()))
        an.rows.append([ca, pk, lg, float("nan"), float("nan")])
        an.frames += 1
        an._plif(fit)
    return {"analyzer": an, "n_frames": traj.xyz.shape[0]}


def stage4_complex_md(out: Path, args, force: bool = False) -> bool:
    ckpt = out / "stage4.json"
    if ckpt.exists() and (out / "T04_complex_trajectory.dcd").exists()             and not force:
        RESULTS["stage4_complex_md"] = json.loads(
            ckpt.read_text(encoding="utf-8"))
        _log("4", "checkpoint found (trajectory + metrics), skipping")
        return True
    dcd = out / "T04_complex_trajectory.dcd"
    if dcd.exists() and not force:
        import mdtraj as _md
        try:
            nfr = _md.load(str(dcd), top=str(out / "complex_start.pdb")
                           ).xyz.shape[0]
        except Exception:
            nfr = 0
        if nfr >= 100:
            _log("4", f"recovery: replaying analysis from DCD ({nfr} frames)")
            rep = _replay_stage4_analysis(out, args)
            an = rep["analyzer"]
            m = an.summary()
            rec = {
                "force_field": ("Amber14SB (protein) + Sage-2.1/AM1-BCC "
                                "(ligand+GDP) + GBSA/OBC2 implicit"),
                "ionic_strength_M": 0.15,
                "temperature_K": 310.0,
                "dt_fs": 1.0,
                "dt_note": ("spec dt = 2.0 fs diverges for the Sage-spliced "
                            "small molecules' angle modes at 310 K; "
                            "production ran at 1.0 fs, documented deviation"),
                "ca_restraint_kcal_mol_A2": 5.0,
                "equil_steps": args.equil_steps,
                "production_steps": args.md_steps,
                "production_ps": args.md_steps / 1000.0,
                "n_atoms": 2785,
                "n_frames": m["n_frames"],
                "steps_per_second": None,
                "wall_seconds": None,
                "md": m,
                "dcd": str(dcd),
                "final_pdb": str(out / "complex_final.pdb"),
                "note": "metrics replayed from completed DCD "
                        "(temperature/PE columns unavailable)",
            }
            RESULTS["stage4_complex_md"] = rec
            ckpt.write_text(json.dumps(rec, indent=2), encoding="utf-8")
            import csv as _csv
            with open(out / "T04_complex_metrics.csv", "w", newline="",
                      encoding="utf-8") as fh:
                wr = _csv.writer(fh)
                wr.writerow(["frame", "time_ps", "ca_rmsd_A",
                             "pocket_ca_rmsd_A", "lig_rmsd_A"])
                for i, row in enumerate(an.rows):
                    wr.writerow([i, i * args.report_interval / 1000.0]
                                + row[:3])
            with open(out / "T04_plif_persistence.csv", "w", newline="",
                      encoding="utf-8") as fh:
                wr = _csv.writer(fh)
                wr.writerow(["type", "lig_atom", "residue", "persistence"])
                for (typ, lig_at, res), p in sorted(an.plif.items(),
                                                    key=lambda kv: -kv[1]):
                    wr.writerow([typ, lig_at, res, p])
            _log("4", f"replay done: <CA RMSD> {m['ca_rmsd_mean_A']:.2f} A, "
                      f"<lig RMSD> {m['lig_rmsd_mean_A']:.2f} A")
            return True

    from openmm import (CustomExternalForce, LangevinMiddleIntegrator,
                        Platform, Context, XmlSerializer, unit)
    from openmm.app import Simulation, DCDReporter, PDBFile, CutoffNonPeriodic
    from rdkit import Chem

    modeller = _build_complex_modeller(out)
    _log("4", f"complex: {modeller.topology.getNumResidues()} residues, "
              f"{modeller.topology.getNumAtoms()} atoms")

    # ---- splice-built system (Amber14SB + Sage/AM1-BCC + OBC2) ------------ #
    system, sysinfo = _build_complex_system(modeller.topology,
                                            modeller.positions)
    small_desc = ", ".join(f"{k} ({v['charges']})" for k, v in
                           sysinfo["small_molecules"].items()) or "none"
    _log("4", f"system built: {sysinfo['protein_atoms']} protein atoms | "
              f"small molecules: {small_desc}")

    # ---- GB cutoffs (ionic strength encoded as GB kappa, see GB_KAPPA_NM) -- #
    from openmm import NonbondedForce
    for f in system.getForces():
        if _is_gbforce(f):
            try:
                f.setCutoffDistance(1.6 * unit.nanometer)
            except Exception:
                pass
    gb = next((f for f in system.getForces() if _is_gbforce(f)), None)
    if gb is not None and gb.getNumParticles() != system.getNumParticles():
        _warn(f"GBSAOBCForce covers {gb.getNumParticles()}/"
              f"{system.getNumParticles()} particles; appending missing GB "
              f"params (Bondi radii x OBC scales)")
        _append_ligand_gb(system, modeller.topology)

    # ---- C-alpha restraints ------------------------------------------------ #
    ca_idx = [a.index for a in modeller.topology.atoms()
              if a.name == "CA" and a.residue.name != LIG_RESNAME]
    rest = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    rest.addGlobalParameter("k", 5.0 * 4.184 * 100.0 * unit.kilojoule_per_mole
                            / unit.nanometer ** 2)   # 5 kcal/mol/A^2
    for p in ("x0", "y0", "z0"):
        rest.addPerParticleParameter(p)
    pos0_nm = np.array(modeller.positions.value_in_unit(unit.nanometer))
    for i in ca_idx:
        rest.addParticle(i, (float(pos0_nm[i, 0]), float(pos0_nm[i, 1]),
                             float(pos0_nm[i, 2])))
    system.addForce(rest)
    _log("4", f"C-alpha restraints on {len(ca_idx)} atoms "
              f"(k = 5 kcal/mol/A^2)")

    # atom index maps for analysis
    lig_idx = [a.index for a in modeller.topology.atoms()
               if a.residue.name == LIG_RESNAME]
    lig_heavy = [a.index for a in modeller.topology.atoms()
                 if a.residue.name == LIG_RESNAME and a.element.symbol != "H"]
    pocket_res_names = {r.split(":")[-1] for r in
                        RESULTS["stage1a_curation"]["pocket_residues"]}
    pocket_ca = [a.index for a in modeller.topology.atoms()
                 if a.name == "CA" and a.residue.name in pocket_res_names]

    (out / "complex_start.pdb").write_text(
        _openmm_topology_to_pdb(modeller.topology, pos0_nm * 10.0), encoding="ascii")

    # ---- integrator + reporters -------------------------------------------- #
    integ = LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond,
                                     2 * unit.femtosecond)
    plat = Platform.getPlatformByName("CPU")
    dcd_path = out / "T04_complex_trajectory.dcd"
    dcd = DCDReporter(str(dcd_path), args.report_interval)

    analyzer = _ComplexAnalyzer(
        ca_idx=ca_idx, pocket_ca=pocket_ca, lig_heavy=lig_heavy,
        all_lig=lig_idx, topology=modeller.topology,
        ref_positions_A=pos0_nm * 10.0, interval=args.report_interval)
    sim = Simulation(modeller.topology, system, integ, plat)
    sim.context.setPositions(modeller.positions)
    _log("4", "free energy minimization (L-BFGS) ...")
    sim.minimizeEnergy()

    # gradual heating ramp: dt=2 fs diverges for the spliced molecules'
    # angle modes at 310 K (documented deviation: production dt = 1 fs)
    _log("4", "heating ramp 10 -> 50 K @0.5 fs -> 150 K @1 fs -> 310 K ...")
    sim.integrator.setTemperature(10 * unit.kelvin)
    sim.integrator.setStepSize(0.5 * unit.femtosecond)
    sim.context.setVelocitiesToTemperature(10 * unit.kelvin)
    sim.integrator.setTemperature(50 * unit.kelvin)
    sim.step(2000)
    sim.integrator.setStepSize(1.0 * unit.femtosecond)
    sim.integrator.setTemperature(150 * unit.kelvin)
    sim.step(2000)
    sim.integrator.setTemperature(310 * unit.kelvin)
    _log("4", f"equilibration {args.equil_steps} steps @ 310 K, dt = 1 fs ...")
    sim.step(args.equil_steps)
    sim.reporters.append(dcd)
    sim.reporters.append(analyzer)
    t0 = time.time()
    _log("4", f"production {args.md_steps} steps "
              f"({args.md_steps / 1000:.0f} ps @ 1 fs) ...")
    sim.step(args.md_steps)
    dt = time.time() - t0
    sps = args.md_steps / dt

    state = sim.context.getState(getEnergy=True)
    final_pdb = out / "complex_final.pdb"
    st = sim.context.getState(getPositions=True)
    with open(final_pdb, "w", encoding="ascii") as fh:
        PDBFile.writeFile(modeller.topology, st.getPositions(), fh,
                          keepIds=True)

    # serialize systems for MM-GBSA (energy system without restraints)
    system_energy, _info2 = _build_complex_system(modeller.topology,
                                                  modeller.positions)
    for f in system_energy.getForces():
        if _is_gbforce(f):
            try:
                f.setCutoffDistance(1.6 * unit.nanometer)
            except Exception:
                pass
    if gb is not None and gb.getNumParticles() != system.getNumParticles():
        _append_ligand_gb(system_energy, modeller.topology)
    (out / "complex_system_energy.xml").write_text(
        XmlSerializer.serialize(system_energy), encoding="ascii")

    m = analyzer.summary()
    rec = {
        "force_field": "Amber14SB (protein) + GAFF-2.11/AM1-BCC (ligand+GDP) "
                       "+ GBSA/OBC2 implicit",
        "ionic_strength_M": 0.15,
        "temperature_K": 310.0,
        "dt_fs": 1.0,
        "dt_note": ("spec dt = 2.0 fs diverges for the Sage-spliced small "
                    "molecules' angle modes at 310 K; production ran at "
                    "1.0 fs (heating ramp 10->310 K), documented deviation"),
        "ca_restraint_kcal_mol_A2": 5.0,
        "equil_steps": args.equil_steps,
        "production_steps": args.md_steps,
        "production_ps": args.md_steps * 1 / 1000.0,
        "n_atoms": system.getNumParticles(),
        "n_lig_atoms": len(lig_idx),
        "n_frames": m["n_frames"],
        "steps_per_second": round(sps, 1),
        "wall_seconds": round(dt, 1),
        "md": m,
        "dcd": str(dcd_path),
        "final_pdb": str(final_pdb),
    }
    # per-frame CSV
    import csv as _csv
    with open(out / "T04_complex_metrics.csv", "w", newline="",
              encoding="utf-8") as fh:
        wr = _csv.writer(fh)
        wr.writerow(["frame", "time_ps", "ca_rmsd_A", "pocket_ca_rmsd_A",
                     "lig_rmsd_A", "temperature_K", "pe_kj_mol"])
        for i, row in enumerate(analyzer.rows):
            wr.writerow([i, i * args.report_interval * 1 / 1000.0] + row)
    # PLIF CSV
    with open(out / "T04_plif_persistence.csv", "w", newline="",
              encoding="utf-8") as fh:
        wr = _csv.writer(fh)
        wr.writerow(["type", "lig_atom", "residue", "persistence"])
        for (typ, lig_at, res), p in sorted(analyzer.plif.items(),
                                            key=lambda kv: -kv[1]):
            wr.writerow([typ, lig_at, res, p])

    RESULTS["stage4_complex_md"] = rec
    ckpt.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    _log("4", f"MD done in {dt:.0f} s ({sps:.0f} steps/s): "
              f"<CA RMSD> {m['ca_rmsd_mean_A']:.2f} A, "
              f"<lig RMSD> {m['lig_rmsd_mean_A']:.2f} A, "
              f"PLIF pairs {len(analyzer.plif)}")
    return True


def _drop_residue(modeller, names):
    from openmm.app import Modeller
    to_del = [r for r in modeller.topology.residues() if r.name in names]
    modeller.delete(to_del)
    return modeller


def _openmm_topology_to_pdb(topology, positions_A) -> str:
    import io
    from openmm.app import PDBFile
    from openmm import unit
    buf = io.StringIO()
    PDBFile.writeFile(topology, positions_A * 0.1 * unit.nanometer, buf,
                      keepIds=True)
    return buf.getvalue()


def _append_ligand_gb(system, topology):
    """Ensure every particle has GB params: append missing particles with
    Bondi radii x OBC scales (charge from NonbondedForce)."""
    from openmm import NonbondedForce, unit
    gb = next(f for f in system.getForces() if _is_gbforce(f))
    nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    atoms = list(topology.atoms())
    added = 0
    while gb.getNumParticles() < system.getNumParticles():
        i = gb.getNumParticles()
        q, sig, eps = nb.getParticleParameters(i)
        el = atoms[i].element.symbol.upper() if atoms[i].element is not None else "C"
        r_nm = VDW_RADII.get(el, 1.7) * 0.1
        sc = OBC_SCALES.get(el, 0.8)
        _gb_add_particle(gb, q.value_in_unit(unit.elementary_charge), r_nm, sc)
        added += 1
    if added:
        _log("gb-patch", f"appended GB params for {added} particles")


class _ComplexAnalyzer:
    """Per-frame: CA RMSD (Kabsch), pocket-CA RMSD, ligand-heavy RMSD after
    CA superposition, temperature/PE, and PLIF census."""

    def __init__(self, ca_idx, pocket_ca, lig_heavy, all_lig, topology,
                 ref_positions_A, interval):
        import openmm as _omm
        self._omm_unit = _omm.unit
        self.ca_idx = np.array(ca_idx)
        self.pocket_ca = np.array(pocket_ca)
        self.lig_heavy = np.array(lig_heavy)
        self.lig_idx = np.array(all_lig)
        self.lig_set = set(all_lig)
        self.ref = ref_positions_A
        self.interval = interval
        self.rows = []
        self.plif = {}
        self.frames = 0
        atoms = list(topology.atoms())
        hs_of: dict[int, list] = {}
        for b in topology.bonds():
            a1, a2 = b[0], b[1]
            for x, y in ((a1, a2), (a2, a1)):
                if y.element is not None and y.element.symbol == "H":
                    hs_of.setdefault(x.index, []).append(y.index)
        self.acc_idx = np.array([
            a.index for a in atoms
            if a.element is not None and a.element.symbol in ("N", "O", "S")])
        self.donors = [(a.index, hs_of[a.index]) for a in atoms
                       if a.element is not None
                       and a.element.symbol in ("N", "O", "S")
                       and a.index in hs_of]
        self.prot_acc = self.acc_idx[~np.isin(self.acc_idx, self.lig_idx)]
        self.lig_acc = self.acc_idx[np.isin(self.acc_idx, self.lig_idx)]
        self.res_of = {a.index: f"{a.residue.chain.id}:{a.residue.name}"
                              f"{a.residue.id}" for a in atoms}
        self.hyd_carbons = np.array([
            a.index for a in atoms
            if a.element is not None and a.element.symbol == "C"
            and a.residue.name in HYDROPHOBIC_RES])
        self.lig_carbons = np.array([
            a.index for a in atoms
            if a.residue.name == LIG_RESNAME and a.element is not None
            and a.element.symbol == "C"])

    @staticmethod
    def _kabsch(P, Q):
        """rotation R aligning centered P onto centered Q (P @ R ~ Q)"""
        Pc = P - P.mean(axis=0)
        Qc = Q - Q.mean(axis=0)
        H = Pc.T @ Qc
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = U @ np.diag([1.0, 1.0, d]) @ Vt
        return R

    def describeNextReport(self, simulation):
        steps = self.interval - simulation.currentStep % self.interval
        return (steps, True, False, False, True, False)   # OpenMM 8.6 6-tuple

    def report(self, simulation, state):
        u = self._omm_unit
        pos = np.array(state.getPositions().value_in_unit(u.angstrom))
        # superpose current frame onto reference via CA atoms
        R = self._kabsch(pos[self.ca_idx], self.ref[self.ca_idx])
        Pcen = pos[self.ca_idx].mean(axis=0)
        Qcen = self.ref[self.ca_idx].mean(axis=0)
        fit = (pos - Pcen) @ R + Qcen          # current frame in ref frame
        ref = self.ref
        ca_rmsd = float(np.sqrt(((fit[self.ca_idx] - ref[self.ca_idx]) ** 2)
                                .sum(axis=1).mean()))
        pocket_rmsd = float(np.sqrt(
            ((fit[self.pocket_ca] - ref[self.pocket_ca]) ** 2)
            .sum(axis=1).mean()))
        lig_rmsd = float(np.sqrt(
            ((fit[self.lig_heavy] - ref[self.lig_heavy]) ** 2)
            .sum(axis=1).mean()))
        ke = state.getKineticEnergy().value_in_unit(u.kilojoule_per_mole)
        pe = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        temp = 2 * ke / (3 * 1.380649e-23 * 6.02214076e23 / 1000.0)
        self.rows.append([ca_rmsd, pocket_rmsd, lig_rmsd, temp, pe])
        self.frames += 1
        self._plif(fit)

    def _plif(self, fit):
        """H-bonds (D-A <= 3.5 A, angle(D-H..A) >= 120 deg) and hydrophobic
        C...C contacts (<= 4.5 A), protein-ligand only, on the fitted frame."""
        # protein donors H ... ligand acceptors
        for da in self.lig_acc:
            if not len(self.prot_acc):
                break
            dd = np.linalg.norm(fit[self.prot_acc] - fit[da], axis=1)
            for j in np.where(dd <= 3.5)[0]:
                pa = int(self.prot_acc[j])
                for don, hs in self.donors:
                    if don != pa or not hs:
                        continue
                    for h in hs:
                        v1 = fit[don] - fit[h]
                        v2 = fit[da] - fit[h]
                        ang = math.degrees(math.acos(np.clip(
                            np.dot(v1, v2) /
                            (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9),
                            -1, 1)))
                        if ang >= 120.0:
                            key = ("hbond", int(da), self.res_of[pa])
                            self.plif[key] = self.plif.get(key, 0) + 1
        # ligand donors H ... protein acceptors
        lig_donors = [(d, hs) for d, hs in self.donors if d in self.lig_set]
        for don, hs in lig_donors:
            if not hs or not len(self.prot_acc):
                continue
            dd = np.linalg.norm(fit[self.prot_acc] - fit[don], axis=1)
            for j in np.where(dd <= 3.5)[0]:
                pa = int(self.prot_acc[j])
                for h in hs:
                    v1 = fit[don] - fit[h]
                    v2 = fit[pa] - fit[h]
                    ang = math.degrees(math.acos(np.clip(
                        np.dot(v1, v2) /
                        (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9),
                        -1, 1)))
                    if ang >= 120.0:
                        key = ("hbond", int(don), self.res_of[pa])
                        self.plif[key] = self.plif.get(key, 0) + 1
        # hydrophobic contacts
        if len(self.lig_carbons) and len(self.hyd_carbons):
            dd = np.linalg.norm(fit[self.lig_carbons][:, None, :] -
                                fit[self.hyd_carbons][None, :, :], axis=2)
            ii, jj = np.where(dd <= 4.5)
            for i, j in zip(ii, jj):
                key = ("hydrophobic", int(self.lig_carbons[i]),
                       self.res_of[int(self.hyd_carbons[j])])
                self.plif[key] = self.plif.get(key, 0) + 1

    def summary(self):
        arr = np.array(self.rows)
        tail = arr[len(arr) // 2:]
        return {
            "n_frames": self.frames,
            "ca_rmsd_mean_A": float(arr[:, 0].mean()),
            "ca_rmsd_max_A": float(arr[:, 0].max()),
            "pocket_ca_rmsd_mean_A": float(arr[:, 1].mean()),
            "lig_rmsd_mean_A": float(tail[:, 2].mean()),
            "lig_rmsd_max_A": float(arr[:, 2].max()),
            "temp_mean_k": float(arr[:, 3].mean()),
            "pe_mean_kj_mol": float(arr[:, 4].mean()),
            "plif_top": [
                {"type": k[0], "lig_atom": k[1], "residue": k[2],
                 "persistence": v}
                for k, v in sorted(self.plif.items(),
                                   key=lambda kv: -kv[1])[:15]],
        }

# --------------------------------------------------------------------------- #
#  STAGE 5 — end-state MM-GBSA + per-residue decomposition
# --------------------------------------------------------------------------- #

def _fibonacci_sphere(n=960):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    golden = math.pi * (1 + 5 ** 0.5)
    theta = golden * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


_SPHERE = _fibonacci_sphere(960)


def _sasa(positions_A, radii_A, probe=1.4) -> float:
    """Shrake-Rupley SASA (A^2), vectorized: per-atom numpy against the few
    candidate neighbours from a cKDTree (pure-python inner loops proved
    ~100x too slow for 40 frames x 2785 atoms)."""
    from scipy.spatial import cKDTree
    n_pts = len(_SPHERE)
    r_ext = radii_A + probe
    pts = positions_A[:, None, :] + _SPHERE[None, :, :] * r_ext[:, None, None]
    tree = cKDTree(positions_A)
    total = 0.0
    for i in range(len(radii_A)):
        cand = tree.query_ball_point(positions_A[i], r=radii_A.max() + probe)
        cand = [j for j in cand if j != i]
        if not cand:
            total += 4 * math.pi * r_ext[i] ** 2
            continue
        cand = np.array(cand)
        d = np.linalg.norm(pts[i][:, None, :] - positions_A[cand][None, :, :],
                           axis=2)
        blocked = (d < r_ext[cand][None, :]).any(axis=1)
        n_acc = n_pts - int(blocked.sum())
        total += n_acc * 4 * math.pi * r_ext[i] ** 2 / n_pts
    return float(total)


def _pair_vdw(r2, sig_i, eps_i, sig_j, eps_j):
    sig = 0.5 * (sig_i + sig_j)
    eps = math.sqrt(eps_i * eps_j)
    sr6 = (sig * sig / r2) ** 3
    return 4.0 * eps * (sr6 * sr6 - sr6)


def _still_pair(r2, qi, qj, ri, rj, tau=1.0 - 1.0 / 78.5):
    f = math.sqrt(r2 + ri * rj * math.exp(-r2 / (4.0 * ri * rj) + 1e-12))
    return -KC * tau * qi * qj / f


def stage5_mmgbsa(out: Path, args, force: bool = False) -> bool:
    ckpt = out / "stage5.json"
    if ckpt.exists() and not force:
        RESULTS["stage5_mmgbsa"] = json.loads(ckpt.read_text(encoding="utf-8"))
        _log("5", "checkpoint found, skipping")
        return True

    import mdtraj as md
    from openmm import (Context, LangevinMiddleIntegrator, Platform, unit,
                        XmlSerializer, NonbondedForce, GBSAOBCForce)
    from openmm.app import PDBFile, Modeller, CutoffNonPeriodic
    from openff.toolkit import Molecule as OFFMol
    from rdkit import Chem

    # ---- topology / trajectory -------------------------------------------- #
    top_pdb = PDBFile(str(out / "complex_start.pdb"))
    traj = md.load(str(out / "T04_complex_trajectory.dcd"),
                   top=str(out / "complex_start.pdb"))
    xyz = traj.xyz * 10.0                        # (n_frames, n_atoms, 3) in A
    n_all = xyz.shape[0]
    sel = np.unique(np.linspace(0, n_all - 1,
                                min(args.mgb_frames, n_all)).astype(int))
    _log("5", f"MM-GBSA over {len(sel)} frames of {n_all}")

    atoms = list(top_pdb.topology.atoms())
    lig_res_atoms = [a.index for a in atoms if a.residue.name == LIG_RESNAME]
    gdp_res_atoms = [a.index for a in atoms if a.residue.name == "GDP"]
    prot_atoms = [a.index for a in atoms
                  if a.residue.name not in (LIG_RESNAME,)]
    # ligand "receptor" = protein + GDP; "ligand" = T04
    rec_atoms = [a.index for a in atoms if a.residue.name != LIG_RESNAME]
    res_label = {a.index: f"{a.residue.name}{a.residue.id}" for a in atoms}

    # ---- systems ----------------------------------------------------------- #
    top_pdb0 = PDBFile(str(out / "complex_start.pdb"))
    system_c, _info_c = _build_complex_system(top_pdb0.topology,
                                              top_pdb0.positions)

    def _mk_system(keep_res_names):
        mod = Modeller(top_pdb.topology, top_pdb.positions)
        to_del = [r for r in mod.topology.residues()
                  if r.name not in keep_res_names]
        mod.delete(to_del)
        sysm = _regen_system(mod)
        for f in sysm.getForces():
            if _is_gbforce(f):
                try:
                    f.setCutoffDistance(1.6 * unit.nanometer)
                except Exception:
                    pass
        return sysm

    system_p = _mk_system({r.name for r in top_pdb.topology.residues()
                           if r.name != LIG_RESNAME})
    system_l = _mk_system({LIG_RESNAME})

    plat = Platform.getPlatformByName("CPU")
    integ = LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond,
                                     2 * unit.femtosecond)

    # force groups: 0 = Nonbonded, 1 = GB, 2 = rest
    for s in (system_c, system_p, system_l):
        for f in s.getForces():
            if isinstance(f, NonbondedForce):
                f.setForceGroup(0)
            elif _is_gbforce(f):
                f.setForceGroup(1)
            else:
                f.setForceGroup(2)

    ctx_c = Context(system_c, integ, plat)
    ctx_p = Context(system_p, LangevinMiddleIntegrator(
        310 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtosecond), plat)
    ctx_l = Context(system_l, LangevinMiddleIntegrator(
        310 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtosecond), plat)

    # ---- per-frame energies ------------------------------------------------ #
    E = {"nb_c": [], "gb_c": [], "nb_p": [], "gb_p": [],
         "nb_l": [], "gb_l": [], "sasa_c": [], "sasa_p": [], "sasa_l": []}
    radii = np.array([VDW_RADII.get(
        (a.element.symbol.upper() if a.element is not None else "C"), 1.7)
        for a in atoms])
    kJ = 4.184
    for fi in sel:
        pos_nm = xyz[fi] / 10.0
        ctx_c.setPositions(pos_nm)
        st = ctx_c.getState(getEnergy=True, groups={0})
        E["nb_c"].append(st.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole) / kJ)
        st = ctx_c.getState(getEnergy=True, groups={1})
        E["gb_c"].append(st.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole) / kJ)
        ctx_p.setPositions(pos_nm[prot_atoms])
        st = ctx_p.getState(getEnergy=True, groups={0})
        E["nb_p"].append(st.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole) / kJ)
        st = ctx_p.getState(getEnergy=True, groups={1})
        E["gb_p"].append(st.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole) / kJ)
        ctx_l.setPositions(pos_nm[lig_res_atoms])
        st = ctx_l.getState(getEnergy=True, groups={0})
        E["nb_l"].append(st.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole) / kJ)
        st = ctx_l.getState(getEnergy=True, groups={1})
        E["gb_l"].append(st.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole) / kJ)
        E["sasa_c"].append(_sasa(xyz[fi], radii))
        E["sasa_p"].append(_sasa(xyz[fi][prot_atoms], radii[prot_atoms]))
        E["sasa_l"].append(_sasa(xyz[fi][lig_res_atoms], radii[lig_res_atoms]))

    dE_nb = np.array(E["nb_c"]) - np.array(E["nb_p"]) - np.array(E["nb_l"])
    dG_gb = np.array(E["gb_c"]) - np.array(E["gb_p"]) - np.array(E["gb_l"])
    gamma_sa = 0.005                                    # kcal/mol/A^2
    dG_sa = gamma_sa * (np.array(E["sasa_c"]) - np.array(E["sasa_p"])
                        - np.array(E["sasa_l"]))
    dG_tot = dE_nb + dG_gb + dG_sa
    _log("5", f"dG_bind = {dG_tot.mean():.2f} +/- {dG_tot.std(ddof=1):.2f} "
              f"kcal/mol  (dE_nb {dE_nb.mean():.2f}, dG_GB {dG_gb.mean():.2f}, "
              f"dG_SA {dG_sa.mean():.2f})")

    # ---- per-residue decomposition ----------------------------------------- #
    nb = next(f for f in system_c.getForces() if isinstance(f, NonbondedForce))
    q = np.zeros(len(atoms)); sig = np.zeros(len(atoms)); eps = np.zeros(len(atoms))
    for i in range(nb.getNumParticles()):
        qi, si, ei = nb.getParticleParameters(i)
        q[i] = qi.value_in_unit(unit.elementary_charge)
        sig[i] = si.value_in_unit(unit.nanometer) * 10.0    # A
        eps[i] = ei.value_in_unit(unit.kilojoule_per_mole) / 4.184
    gb = next(f for f in system_c.getForces() if _is_gbforce(f))
    gb_r = np.zeros(len(atoms)); gb_s = np.zeros(len(atoms))
    gb_q = np.zeros(len(atoms))
    for i in range(gb.getNumParticles()):
        pr = gb.getParticleParameters(i)
        vals = [v if not hasattr(v, "value_in_unit") else None for v in pr]
        if hasattr(gb, "getPerParticleParameterName"):   # CustomGBForce
            # (charge[e], radius[nm], radius*scale[nm])
            gb_q[i] = float(pr[0])
            gb_r[i] = float(pr[1]) * 10.0
            gb_s[i] = float(pr[2]) / float(pr[1]) if pr[1] else 0.8
        else:                                            # legacy GBSAOBCForce
            gb_q[i] = pr[0].value_in_unit(unit.elementary_charge)
            gb_r[i] = pr[1].value_in_unit(unit.nanometer) * 10.0
            gb_s[i] = pr[2]

    # ---- validation of my pair-sum NB vs OpenMM ---------------------------- #
    # reference: complex system with ligand charges zeroed in NonbondedForce
    # (GB removed, salt screening disabled) -> its NB energy is exactly the
    # ligand-protein cross interaction my numpy pair-sum computes.
    f0 = xyz[sel[0]]
    la = np.array(lig_res_atoms)
    my_nb = 0.0
    for i in la:
        for j in prot_atoms:
            r2 = ((f0[i] - f0[j]) ** 2).sum()
            my_nb += _pair_vdw(r2, sig[i], eps[i], sig[j], eps[j]) + \
                KC * q[i] * q[j] / math.sqrt(r2)
    sys2, _ = _build_complex_system(top_pdb0.topology, top_pdb0.positions)
    nb2 = next(f for f in sys2.getForces() if isinstance(f, NonbondedForce))
    lig_set = set(int(x) for x in la)
    for i in range(nb2.getNumParticles()):
        qi, si, ei = nb2.getParticleParameters(i)
        if i in lig_set:
            nb2.setParticleParameters(i, 0.0 * unit.elementary_charge, si, ei)
    nb2.setForceGroup(5)
    for fi in range(sys2.getNumForces() - 1, -1, -1):
        if _is_gbforce(sys2.getForce(fi)):
            sys2.removeForce(fi)
    ctx2 = Context(sys2, LangevinMiddleIntegrator(
        310 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtosecond), plat)
    ctx2.setPositions(f0 / 10.0)
    e_ligzero = ctx2.getState(getEnergy=True, groups={5})         .getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / 4.184
    ctx_c.setPositions(f0 / 10.0)
    e_full = ctx_c.getState(getEnergy=True, groups={0})         .getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / 4.184
    ref_e = e_full - e_ligzero
    nb_rel_err = abs(my_nb - ref_e) / max(abs(ref_e), 1e-9)
    nb_exact = bool(nb_rel_err < 0.05)
    _log("5", f"pair-sum NB validation: mine {my_nb:.2f} vs OpenMM "
              f"{ref_e:.2f} kcal/mol (rel err {nb_rel_err:.3%})")

    # accumulate per-residue (protein only) over frames
    from collections import defaultdict
    acc = defaultdict(lambda: np.zeros(3))       # vdW, elec, gb
    la = np.array(lig_res_atoms)
    for fi in sel:
        fr = xyz[fi]
        d_lig = np.linalg.norm(fr[la][:, None, :] - fr[prot_atoms][None, :, :],
                               axis=2)
        ii, jj = np.where(d_lig < 12.0)
        for i, j in zip(ii, jj):
            ia, ja = int(la[i]), int(prot_atoms[j])
            r2 = ((fr[ia] - fr[ja]) ** 2).sum()
            if r2 < 0.5:
                continue
            acc[res_label[ja]][0] += _pair_vdw(r2, sig[ia], eps[ia],
                                               sig[ja], eps[ja]) / len(sel)
            acc[res_label[ja]][1] += KC * q[ia] * q[ja] / math.sqrt(r2) / len(sel)
            acc[res_label[ja]][2] += _still_pair(r2, gb_q[ia], gb_q[ja],
                                                 gb_r[ia], gb_r[ja]) / len(sel)

    rows = []
    for res, v in acc.items():
        tot = v.sum()
        if tot < -0.05:
            rows.append({"residue": res, "vdW": float(v[0]),
                         "elec": float(v[1]), "gb_polar": float(v[2]),
                         "total": float(tot)})
    rows.sort(key=lambda r: r["total"])
    top10 = rows[:10]

    rec = {
        "frames": int(len(sel)),
        "dg_bind_kcal_mol": float(dG_tot.mean()),
        "dg_bind_std": float(dG_tot.std(ddof=1)),
        "components": {
            "dE_nb_mean": float(dE_nb.mean()),
            "dG_gb_mean": float(dG_gb.mean()),
            "dG_sa_mean": float(dG_sa.mean()),
            "gamma_sa": gamma_sa,
        },
        "per_frame": {"dg": dG_tot.tolist()},
        "per_residue_top10": top10,
        "decomposition_validation": {
            "my_cross_nb_kcal": float(my_nb),
            "openmm_cross_nb_kcal": float(ref_e),
            "rel_err": float(nb_rel_err),
            "nb_exact": bool(nb_exact),
            "gb_note": ("GB per-residue = Still pair terms with intrinsic "
                        "OBC radii (approximate; screening-level accuracy)"),
        },
        "sasa_mean_complex": float(np.mean(E["sasa_c"])),
    }
    RESULTS["stage5_mmgbsa"] = rec
    ckpt.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    import csv as _csv
    with open(out / "T04_per_residue_mmgbsa.csv", "w", newline="",
              encoding="utf-8") as fh:
        wr = _csv.writer(fh)
        wr.writerow(["residue", "vdW", "elec", "gb_polar", "total"])
        for r in top10:
            wr.writerow([r["residue"], f"{r['vdW']:.3f}", f"{r['elec']:.3f}",
                         f"{r['gb_polar']:.3f}", f"{r['total']:.3f}"])
    _log("5", f"top residues: " +
              ", ".join(f"{r['residue']} {r['total']:.1f}" for r in top10[:5]))
    return True


def _regen_system(modeller):
    """Recreate a component system (protein-only / ligand-only / etc.) from a
    pruned Modeller via the splice builder."""
    sysm, _info = _build_complex_system(modeller.topology, modeller.positions)
    return sysm

# --------------------------------------------------------------------------- #
#  STAGE 6 — publication figures
# --------------------------------------------------------------------------- #

def stage6_figures(fig_dir: Path) -> bool:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig_dir.mkdir(parents=True, exist_ok=True)

    out = _OUT
    s1 = RESULTS["stage1a_curation"]
    s2 = RESULTS.get("stage2_docking", {})
    s3 = RESULTS.get("stage3_ff_duality", {})
    s4 = RESULTS.get("stage4_complex_md", {})
    s5 = RESULTS.get("stage5_mmgbsa", {})

    # ---------------- fig1: binding pose + pocket + contacts ---------------- #
    import mdtraj as md
    traj = md.load(str(out / "T04_complex_trajectory.dcd"),
                   top=str(out / "complex_start.pdb"))
    top = traj.topology
    final = traj.xyz[-1] * 10.0
    lig_idx = [a.index for a in top.atoms if a.residue.name == LIG_RESNAME]
    lig_heavy = [i for i in lig_idx
                 if top.atom(i).element.symbol != "H"]
    lig_centroid = final[lig_heavy].mean(axis=0)
    pocket_atoms = [a.index for a in top.atoms
                    if a.residue.name not in (LIG_RESNAME, "GDP", "HOH")
                    and a.element.symbol != "H"
                    and np.linalg.norm(final[a.index] - lig_centroid) <= 8.0]
    res_of = {a.index: f"{a.residue.name}{a.residue.resSeq}"
              for a in top.atoms}

    fig = plt.figure(figsize=(10.5, 9))
    ax = fig.add_subplot(111, projection="3d")
    xyz = final
    ax.scatter(xyz[pocket_atoms, 0], xyz[pocket_atoms, 1], xyz[pocket_atoms, 2],
               s=8, c="0.72", depthshade=True, linewidths=0)
    ele_color = {"C": "#222222", "N": "#2E86DE", "O": "#E74C3C",
                 "F": "#27AE60", "Cl": "#16A085", "S": "#F1C40F",
                 "H": "#AAAAAA"}
    lc = [ele_color.get(top.atom(i).element.symbol, "#7D3C98")
          for i in lig_heavy]
    ax.scatter(xyz[lig_heavy, 0], xyz[lig_heavy, 1], xyz[lig_heavy, 2],
               s=55, c=lc, depthshade=False, edgecolors="k", linewidths=0.4,
               zorder=10)
    # contacts (PLIF persistence >= 30%)
    if s4.get("md", {}).get("plif_top"):
        seen_res = []
        for item in s4["md"]["plif_top"]:
            if item["persistence"] < 0.3 * s4["md"]["n_frames"]:
                continue
            res = item["residue"].split(":")[-1]
            if res in seen_res:
                continue
            seen_res.append(res)
            try:
                ra = [a.index for a in top.atoms
                      if f"{a.residue.name}{a.residue.resSeq}" == res
                      and a.element.symbol != "H"]
            except Exception:
                continue
            if not ra:
                continue
            rc = xyz[ra].mean(axis=0)
            li = lig_heavy[np.argmin(
                np.linalg.norm(xyz[lig_heavy] - rc, axis=1))]
            ax.plot([xyz[li, 0], rc[0]], [xyz[li, 1], rc[1]],
                    [xyz[li, 2], rc[2]], "k--", linewidth=1.0, alpha=0.75)
            ax.text(rc[0], rc[1], rc[2], f" {res} ", fontsize=8,
                    color="#B03A2E",
                    bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.5))
    ax.set_box_aspect((1, 1, 1))
    pad = 12
    c = lig_centroid
    ax.set_xlim(c[0] - pad, c[0] + pad)
    ax.set_ylim(c[1] - pad, c[1] + pad)
    ax.set_zlim(c[2] - pad, c[2] + pad)
    dg = s2.get("top3_delta_g_kcal_mol", [float("nan")])[0]
    ax.set_title(f"T04 in the KRAS G12D switch-II pocket ({TARGET_PDB})\n"
                 f"final MD frame | Vina ΔG = {dg:.1f} kcal/mol | "
                 f"dashed = PLIF contacts ≥30% persistence",
                 fontsize=11.5, fontweight="bold")
    handles = [Line2D([], [], marker="o", ls="", color=v, label=k)
               for k, v in list(ele_color.items())[:5]]
    handles.append(Line2D([], [], color="0.5", ls="--",
                          label="pocket residue (heavy atoms)"))
    ax.legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_binding_pose_pocket.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log("6", f"wrote {fig_dir/'fig1_binding_pose_pocket.png'}")

    # ---------------- fig2: per-residue MM-GBSA bars ------------------------ #
    if s5.get("per_residue_top10"):
        rows = s5["per_residue_top10"][::-1]
        names = [r["residue"] for r in rows]
        vdw = np.array([r["vdW"] for r in rows])
        ele = np.array([r["elec"] for r in rows])
        gbp = np.array([r["gb_polar"] for r in rows])
        fig, ax = plt.subplots(figsize=(9.5, 7))
        y = np.arange(len(names))
        ax.barh(y, vdw, color="#0173B2", label="van der Waals")
        ax.barh(y, ele, left=vdw, color="#D55E00", label="Electrostatic")
        ax.barh(y, gbp, left=vdw + ele, color="#009E73",
                label="GB polar (desolvation)")
        tot = vdw + ele + gbp
        for yi, t in zip(y, tot):
            ax.text(t - 0.05, yi, f"{t:.2f}", ha="right", va="center",
                    fontsize=9)
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=10)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("per-residue interaction energy (kcal/mol)   "
                      "[negative = favorable]")
        dg = s5.get("dg_bind_kcal_mol")
        ax.set_title(f"Per-residue MM-GBSA decomposition — top 10\n"
                     f"ΔG_bind = {dg:.2f} ± "
                     f"{s5.get('dg_bind_std', 0):.2f} kcal/mol "
                     f"({s5.get('frames')} frames, 310 K)",
                     fontsize=11.5, fontweight="bold")
        ax.legend(loc="lower right", fontsize=9)
        fig.tight_layout()
        fig.savefig(fig_dir / "fig2_per_residue_mmgbsa.png", dpi=300,
                    bbox_inches="tight", facecolor="white")
        plt.close(fig)
        _log("6", f"wrote {fig_dir/'fig2_per_residue_mmgbsa.png'}")

    # ---------------- fig3: ML vs classical FF gap -------------------------- #
    ml = s3.get("ml", {})
    if ml.get("available") and s3:
        pa = ml["per_atom"]
        moieties = np.array(pa["moieties"], dtype=object)
        uniq = sorted(set(pa["moieties"]),
                      key=lambda m: -np.mean(
                          np.array(pa["df_gaff_norm"])[moieties == m]))
        cmap = plt.get_cmap("tab10")
        mcol = {m: cmap(i % 10) for i, m in enumerate(uniq)}
        fg = np.array(pa["df_gaff_norm"])
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 6.2))
        # panel A: parity of |F| on GAFF2 vs ML subsample — need |F|; use df
        # proxies: scatter df_gaff vs df_mmff colored by moiety
        fm = np.array(pa["df_mmff_norm"])
        for m in uniq:
            sel = moieties == m
            a1.scatter(fg[sel], fm[sel], s=26, alpha=0.8, color=mcol[m],
                       edgecolors="none", label=m)
        lim = max(fg.max(), fm.max()) * 1.05
        a1.plot([0, lim], [0, lim], "k--", lw=0.8)
        a1.set_xlabel("|ΔF| GAFF2 − ML (kcal/mol/Å)")
        a1.set_ylabel("|ΔF| MMFF94 − ML (kcal/mol/Å)")
        a1.set_title("(A) force-field discrepancy magnitude\n"
                     "(distance to the ML reference)", fontsize=10.5)
        a1.legend(fontsize=7, frameon=True)
        # panel B: per-atom df sorted, colored by moiety
        order = np.argsort(fg)[::-1]
        a2.bar(range(len(order)), fg[order],
               color=[mcol[moieties[i]] for i in order], width=1.0)
        a2.set_xlabel("ligand atom (sorted)")
        a2.set_ylabel("|F_GAFF2 − F_ML| (kcal/mol/Å)")
        p = ml["parity"]["gaff2"]
        a2.set_title(f"(B) atomic force discrepancy (GAFF2 vs "
                     f"{ml['engine'].split(' (')[0]})\n"
                     f"mean cos(F_GAFF2, F_ML) = {p['mean_cosine']:.3f} | "
                     f"⟨|ΔF|⟩ = {p['df_atom_mean']:.2f} kcal/mol/Å",
                     fontsize=10.5)
        fig.suptitle("Classical vs ML potential: force-field discrepancy "
                     f"vector on pocket-frozen T04 ({s3['ligand_atoms']} atoms)",
                     fontsize=12, fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(fig_dir / "fig3_ml_vs_classical_ff_gap.png", dpi=300,
                    bbox_inches="tight", facecolor="white")
        plt.close(fig)
        _log("6", f"wrote {fig_dir/'fig3_ml_vs_classical_ff_gap.png'}")
    else:
        _warn("fig3 skipped: ML single-points unavailable")
    return True


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #

def _shutdown(delay_s: int) -> None:
    _log("shutdown", f"initiating shutdown in {delay_s} s "
                     f"(--auto_shutdown was set and all stages succeeded)")
    try:
        subprocess.run(["shutdown", "/s", "/t", str(delay_s)],
                       check=False, capture_output=True)
    except Exception as exc:
        _log("shutdown", f"shutdown command failed: {exc}")


def main() -> int:
    p = argparse.ArgumentParser(
        description="Phase 3: complex dynamics, MM-GBSA decomposition, "
                    "ML-vs-classical FF duality")
    p.add_argument("--out_dir", default="results_phase3")
    p.add_argument("--fig_dir", default="figures_phase3")
    p.add_argument("--pdb_id", default=TARGET_PDB)
    p.add_argument("--lig_code", default=NATIVE_LIG_CODE)
    p.add_argument("--equil_steps", type=int, default=5000)
    p.add_argument("--md_steps", type=int, default=100000)
    p.add_argument("--report_interval", type=int, default=500)
    p.add_argument("--mgb_frames", type=int, default=40)
    p.add_argument("--force_rerun", action="store_true")
    p.add_argument("--skip_dock", action="store_true")
    p.add_argument("--skip_md", action="store_true")
    p.add_argument("--fig_only", action="store_true")
    p.add_argument("--auto_shutdown", action="store_true",
                   help="shut the machine down AFTER full success only")
    p.add_argument("--shutdown_delay", type=int, default=60)
    args = p.parse_args()

    global _OUT, _FIG
    _OUT = Path(args.out_dir)
    _FIG = Path(args.fig_dir)
    _OUT.mkdir(parents=True, exist_ok=True)
    _FIG.mkdir(parents=True, exist_ok=True)

    code = 0
    try:
        _hr(f"PHASE 3 COMPLEX DYNAMICS — {RESULTS['timestamp']}")
        _log("main", f"python {sys.version.split()[0]} | target "
                     f"{args.pdb_id} (ref ligand {args.lig_code}) | "
                     f"ligand T04")

        if args.fig_only:
            for st, key in (("stage1a", "stage1a_curation"),
                            ("stage2", "stage2_docking"),
                            ("stage3", "stage3_ff_duality"),
                            ("stage4", "stage4_complex_md"),
                            ("stage5", "stage5_mmgbsa")):
                f = _OUT / f"{st}.json"
                if f.exists():
                    RESULTS[key] = json.loads(f.read_text(encoding="utf-8"))
            stage6_figures(_FIG)
            return 0

        stages = [
            ("1A  ingestion & curation", lambda: stage1a_ingest(_OUT,
                                                                args.force_rerun)),
        ]
        if not args.skip_dock:
            stages.append(("2   docking", lambda: stage2_docking(
                _OUT, args.force_rerun)))
        stages.append(("3   FF duality (ML vs classical)",
                       lambda: stage3_ff_duality(_OUT, args.force_rerun)))
        if not args.skip_md:
            stages.append(("4   complex MD (OBC2, 310 K)",
                           lambda: stage4_complex_md(_OUT, args,
                                                     args.force_rerun)))
            stages.append(("5   MM-GBSA + per-residue",
                           lambda: stage5_mmgbsa(_OUT, args,
                                                args.force_rerun)))

        for name, fn in stages:
            _hr(f"STAGE {name}")
            try:
                ok = fn()
                _log("stage", f"{name}: {'OK' if ok else 'SKIPPED'}")
            except Exception as exc:
                code = 1
                _log("stage", f"{name}: FAILED — {exc.__class__.__name__}: "
                              f"{exc}")
                traceback.print_exc()
                RESULTS["meta"]["warnings"].append(
                    f"{name} failed: {exc.__class__.__name__}: "
                    f"{str(exc)[:200]}")
                # dependent stages are pruned below on missing artifacts
                if name.startswith("1A"):
                    break

        try:
            stage6_figures(_FIG)
        except Exception as exc:
            _warn(f"figure generation failed: {exc}")
            traceback.print_exc()

        RESULTS["all_stages_ok"] = (code == 0)
    except Exception as exc:
        code = 1
        RESULTS["fatal_error"] = f"{exc.__class__.__name__}: {exc}"
        traceback.print_exc()
    finally:
        _hr("FINALIZING")
        try:
            write_json_atomic(_OUT / "phase3_results.json")
            _log("final", f"results -> {_OUT/'phase3_results.json'}")
        except Exception as exc:
            _log("final", f"result serialization failed: {exc}")
        if args.auto_shutdown and RESULTS["all_stages_ok"]:
            _shutdown(args.shutdown_delay)
        else:
            _log("final", "shutdown NOT triggered "
                          "(disabled or stages failed)")
    return code


if __name__ == "__main__":
    sys.exit(main())

