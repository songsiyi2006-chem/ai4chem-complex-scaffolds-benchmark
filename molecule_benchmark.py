#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
molecule_benchmark.py
=====================
Hardcore 10-molecule benchmark for novel, structurally complex entities.

Per molecule:
  1. Robust SMILES parse + sanitization (stepwise fallback pinpoints the
     failing stage), stereochemistry perception, Bemis-Murcko scaffold
     (native + generic), and atom/bond feature tensors for GNN input.
  2. 50-conformer ETKDGv3 ensemble (small-ring + macrocycle torsions on),
     MMFF94 optimization with UFF fallback, lowest-energy conformer
     identification (delta-E = 0 kcal/mol reference).
  3. Physicochemical profiling: MW, cLogP, TPSA, Fsp3, rotatable bonds,
     HBD/HBA, stereocenters, formal charges, and intramolecular
     hydrogen-bond (IMHB) geometry (d(H...A) <= 2.5 A, angle >= 120 deg,
     >= 3 bonds apart) in the lowest-energy conformer.
  4. Failure analysis for 2D GNNs: steric clash scan (min non-bonded
     heavy-atom distance), ring-strain inventory (3/4-rings, macrocycles),
     atropisomerism heuristic (ortho-blocked biaryl axes), zwitterion /
     perfluoro / flexibility flags, conformational entropy (ensemble
     delta-E).

Concurrency: ProcessPoolExecutor over molecules (RDKit embedding is CPU-bound).

Outputs (default ./bench_results):
  sdf/<ID>_ensemble.sdf    all optimized conformers (E as SD property)
  sdf/<ID>_min.sdf         lowest-energy conformer
  features/<ID>.npz        atom/bond feature tensors + edge index
  benchmark_results.json   full machine-readable records
  benchmark_report.md      Markdown benchmark table + per-molecule analysis
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

if sys.platform.startswith("win"):
    for _stream in (sys.stdout, sys.stderr):
        try:
            _stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

SEED = 0x2026

REGISTRY: List[Dict[str, str]] = [
    {"id": "M01", "name": "Macrocyclic chameleon peptidomimetic",
     "name_cn": "大环变色龙拟肽",
     "challenge": "Transannular H-bond networks, solvent-dependent conformational switching, extreme conformational entropy.",
     "smiles": "O=C1N[C@H](C(C)C)C(=O)N[C@@H](Cc2ccccc2)C(=O)N(C)[C@H](CC(=O)NC1)C(=O)N2CCC[C@H]2C(=O)O"},
    {"id": "M02", "name": "Strained azaspiro-cubane bioisostere",
     "name_cn": "稠合氮杂螺环-立方烷生物电子等排体",
     "challenge": "Extreme C-C-C angle distortion, bridgehead nitrogen strain, hard distance geometry.",
     "smiles": "O=C(N1CC23C4C1C2C34)c1c(F)c(F)c(N5CC6(COC6)C5)c(F)c1F"},
    {"id": "M03", "name": "Bivalent heterobifunctional degron (PROTAC prototype)",
     "name_cn": "新型刚柔偶联PROTAC分子",
     "challenge": "MW > 850 Da, 16+ rotatable linker bonds, inter-domain steric clashes on conformer collapse; input is a 2-fragment disconnected assembly.",
     "smiles": "O=C1c2ccccc2C(=O)N1C3CCC(=O)NC3=O.O=C(NCCOCCOCCOc4ccc(NC(=O)c5cnc(Nc6ccc(S(=O)(=O)C)cc6)nc5)cc4)C"},
    {"id": "M04", "name": "Axially + centrally chiral atropisomeric biaryl",
     "name_cn": "轴手性与点手性结合的双芳基抑制剂",
     "challenge": "High biaryl rotation barrier (>30 kcal/mol), perpendicular aromatic planes, non-planar conjugation.",
     "smiles": "CC(=O)Oc1c(C)c(c2c(OC(=O)C)c(C)ccc2[C@@H](C)NC(=O)CF)ccc1[C@H](C)O"},
    {"id": "M05", "name": "Perfluorinated polycyclic cage vector",
     "name_cn": "全氟代多环笼状分子",
     "challenge": "Dense fluorine electrostatics, extreme hydrophobicity, abnormal surface/volume ratio, electron-deficient core.",
     "smiles": "FC1(F)C2(F)C3(F)C1(F)C4(F)C2(F)C3(F)C4(C(=O)NC5(CC5)C(=O)O)(F)"},
    {"id": "M06", "name": "Oxetane-fused polyketide NP mimic",
     "name_cn": "含氧杂四元环的多手性中心天然产物类似物",
     "challenge": "7 stereocenters, constrained oxygenated small ring, dense intramolecular steric repulsion.",
     "smiles": "C[C@H]1O[C@@]2(CO2)[C@@H](O)[C@H](C)[C@@H](OC(=O)c3ccccc3)[C@H]1C(=O)N[C@H](C)c4nc(C)cs4"},
    {"id": "M07", "name": "Zwitterionic B-N dative macrocycle",
     "name_cn": "含B-N配位内盐大环化合物",
     "challenge": "Dative B-N bond, formal B(-)/N(+) charges, tetrahedral borate coordination; MMFF lacks B params (UFF fallback).",
     "smiles": "c1ccc2c(c1)[B-]3(c4ccccc42)OCCN[N+]3=Cc5cccc(O)c5"},
    {"id": "M08", "name": "Strained bicyclo-acrylamide covalent warhead",
     "name_cn": "张力双环丙烷共价弹头分子",
     "challenge": "Electrophilic warhead alignment, localized ring strain, reactive dihedral vectoring for cysteine trapping.",
     "smiles": "C=CC(=O)N1CC2(CC12)C(=O)Nc3ccc(OC(F)(F)F)c(Cl)c3"},
    {"id": "M09", "name": "Heteroatom-doped [5]-carbohelicene",
     "name_cn": "杂原子掺杂五螺烯手性发光分子",
     "challenge": "Helical non-planar aromatic distortion, overlapping terminal rings, optical asymmetry.",
     "smiles": "c1cc2c(s1)c3ccc4c(c3c2)c5ccc6ncccc6c5c7cccnc47"},
    {"id": "M09R", "name": "Aza-[5]-carbohelicene (repaired reference for M09)",
     "name_cn": "氮杂[5]螺烯（M09修复参照）",
     "challenge": "Reference structure: the M09 SMILES as supplied is unkekulizable (invalid ring-fusion "
                  "aromaticity); this program-generated aza-doped [5]-helicene (C21H13N) demonstrates the "
                  "intended helical topology and failure modes.",
     "smiles": "c1ccc2c(c1)ccc1c2ccc2c3ccncc3ccc21"},
    {"id": "M10", "name": "Tetra-ortho hindered peptoid core",
     "name_cn": "四邻位超拥挤拟肽骨架",
     "challenge": "Severe steric clash blocking amide cis/trans planarization; hindered-rotation energy barriers.",
     "smiles": "CC1=C(C)C(=C(C)C(=C1C)N(C)C(=O)CN(C)C(=O)c2c(C)c(C)c(C)c(C)c2C)C"},
]

ATOM_FEATURES = ["atomic_num", "heavy_degree", "formal_charge", "chiral_tag",
                 "total_num_Hs", "hybridization", "is_aromatic", "is_in_ring"]
BOND_FEATURES = ["bond_type", "is_conjugated", "is_in_ring", "stereo"]

FLAG_LEGEND = {
    "MAC": "macrocycle (>=12-ring): transannular effects invisible to 2D GNNs",
    "STRAIN": "contains 3/4-membered rings: high angle strain, distance-geometry edge case",
    "ATRO": "ortho-blocked biaryl axis (>=3/4 positions): atropisomerism, axial chirality lost in 2D",
    "FLUORO": ">=6 fluorines: dense C-F electrostatics, LogP/vdW model extrapolation",
    "ZWIT": "formal charges / zwitterion / dative bonding: OOD for typical pretrained 2D GNNs",
    "BIG": "MW > 800 Da: outside Ro5 / typical GNN pretraining domain",
    "FLEX": ">=15 rotatable bonds: conformational entropy, single-graph under-sampling",
    "FRAG": "disconnected multi-fragment input: 3D run on largest fragment",
    "ENTROPY": "ensemble dE >= 10 kcal/mol: extreme conformational polymorphism",
    "CLASH": "non-bonded heavy-atom pair < 2.5 A in Emin conformer: steric clash",
    "UFF": "MMFF94 params unavailable -> UFF fallback (lower accuracy)",
    "UST": "unassigned stereocenters detected",
}


# --------------------------------------------------------------------------- #
# chemistry helpers
# --------------------------------------------------------------------------- #
def robust_parse(smiles: str) -> Tuple[Optional[Any], str]:
    """Parse + sanitize with a stepwise fallback that pinpoints failures.

    Args:
        smiles: input SMILES.

    Returns:
        (mol or None, human-readable status).
    """
    from rdkit import Chem

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return mol, "sanitized"
    except Exception:
        pass

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None, "SMILES parse error"
    from rdkit import Chem
    from rdkit.Chem import rdmolops

    # diagnostic stepwise sanitization; functions that no longer exist in
    # newer RDKit builds are skipped automatically
    step_specs = [
        ("UpdatePropertyCache", None, None),
        ("SymmetrizeSSSR", "SymmetrizeSSSR", None),
        ("AssignRadicals", "AssignRadicals", None),
        ("Kekulize", "Kekulize", {"clearAromaticFlags": True}),
        ("SetAromaticity", "SetAromaticity", None),
        ("SetConjugation", "SetConjugation", None),
        ("SetHybridization", "SetHybridization", None),
        ("AdjustHs", "AdjustHs", None),
    ]
    steps: List[Tuple[str, Any]] = []
    for label, fname, kwargs in step_specs:
        if fname is None:
            steps.append((label, lambda m: m.UpdatePropertyCache(strict=False)))
            continue
        fn = getattr(rdmolops, fname, None)
        if fn is None:
            continue
        steps.append((label, (lambda m, fn=fn, kw=kwargs: fn(m, **kw)) if kwargs else fn))
    for name, fn in steps:
        try:
            fn(mol)
        except Exception as exc:
            return None, f"sanitization failed at {name}: {exc}"
    return mol, "sanitized (stepwise fallback)"


def murcko_analysis(mol: Any) -> Dict[str, Any]:
    """Native + generic Bemis-Murcko scaffolds and scaffold atom fraction."""
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    return {
        "scaffold_smiles": Chem.MolToSmiles(scaffold),
        "generic_smiles": Chem.MolToSmiles(generic),
        "scaffold_atom_fraction": round(scaffold.GetNumAtoms() / max(mol.GetNumAtoms(), 1), 3),
    }


def build_gnn_features(mol: Any) -> Dict[str, Any]:
    """Integer-encoded atom/bond feature tensors + edge index (both directions).

    Saved as .npz; directly loadable by torch/PyG
    (Data(x=atom_feat.float(), edge_index=edge_index, edge_attr=bond_feat.float())).
    """
    import numpy as np
    from rdkit import Chem

    atom_feat = []
    for a in mol.GetAtoms():
        atom_feat.append([
            float(a.GetAtomicNum()),
            float(a.GetDegree()),
            float(a.GetFormalCharge()),
            float(int(a.GetChiralTag())),
            float(a.GetTotalNumHs()),
            float(int(a.GetHybridization())),
            float(a.GetIsAromatic()),
            float(a.IsInRing()),
        ])
    edge_index, bond_feat = [], []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        feats = [float(int(b.GetBondType())), float(b.GetIsConjugated()),
                 float(b.IsInRing()), float(int(b.GetStereo()))]
        edge_index += [[i, j], [j, i]]      # undirected message passing
        bond_feat += [feats, feats]
    return {
        "atom_feat": np.asarray(atom_feat, dtype=np.float32),
        "edge_index": np.asarray(edge_index, dtype=np.int64).T.reshape(2, -1) if edge_index
                      else np.zeros((2, 0), dtype=np.int64),
        "bond_feat": np.asarray(bond_feat, dtype=np.float32) if bond_feat
                     else np.zeros((0, 4), dtype=np.float32),
    }


def gnn_readiness_smoke(atom_feat, edge_index, bond_feat) -> Dict[str, Any]:
    """One GCNConv forward pass through the molecule graph if torch/PyG exist."""
    try:
        import torch
        import torch_geometric
        from torch_geometric.data import Data
        from torch_geometric.nn import GCNConv

        data = Data(x=torch.from_numpy(atom_feat).float(),
                    edge_index=torch.from_numpy(edge_index),
                    edge_attr=torch.from_numpy(bond_feat).float())
        out = GCNConv(data.num_node_features, 64)(data.x, data.edge_index)
        return {"ok": True, "torch": torch.__version__,
                "pyg": torch_geometric.__version__,
                "embedding_shape": list(out.shape)}
    except Exception as exc:
        return {"ok": False, "reason": f"{type(exc).__name__}: {exc}"}


def embed_and_optimize(mol: Any, num_confs: int) -> Dict[str, Any]:
    """ETKDGv3 ensemble + MMFF94 (UFF fallback) optimization.

    Returns dict with mol_h (H-added mol), conformer ids, FF name, energies
    (kcal/mol), and convergence stats.
    """
    from rdkit.Chem import AllChem

    mol_h = AllChem.AddHs(mol)
    params = AllChem.ETKDGv3()
    for attr, val in (("randomSeed", SEED), ("useSmallRingTorsions", True),
                      ("useMacrocycleTorsions", True), ("pruneRmsThresh", -1.0),
                      ("numThreads", 1), ("timeout", 300)):
        try:
            setattr(params, attr, val)
        except Exception:
            pass

    cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=params))
    used_random_coords = False
    if not cids:  # desperate retry for pathological ring systems
        used_random_coords = True
        try:
            params.useRandomCoords = True
            params.randomSeed = SEED + 1
        except Exception:
            pass
        cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=params))
    if not cids:
        return {"ok": False, "reason": "ETKDGv3 embedding failed (both deterministic and random-coordinate)"}

    if AllChem.MMFFHasAllMoleculeParams(mol_h):
        ff = "MMFF94"
        res = AllChem.MMFFOptimizeMoleculeConfs(mol_h, numThreads=1, maxIters=2000)
    else:
        ff = "UFF"
        res = AllChem.UFFOptimizeMoleculeConfs(mol_h, numThreads=1, maxIters=2000)

    energies = [float(e) for _, e in res]
    not_converged = sum(1 for flag, _ in res if flag != 0)
    ids = [conf.GetId() for conf in mol_h.GetConformers()]
    best_i = min(range(len(energies)), key=lambda k: energies[k])
    return {
        "ok": True, "mol_h": mol_h, "ff": ff, "n_embedded": len(cids),
        "conf_ids": ids, "energies": energies, "not_converged": not_converged,
        "e_min": energies[best_i], "e_min_conf_id": ids[best_i],
        "e_max": max(energies), "delta_e": max(energies) - energies[best_i],
        "used_random_coords": used_random_coords,
    }


def _dist(p, q) -> float:
    return math.sqrt((p.x - q.x) ** 2 + (p.y - q.y) ** 2 + (p.z - q.z) ** 2)


def _angle(p1, p2, p3) -> float:
    """Angle p1-p2-p3 in degrees (vertex at p2)."""
    v1 = (p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)
    v2 = (p3.x - p2.x, p3.y - p2.y, p3.z - p2.z)
    n1 = math.sqrt(sum(c * c for c in v1))
    n2 = math.sqrt(sum(c * c for c in v2))
    cosang = max(-1.0, min(1.0, sum(a * b for a, b in zip(v1, v2)) / (n1 * n2 + 1e-12)))
    return math.degrees(math.acos(cosang))


def imhb_analysis(mol_h: Any, conf: Any) -> Dict[str, Any]:
    """Intramolecular H-bonds: d(H...A) <= 2.5 A, angle(D-H...A) >= 120 deg,
    donor/acceptor separated by >= 3 heavy-atom bonds."""
    from rdkit import Chem

    donors, acceptors = [], []
    for a in mol_h.GetAtoms():
        if a.GetSymbol() not in ("N", "O"):
            continue
        if any(nb.GetAtomicNum() == 1 for nb in a.GetNeighbors()):
            donors.append(a.GetIdx())
        if a.GetFormalCharge() <= 0:
            acceptors.append(a.GetIdx())
    if not donors or not acceptors:
        return {"count": 0, "bonds": []}

    dm = Chem.GetDistanceMatrix(mol_h)
    pos = {i: conf.GetAtomPosition(i) for i in range(mol_h.GetNumAtoms())}
    hits = []
    for d_idx in donors:
        h_ids = [nb.GetIdx() for nb in mol_h.GetAtomWithIdx(d_idx).GetNeighbors()
                 if nb.GetAtomicNum() == 1]
        for a_idx in acceptors:
            if a_idx == d_idx or dm[d_idx][a_idx] < 3 or dm[d_idx][a_idx] > 1e6:
                continue
            for h_idx in h_ids:
                d_ha = _dist(pos[h_idx], pos[a_idx])
                if d_ha > 2.5:
                    continue
                angle = _angle(pos[d_idx], pos[h_idx], pos[a_idx])
                if angle >= 120.0:
                    hits.append({"donor": int(d_idx), "h": int(h_idx), "acceptor": int(a_idx),
                                 "dHA_A": round(d_ha, 2), "angle_deg": round(angle, 1)})
    # de-duplicate on (donor, acceptor), keep the shortest H...A contact
    best: Dict[Tuple[int, int], Dict[str, Any]] = {}
    for hb in hits:
        key = (hb["donor"], hb["acceptor"])
        if key not in best or hb["dHA_A"] < best[key]["dHA_A"]:
            best[key] = hb
    return {"count": len(best), "bonds": sorted(best.values(), key=lambda h: h["dHA_A"])}


def steric_scan(mol_h: Any, conf: Any) -> Dict[str, Any]:
    """Minimum heavy-atom distance for topologically separated (>=4 bonds) pairs."""
    from rdkit import Chem

    heavy = [a.GetIdx() for a in mol_h.GetAtoms() if a.GetAtomicNum() > 1]
    dm = Chem.GetDistanceMatrix(mol_h)
    pos = {i: conf.GetAtomPosition(i) for i in range(mol_h.GetNumAtoms())}
    min_dist, pair, clashes = None, None, 0
    for ii, i in enumerate(heavy):
        for j in heavy[ii + 1:]:
            if dm[i][j] < 4 or dm[i][j] > 1e6:
                continue
            d = _dist(pos[i], pos[j])
            if min_dist is None or d < min_dist:
                min_dist, pair = d, (i, j)
            if d < 2.5:
                clashes += 1
    return {"min_dist_A": round(min_dist, 2) if min_dist is not None else None,
            "min_pair": list(pair) if pair else None,
            "n_pairs_under_2.5A": clashes}


def atropisomer_scan(mol: Any) -> Dict[str, Any]:
    """Rotatable aromatic C-C bonds with ortho blocking (fused bonds excluded)."""
    from rdkit import Chem

    axes = []
    ri = mol.GetRingInfo()
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if not (a1.GetIsAromatic() and a2.GetIsAromatic()
                and a1.GetSymbol() == "C" and a2.GetSymbol() == "C"
                and b.GetBondType() == Chem.BondType.SINGLE):
            continue
        if ri.AreAtomsInSameRing(a1.GetIdx(), a2.GetIdx()):
            continue  # fused biaryl bond is not rotatable

        def ortho_blocked(a: Any, other: Any) -> int:
            return sum(1 for nb in a.GetNeighbors()
                       if nb.GetIdx() != other.GetIdx() and nb.GetAtomicNum() > 1)

        blocked = ortho_blocked(a1, a2) + ortho_blocked(a2, a1)
        axes.append({"bond": (a1.GetIdx(), a2.GetIdx()), "ortho_blocked": blocked})
    axes.sort(key=lambda ax: -ax["ortho_blocked"])
    return {"n_biaryl_axes": len(axes),
            "max_ortho_blocked": axes[0]["ortho_blocked"] if axes else 0,
            "axes": axes}


# --------------------------------------------------------------------------- #
# per-molecule pipeline (runs in worker processes)
# --------------------------------------------------------------------------- #
def process_molecule(entry: Dict[str, str], num_confs: int, out_dir: str) -> Dict[str, Any]:
    """Full-stack pipeline for one molecule; never raises (records failures)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

    t0 = time.time()
    mid = entry["id"]
    rec: Dict[str, Any] = {"id": mid, "name": entry["name"], "name_cn": entry["name_cn"],
                           "smiles": entry["smiles"], "challenge": entry["challenge"],
                           "status": "ok", "risk_flags": [], "notes": []}
    try:
        # ---- 1. parse / sanitize / stereo / fragments --------------------- #
        mol_full, parse_status = robust_parse(entry["smiles"])
        rec["parse_status"] = parse_status
        if mol_full is None:
            rec["status"] = "failed_parse"
            rec["notes"].append(f"PARSE FAILURE: {parse_status}")
            return rec
        Chem.AssignStereochemistry(mol_full, cleanIt=True, force=True)

        try:
            frags = Chem.GetMolFrags(mol_full, asMols=True)
        except Exception:
            frags = [mol_full]
        primary = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        rec["n_fragments"] = len(frags)
        if len(frags) > 1:
            rec["notes"].append(
                f"input is a {len(frags)}-fragment disconnected assembly; bulk properties "
                f"reported for the full assembly, graph/3D on the largest fragment "
                f"({primary.GetNumHeavyAtoms()} heavy atoms)")
            rec["risk_flags"].append("FRAG")

        # ---- 2. 2D graph verification ------------------------------------- #
        centers = Chem.FindMolChiralCenters(mol_full, includeUnassigned=True,
                                            useLegacyImplementation=False)
        n_assigned = sum(1 for _, tag in centers if tag in ("R", "S"))
        n_unassigned = len(centers) - n_assigned
        ri = primary.GetRingInfo()
        ring_sizes = sorted(len(r) for r in ri.AtomRings())
        rec["graph"] = {
            "heavy_atoms": primary.GetNumHeavyAtoms(),
            "bonds": primary.GetNumBonds(),
            "aromatic_atoms": sum(1 for a in primary.GetAtoms() if a.GetIsAromatic()),
            "max_ring": ring_sizes[-1] if ring_sizes else 0,
            "n_rings": len(ring_sizes),
            "n_3_4_rings": sum(1 for s in ring_sizes if s <= 4),
            "n_macrocycles": sum(1 for s in ring_sizes if s >= 12),
        }
        rec["stereo"] = {"assigned": n_assigned, "unassigned": n_unassigned,
                         "labels": [tag for _, tag in centers if tag != "?"][:12]}
        if n_unassigned > 0:
            rec["risk_flags"].append("UST")

        scaffold = murcko_analysis(primary)
        rec["scaffold"] = scaffold

        feats = build_gnn_features(primary)
        npz_path = Path(out_dir) / "features" / f"{mid}.npz"
        import numpy as np
        np.savez_compressed(npz_path, feature_names=np.array(ATOM_FEATURES),
                            bond_feature_names=np.array(BOND_FEATURES), **feats)
        rec["gnn"] = {
            "n_nodes": int(feats["atom_feat"].shape[0]),
            "n_directed_edges": int(feats["edge_index"].shape[1]),
            "atom_feature_dim": int(feats["atom_feat"].shape[1]),
            "features_file": str(npz_path),
            "smoke": gnn_readiness_smoke(feats["atom_feat"], feats["edge_index"], feats["bond_feat"]),
        }

        # ---- 3. physicochemical profile (full input as written) ----------- #
        n_f = sum(1 for a in mol_full.GetAtoms() if a.GetSymbol() == "F")
        n_pos = sum(1 for a in mol_full.GetAtoms() if a.GetFormalCharge() > 0)
        n_neg = sum(1 for a in mol_full.GetAtoms() if a.GetFormalCharge() < 0)
        rec["props"] = {
            "mw_full_input": round(Descriptors.MolWt(mol_full), 2),
            "mw_primary": round(Descriptors.MolWt(primary), 2),
            "clogp": round(Crippen.MolLogP(mol_full), 2),
            "tpsa": round(rdMolDescriptors.CalcTPSA(mol_full), 2),
            "fsp3": round(Lipinski.FractionCSP3(mol_full), 3),
            "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol_full)),
            "hbd": int(Lipinski.NumHDonors(mol_full)),
            "hba": int(Lipinski.NumHAcceptors(mol_full)),
            "n_f": n_f,
            "formal_charge": int(Chem.GetFormalCharge(mol_full)),
            "n_pos_atoms": n_pos, "n_neg_atoms": n_neg,
        }

        # ---- 4. 3D ensemble ----------------------------------------------- #
        ens = embed_and_optimize(primary, num_confs)
        if not ens.get("ok"):
            rec["status"] = "partial_2d_only"
            rec["notes"].append(f"3D FAILURE: {ens.get('reason')}")
            rec["conformers"] = {"status": "failed", "n_requested": num_confs}
        else:
            mol_h = ens["mol_h"]
            conf = mol_h.GetConformer(ens["e_min_conf_id"])
            imhb = imhb_analysis(mol_h, conf)
            steric = steric_scan(mol_h, conf)
            rec["conformers"] = {
                "status": "ok", "ff": ens["ff"],
                "n_requested": num_confs, "n_embedded": ens["n_embedded"],
                "not_converged": ens["not_converged"],
                "e_min": round(ens["e_min"], 2),
                "e_min_conf_id": ens["e_min_conf_id"],
                "e_max": round(ens["e_max"], 2),
                "delta_e": round(ens["delta_e"], 2),
                "used_random_coords": ens["used_random_coords"],
                "imhb": imhb, "steric": steric,
            }

            sdf_dir = Path(out_dir) / "sdf"
            title = f"{mid} | {entry['name']}"
            ens_writer = Chem.SDWriter(str(sdf_dir / f"{mid}_ensemble.sdf"))
            try:
                mol_h.SetProp("_Name", title)
                for cid, energy in zip(ens["conf_ids"], ens["energies"]):
                    mol_h.SetProp("E_kcal_per_mol", f"{energy:.4f}")
                    mol_h.SetProp("conf_id", str(cid))
                    ens_writer.write(mol_h, confId=cid)
            finally:
                ens_writer.close()

            min_mol = Chem.Mol(mol_h)
            min_mol.RemoveAllConformers()
            min_mol.AddConformer(Chem.Conformer(conf), assignId=True)
            min_mol.SetProp("_Name", f"{title} | EMIN")
            min_mol.SetProp("E_min_kcal_per_mol", f"{ens['e_min']:.4f}")
            min_mol.SetProp("force_field", ens["ff"])
            min_mol.SetProp("delta_E_ensemble", f"{ens['delta_e']:.4f}")
            min_mol.SetProp("n_IMHB", str(imhb["count"]))
            writer = Chem.SDWriter(str(sdf_dir / f"{mid}_min.sdf"))
            try:
                writer.write(min_mol)
            finally:
                writer.close()

        # ---- 5. failure-analysis flags ------------------------------------- #
        g, p, c = rec["graph"], rec["props"], rec.get("conformers", {})
        if g["max_ring"] >= 12 or g["n_macrocycles"] > 0:
            rec["risk_flags"].append("MAC")
        if g["n_3_4_rings"] > 0:
            rec["risk_flags"].append("STRAIN")
        rec["atropisomerism"] = atropisomer_scan(primary)
        if rec["atropisomerism"]["max_ortho_blocked"] >= 3:
            rec["risk_flags"].append("ATRO")
        if p["n_f"] >= 6:
            rec["risk_flags"].append("FLUORO")
        if (p["n_pos_atoms"] > 0 and p["n_neg_atoms"] > 0) or p["formal_charge"] != 0:
            rec["risk_flags"].append("ZWIT")
        if p["mw_full_input"] > 800:
            rec["risk_flags"].append("BIG")
        if p["rotatable_bonds"] >= 15:
            rec["risk_flags"].append("FLEX")
        if c.get("status") == "ok":
            if c["delta_e"] >= 10:
                rec["risk_flags"].append("ENTROPY")
            if c.get("steric", {}).get("min_dist_A") is not None \
                    and c["steric"]["min_dist_A"] < 2.5:
                rec["risk_flags"].append("CLASH")
            if c.get("ff") == "UFF":
                rec["risk_flags"].append("UFF")

        # keep order stable for the report
        order = list(FLAG_LEGEND.keys())
        rec["risk_flags"] = sorted(set(rec["risk_flags"]), key=order.index)

    except Exception as exc:
        rec["status"] = "error"
        rec["notes"].append(f"EXCEPTION: {type(exc).__name__}: {exc}")
        rec["_traceback"] = traceback.format_exc()
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


# --------------------------------------------------------------------------- #
# report generation
# --------------------------------------------------------------------------- #
def _fmt(x: Optional[float], nd: int = 1) -> str:
    return "-" if x is None else f"{x:.{nd}f}"


def build_markdown(records: List[Dict[str, Any]], meta: Dict[str, Any]) -> str:
    """Assemble the full benchmark_report.md content."""
    lines: List[str] = []
    lines.append("# Hardcore Benchmark: 10 Novel Complex Molecules")
    lines.append("")
    lines.append(f"*Generated {meta['timestamp']} | RDKit {meta['rdkit_version']} | "
                 f"{meta['num_confs']} conformers/molecule, ETKDGv3 + MMFF94/UFF | "
                 f"wall time {meta['wall_time_s']} s*")
    lines.append("")
    lines.append("## Benchmark Table")
    lines.append("")
    header = ("| # | Molecule | MW (Da) | cLogP | TPSA (Å²) | Fsp³ | RotB | Stereo a/u | "
              "MaxRing | Conf | FF | E_min | ΔE_ens | IMHB | d_min (Å) | Risk |")
    lines.append(header)
    lines.append("|" + "---|" * 16)
    for r in records:
        p, g = r.get("props", {}), r.get("graph", {})
        c = r.get("conformers", {})
        st = r.get("stereo", {})
        if r["status"] == "failed_parse":
            row = (f"| {r['id']} | {r['name']} | parse FAILED: see analysis |" +
                   " - |" * 13 + " |")
        else:
            imhb_n = c.get("imhb", {}).get("count", "-") if c.get("status") == "ok" else "-"
            dmin = c.get("steric", {}).get("min_dist_A") if c.get("status") == "ok" else None
            row = (f"| {r['id']} | {r['name']} | {_fmt(p.get('mw_full_input'))} | "
                   f"{_fmt(p.get('clogp'), 2)} | {_fmt(p.get('tpsa'))} | "
                   f"{_fmt(p.get('fsp3'), 3)} | {p.get('rotatable_bonds', '-')} | "
                   f"{st.get('assigned', '-')}/{st.get('unassigned', '-')} | "
                   f"{g.get('max_ring', '-')} | "
                   f"{c.get('n_embedded', '-') if c.get('status') == 'ok' else 'fail'} | "
                   f"{c.get('ff', '-') if c.get('status') == 'ok' else '-'} | "
                   f"{_fmt(c.get('e_min')) if c.get('status') == 'ok' else '-'} | "
                   f"{_fmt(c.get('delta_e')) if c.get('status') == 'ok' else '-'} | "
                   f"{imhb_n} | {_fmt(dmin, 2)} | "
                   f"{', '.join(r.get('risk_flags', [])) or '-'} |")
        lines.append(row)
    lines.append("")
    lines.append("*E_min and ΔE_ens: total MMFF94/UFF steric energy in kcal/mol of the "
                 "lowest-energy conformer and the ensemble spread (E_max − E_min). "
                 "Energies are relative force-field quantities, not formation enthalpies. "
                 "Stereo a/u = assigned/unassigned chiral centers; d_min = shortest "
                 "non-bonded heavy-atom distance (≥4 bonds apart) in the E_min conformer.*")
    lines.append("")
    lines.append("**Risk flag legend:** " + "; ".join(f"`{k}` {v}" for k, v in FLAG_LEGEND.items()))
    lines.append("")
    lines.append("## GNN Feature Tensor Verification (PyG-ready)")
    lines.append("")
    lines.append("| # | Nodes | Directed edges | Atom feat dim | Feature file | GCN smoke |")
    lines.append("|---|---|---|---|---|---|")
    for r in records:
        gnn = r.get("gnn")
        if not gnn:
            lines.append(f"| {r['id']} | - | - | - | - | n/a (parse failed) |")
            continue
        smoke = gnn.get("smoke", {})
        smoke_txt = (f"OK (torch {smoke.get('torch')}, pyg {smoke.get('pyg')}, "
                     f"out {smoke.get('embedding_shape')})" if smoke.get("ok")
                     else f"unavailable ({smoke.get('reason')})")
        lines.append(f"| {r['id']} | {gnn['n_nodes']} | {gnn['n_directed_edges']} | "
                     f"{gnn['atom_feature_dim']} | `{gnn['features_file']}` | {smoke_txt} |")
    lines.append("")
    lines.append("## Per-Molecule Failure Analysis")
    for r in records:
        lines.append("")
        lines.append(f"### {r['id']} — {r['name']}（{r.get('name_cn', '')}）")
        lines.append(f"- **SMILES:** `{r['smiles']}`")
        if r.get("scaffold"):
            lines.append(f"- **Bemis-Murcko scaffold (generic):** `{r['scaffold']['generic_smiles']}` "
                         f"(scaffold atom fraction {r['scaffold']['scaffold_atom_fraction']})")
        lines.append(f"- **Challenge:** {r['challenge']}")
        if r.get("parse_status"):
            lines.append(f"- **Parse status:** {r['parse_status']}")
        for note in r.get("notes", []):
            lines.append(f"- **Note:** {note}")
        if r.get("atropisomerism", {}).get("n_biaryl_axes"):
            ax = r["atropisomerism"]
            lines.append(f"- **Atropisomer scan:** {ax['n_biaryl_axes']} rotatable biaryl "
                         f"axes, max ortho-blocking {ax['max_ortho_blocked']}/4")
        if r.get("conformers", {}).get("status") == "ok":
            c = r["conformers"]
            hbs = "; ".join(f"D{b['donor']}-H{b['h']}···A{b['acceptor']} "
                            f"({b['dHA_A']} Å, {b['angle_deg']}°)"
                            for b in c["imhb"]["bonds"][:6]) or "none detected"
            lines.append(f"- **IMHB in E_min conformer:** {c['imhb']['count']} — {hbs}")
            lines.append(f"- **Ensemble:** {c['n_embedded']} confs optimized with {c['ff']} "
                         f"({c['not_converged']} not fully converged); "
                         f"E ∈ [{c['e_min']}, {c['e_max']}] kcal/mol")
        flags = r.get("risk_flags", [])
        if flags:
            lines.append("- **2D-GNN failure modes:** " +
                         "; ".join(f"`{f}` {FLAG_LEGEND[f]}" for f in flags))
        lines.append(f"- **Runtime:** {r.get('seconds', '-')} s")
    lines.append("")
    lines.append("## Artifacts")
    lines.append("- `sdf/<ID>_ensemble.sdf` — full optimized ensemble (E_kcal_per_mol per conf)")
    lines.append("- `sdf/<ID>_min.sdf` — lowest-energy conformer (ΔE = 0 reference)")
    lines.append("- `features/<ID>.npz` — atom/bond feature tensors, PyG-loadable")
    lines.append("- `benchmark_results.json` — complete machine-readable records")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# shutdown (fault-tolerant execution protocol)
# --------------------------------------------------------------------------- #
def schedule_shutdown(delay_s: int) -> Dict[str, Any]:
    """Cancellable clean power-off (Windows /s /t, POSIX shutdown -h +m)."""
    info: Dict[str, Any] = {"scheduled": False, "delay_seconds": int(delay_s)}
    minutes = max(1, int(round(delay_s / 60)))
    try:
        if sys.platform.startswith("win"):
            cmd = ["shutdown", "/s", "/t", str(int(delay_s)),
                   "/c", "molecule_benchmark.py finished -- powering off."]
            info["cancel_command"] = "shutdown /a"
        else:
            cmd = ["shutdown", "-h", f"+{minutes}"]
            info["cancel_command"] = "sudo shutdown -c"
        proc = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=30)
        info["command"] = " ".join(cmd)
        info["scheduled"] = proc.returncode == 0
        if not info["scheduled"] and not sys.platform.startswith("win"):
            proc2 = subprocess.run(["sudo", "-n"] + cmd, capture_output=True, text=True,
                                   errors="replace", timeout=30)
            info["scheduled"] = proc2.returncode == 0
        elif not info["scheduled"]:
            info["error"] = (proc.stderr or proc.stdout or "").strip()
    except Exception as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"
    if info["scheduled"]:
        print("\n" + "!" * 62)
        print(f"  SYSTEM SHUTDOWN SCHEDULED -- power-off in ~{delay_s} s")
        print(f"  cancel with : {info.get('cancel_command')}")
        print("!" * 62)
    return info


def write_json_atomic(path: Path, payload: Dict[str, Any]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str),
                   encoding="utf-8")
    os.replace(tmp, path)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="10-molecule hardcore RDKit benchmark "
                                            "(2D graph + 3D ensemble + profiling + failure analysis)")
    p.add_argument("--out_dir", default="bench_results")
    p.add_argument("--conformers", type=int, default=50)
    p.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    p.add_argument("--only", default=None, help="comma-separated molecule ids, e.g. M01,M03")
    p.add_argument("--auto_shutdown", action="store_true",
                   help="power off --shutdown_delay seconds after a completed run")
    p.add_argument("--shutdown_delay", type=int, default=60)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    (out_dir / "sdf").mkdir(parents=True, exist_ok=True)
    (out_dir / "features").mkdir(parents=True, exist_ok=True)

    import rdkit
    registry = REGISTRY
    if args.only:
        wanted = {s.strip().upper() for s in args.only.split(",")}
        registry = [e for e in registry if e["id"] in wanted]

    meta: Dict[str, Any] = {
        "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
        "rdkit_version": rdkit.__version__,
        "num_confs": args.conformers,
        "wall_time_s": None,
    }
    records: List[Dict[str, Any]] = []
    exit_code = 0
    t0 = time.time()
    try:
        print(f"[bench] RDKit {rdkit.__version__} | {len(registry)} molecules | "
              f"{args.conformers} confs each | workers={args.workers}")
        if args.workers > 1 and len(registry) > 1:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures = {pool.submit(process_molecule, e, args.conformers, str(out_dir)): e
                           for e in registry}
                for fut in as_completed(futures):
                    entry = futures[fut]
                    try:
                        rec = fut.result()
                    except Exception as exc:
                        rec = {"id": entry["id"], "name": entry["name"],
                               "status": "error",
                               "notes": [f"worker crashed: {type(exc).__name__}: {exc}"]}
                    records.append(rec)
                    n3d = rec.get("conformers", {}).get("n_embedded", "-")
                    print(f"[bench] {rec['id']} done: status={rec['status']} "
                          f"3d_confs={n3d} flags={','.join(rec.get('risk_flags', [])) or '-'} "
                          f"({rec.get('seconds', '?')} s)")
        else:
            for e in registry:
                rec = process_molecule(e, args.conformers, str(out_dir))
                records.append(rec)
                print(f"[bench] {rec['id']} done: status={rec['status']} "
                      f"({rec.get('seconds', '?')} s)")
        records.sort(key=lambda r: r["id"])
    except Exception as exc:
        print(f"[bench][FATAL] {type(exc).__name__}: {exc}")
        traceback.print_exc()
        exit_code = 1
    finally:
        meta["wall_time_s"] = round(time.time() - t0, 1)
        payload = {"meta": meta, "flag_legend": FLAG_LEGEND, "results": records}
        try:
            write_json_atomic(out_dir / "benchmark_results.json", payload)
            report = build_markdown(records, meta)
            (out_dir / "benchmark_report.md").write_text(report, encoding="utf-8")
            print(f"[bench] artifacts: {out_dir / 'benchmark_report.md'} | "
                  f"{out_dir / 'benchmark_results.json'}")
            print("\n" + report)
        except Exception as exc:
            print(f"[bench][FATAL] report generation failed: {exc}")
            traceback.print_exc()
            exit_code = 1
        if args.auto_shutdown:
            payload["shutdown"] = schedule_shutdown(args.shutdown_delay)
            write_json_atomic(out_dir / "benchmark_results.json", payload)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
