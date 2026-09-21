"""Bounded ETKDGv3/MMFF94s benchmark, NOT a molecular-glue discovery engine."""
from pathlib import Path
import math


def benchmark_smiles():
    # Simple chiral lactam software fixtures, not analogues of macroPROTAC-1.
    return ["O=C1N[C@@H](C)" + "C" * n + "C1" for n in (8, 9, 10)]


def generate(out, smiles=None, seed=31, conformers=6):
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem, rdMolTransforms
    if not 1 <= conformers <= 12:
        raise ValueError("local conformer count must be 1..12")
    smiles = benchmark_smiles() if smiles is None else list(smiles)
    if not 1 <= len(smiles) <= 100:
        raise ValueError("1..100 explicit molecules required")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    seen, records = set(), []
    for index, smi in enumerate(smiles):
        base = Chem.MolFromSmiles(smi)
        if base is None or base.GetNumHeavyAtoms() > 120:
            raise ValueError("invalid SMILES or local size cap exceeded")
        canonical = Chem.MolToSmiles(base, isomericSmiles=True)
        if canonical in seen:
            raise ValueError("duplicate stereochemical identity")
        seen.add(canonical)
        rings = list(base.GetRingInfo().AtomRings())
        if not rings or max(map(len, rings)) < 12:
            raise ValueError("requires a ring of at least 12 atoms")
        if not Chem.FindMolChiralCenters(base, includeUnassigned=False):
            raise ValueError("an explicitly specified stereocentre is required")
        mol = Chem.AddHs(base)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed + index
        params.numThreads = 1
        params.enforceChirality = True
        params.useMacrocycleTorsions = True
        params.useMacrocycle14config = True
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=conformers, params=params))
        if not ids or not AllChem.MMFFHasAllMoleculeParams(mol):
            raise RuntimeError("embedding failed or MMFF parameters missing")
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        entries = []
        ring = max(rings, key=len)
        for cid in ids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
            status = ff.Minimize(maxIts=700)
            energy = float(ff.CalcEnergy())
            if not math.isfinite(energy):
                raise RuntimeError("nonfinite MMFF energy")
            entries.append(dict(conformer=int(cid), converged=status == 0,
                                energy_kcal_mol=energy,
                                ring_dihedral_deg=float(rdMolTransforms.GetDihedralDeg(
                                    mol.GetConformer(cid), *ring[:4]))))
        converged = [e['energy_kcal_mol'] for e in entries if e['converged']]
        for e in entries:
            e['relative_energy_kcal_mol'] = e['energy_kcal_mol'] - min(converged) if e['converged'] and converged else None
        path = out / f"fixture_{index + 1}.sdf"
        writer = Chem.SDWriter(str(path))
        for e in entries:
            mol.SetProp("_Name", f"fixture_{index + 1}_conf_{e['conformer']}")
            mol.SetProp("evidence_class", "LOCAL_MMFF_BENCHMARK_NOT_GLUE")
            mol.SetProp("MMFF94s_energy_kcal_mol", str(e['energy_kcal_mol']))
            mol.SetProp("MMFF_converged", str(e['converged']))
            writer.write(mol, confId=e['conformer'])
        writer.close()
        records.append(dict(id=f"fixture_{index+1}", smiles=canonical, ring_size=len(ring),
                            sdf=path.name, conformers=entries))
    return dict(evidence="EXECUTED_MMFF_BENCHMARK", rdkit_version=rdBase.rdkitVersion,
                seed=seed, records=records, target_binding_evaluated=False,
                absolute_ring_strain_computed=False, proposed_leads=[])
