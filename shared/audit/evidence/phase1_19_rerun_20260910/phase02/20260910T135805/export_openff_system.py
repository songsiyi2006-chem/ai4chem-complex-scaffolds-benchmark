#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
export_openff_system.py
=======================
Serialize an OpenFF-Sage-parameterized OpenMM System for one molecule.

Runs inside the dedicated `phase2ff` conda environment (openff-toolkit +
openff-interchange + openmm + rdkit) because the OpenFF PyPI stack is
currently in a yanked/migration state on the main runtime interpreter.

Usage:
    python export_openff_system.py <smiles> <output_prefix>

Writes:
    <prefix>_system.xml   OpenMM System (Sage 2.1.0, h-bond constrained)
    <prefix>_start.pdb    MMFF-minimized starting geometry (with elements)
    <prefix>_meta.json    particle/element/force inventory
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import openmm
from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule
from rdkit import Chem
from rdkit.Chem import AllChem


def main() -> int:
    smiles, prefix = sys.argv[1], sys.argv[2]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f"RDKit parse failure: {smiles}")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    mh = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 0x2026
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    cid = AllChem.EmbedMolecule(mh, params)
    if cid < 0:
        retry = AllChem.ETKDGv3()
        retry.useRandomCoords = True
        retry.randomSeed = 7
        cid = AllChem.EmbedMolecule(mh, retry)
    if cid < 0:
        raise SystemExit("ETKDGv3 embedding failed")
    if AllChem.MMFFHasAllMoleculeParams(mh):
        AllChem.MMFFOptimizeMolecule(mh, maxIters=2000)
    else:
        AllChem.UFFOptimizeMolecule(mh, maxIters=2000)

    offmol = Molecule.from_rdkit(mh, allow_undefined_stereo=False)
    sage = ForceField("openff-2.1.0.offxml")  # h-bond constrained -> 2 fs stable

    # Charge strategy: if ANY backend can compute AM1-BCC (NAGL/OpenEye/sqc),
    # keep Sage electrostatics intact. Otherwise deregister Electrostatics,
    # let Sage emit valence+vdW only, and assemble MMFF94 charges (RDKit,
    # same atom order) into the NonbondedForce ourselves, applying the
    # SMIRNOFF 1-4 coulomb scale (1/1.2) to every exception.
    try:
        Molecule.from_smiles("CCO").assign_partial_charges("am1bcc")
        charge_source = "am1bcc"
    except Exception:
        for handler_name in ("Electrostatics", "ToolkitAM1BCC",
                             "ChargeIncrementModel", "LibraryCharges"):
            try:
                sage.deregister_parameter_handler(handler_name)
            except Exception:
                pass
        charge_source = "mmff94 (RDKit; Sage valence+vdW; 1-4 coulomb 1/1.2)"

    interchange = Interchange.from_smirnoff(sage, [offmol])
    try:
        interchange.box = None
    except Exception:
        pass

    try:
        system = interchange.to_openmm_system()
    except AttributeError:
        system = interchange.to_openmm(combine_nonbonded_forces=True)[0]

    if charge_source.startswith("mmff94"):
        props = AllChem.MMFFGetMoleculeProperties(mh)
        charges = [props.GetMMFFPartialCharge(i) for i in range(mh.GetNumAtoms())]
        net = sum(charges)
        formal = Chem.GetFormalCharge(mh)
        if abs(net - formal) > 0.02:
            raise SystemExit(f"MMFF94 charge sanity failed: sum={net:.3f} vs formal={formal}")
        nonbonded = [f for f in system.getForces()
                     if isinstance(f, openmm.NonbondedForce)]
        if len(nonbonded) != 1:
            raise SystemExit(f"expected exactly one NonbondedForce, got {len(nonbonded)}")
        nb = nonbonded[0]
        for i in range(system.getNumParticles()):
            _, sigma, epsilon = nb.getParticleParameters(i)
            nb.setParticleParameters(i, charges[i], sigma, epsilon)
        for k in range(nb.getNumExceptions()):
            i, j, qprod, sigma, epsilon = nb.getExceptionParameters(k)
            nb.setExceptionParameters(k, i, j, charges[i] * charges[j] / 1.2,
                                      sigma, epsilon)

    with open(f"{prefix}_system.xml", "w") as fh:
        fh.write(openmm.XmlSerializer.serialize(system))

    offmol.to_file(f"{prefix}_start.pdb", file_format="pdb")

    meta = {
        "smiles": smiles,
        "charge_source": charge_source,
        "n_particles": system.getNumParticles(),
        "n_constraints": system.getNumConstraints(),
        "forces": [system.getForce(i).getName() for i in range(system.getNumForces())],
        "atomic_numbers": [int(a.atomic_number) for a in offmol.atoms],
        "openff_toolkit": __import__("openff.toolkit", fromlist=["__version__"]).__version__,
    }
    Path(f"{prefix}_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"EXPORT_OK particles={meta['n_particles']} constraints={meta['n_constraints']} "
          f"forces={meta['forces']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
