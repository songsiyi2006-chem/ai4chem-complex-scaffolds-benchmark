"""Focused failure-path regressions; no molecular engines or dynamics run."""
import ast
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import types
import unittest
from unittest.mock import Mock, patch

import numpy as np

ROOT = Path(__file__).resolve().parent


def functions(phase, names, **extra):
    if phase == 5 and any(n in names for n in ('run_xtb', 'xtb_opt')):
        names = list(dict.fromkeys([*names, 'XTBFailure']))
    source = next(ROOT.glob(f'run_phase{phase}_*.py'))
    tree = ast.parse(source.read_text(encoding='utf-8'))
    nodes = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.name in names]
    assert len(nodes) == len(names)
    ns = dict(np=np, math=math, Path=Path, re=re, os=os, shutil=shutil,
              tempfile=tempfile, subprocess=subprocess, XTB_EXE='mock-xtb',
              EH_EV=27.211386245988, T_REF=298.15)
    ns.update(extra)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), 'exec'), ns)
    return ns


def probe_phase5_restraint(cache_path):
    """Explicit opt-in diagnostic: two single points on ONE fresh saved geometry.

    Not called by unittest. No optimization, no kinetics or TS validation.
    """
    import hashlib
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    cache = Path(cache_path)
    rec = json.loads(cache.read_text())
    exe = Path(sys.prefix) / 'Library/bin/xtb.exe'
    ns = functions(5, ['_xyz_text', 'xtb_env', 'run_xtb', 'build_R', 'gate', 'formula_of'],
                   XTB_EXE=str(exe), XTB_DIR=str(exe.parent), Chem=Chem,
                   rdMolDescriptors=rdMolDescriptors)
    molecule = ns['build_R']()
    labels = {a.GetProp('p5label'): a.GetIdx()+1 for a in molecule.GetAtoms() if a.HasProp('p5label')}
    outputs = {}
    for constrained in (False, True):
        with tempfile.TemporaryDirectory(prefix='p5restraint_probe_') as td:
            args = ['--sp']
            if constrained:
                inp = ('$constrain\n  force constant=0.8\n'
                       f'  distance: {labels["N1"]}, {labels["C2"]}, {rec["d_star"]:.4f}\n$end\n')
                (Path(td)/'constr.inp').write_text(inp)
                args += ['--input', 'constr.inp']
            stdout, _ = ns['run_xtb'](rec['numbers'], rec['positions'], args,
                                      chrg=rec['chrg'], workdir=Path(td))
            matches = re.findall(r'TOTAL ENERGY\s+(-?\d+\.\d+)', stdout)
            if not matches:
                raise RuntimeError('Single point energy missing')
            outputs['constrained' if constrained else 'unconstrained'] = {
                'E_eh': float(matches[-1]),
                'energy_lines': [s.strip() for s in stdout.splitlines()
                                 if 'energy' in s.lower() and any(k in s.lower() for k in ('total', 'constrain', 'restraint'))]}
    outputs.update(geometry_path=str(cache.resolve()),
                   geometry_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),
                   charge=rec['chrg'], gfn=2, solvent='none', optimized=False,
                   geometry_scope='saved refined TS1; original scan-frame positions not persisted',
                   constraint_atoms_1based=[labels['N1'], labels['C2']],
                   constraint_distance_A=rec['d_star'], force_constant=.8)
    outputs['difference_kcal_mol'] = (outputs['constrained']['E_eh']-outputs['unconstrained']['E_eh'])*627.5094740631
    print(json.dumps(outputs, indent=2))
    return outputs


class RerunTests(unittest.TestCase):
    def test_phase5_cache_keys_distinguish_faces_and_unlocked_state(self):
        tree = ast.parse((ROOT/'run_phase5_chemical_world_model.py').read_text(encoding='utf-8'))
        constants = {n.targets[0].id: ast.literal_eval(n.value) for n in tree.body
                     if isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name)
                     and n.targets[0].id in ('FACE_CACHE_KEYS', 'UNLOCKED_CACHE_KEY')}
        faces = constants['FACE_CACHE_KEYS']
        keys = [faces['M'], faces['m'], constants['UNLOCKED_CACHE_KEY'], 'cx_I1Cat']
        self.assertEqual(len(set(k.casefold() for k in keys)), 4)
        with tempfile.TemporaryDirectory() as td:
            ns = functions(5, ['cached_or_compute'], CACHE=Path(td), json=json,
                           _log=lambda *a: None, _warn=lambda *a: None)
            computes = [Mock(return_value={'state': i}) for i in range(4)]
            for i, key in enumerate(keys):
                self.assertEqual(ns['cached_or_compute'](key, computes[i]), {'state': i})
            for i, key in enumerate(keys):
                self.assertEqual(ns['cached_or_compute'](key, computes[i]), {'state': i})
                computes[i].assert_called_once()
        # Call sites and hydration must use the shared keys, not stale literals.
        stale = {'cx_I1Cat_M', 'cx_I1Cat_m', 'cx_TS2aM', 'cx_TS2am'}
        self.assertFalse(any(isinstance(n, ast.Constant) and isinstance(n.value, str)
                             and n.value in stale for n in ast.walk(tree)))

    def test_deprotonation_pose_uses_global_oxygen_preserves_fragments(self):
        ns = functions(5, ['prepare_deprotonation_pose', 'kabsch_rotate'])
        nums = np.array([6, 1, 6, 8, 15, 8])
        pos = np.array([[0., 0, 0], [1.1, 0, 0], [0, 1.4, 0],
                        [0, 8, 0], [0, 9.5, 0], [1, 9.5, 0]])
        new = ns['prepare_deprotonation_pose'](nums, pos, 0, 1, 3, 4, 3)
        np.testing.assert_array_equal(new[:3], pos[:3])
        np.testing.assert_allclose(new[3], [2.75, 0, 0], atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(new[3:, None]-new[None, 3:], axis=2),
                                   np.linalg.norm(pos[3:, None]-pos[None, 3:], axis=2))
        self.assertGreater(new[4, 0], new[3, 0])
        with self.assertRaisesRegex(ValueError, 'wrong fragments'):
            ns['prepare_deprotonation_pose'](nums, pos, 0, 1, 0, 4, 3)

    def test_scan_failure_retains_engine_evidence_without_weakening_gate(self):
        records = []
        ns = functions(5, ['scan_1d', 'XTBFailure'],
                       _save_scan_record=lambda label, index, record: records.append(record),
                       _warn=lambda *a: None, _log=lambda *a: None, EH_KCAL=627.5)
        ns['xtb_opt'] = Mock(side_effect=[
            ns['XTBFailure']('not converged', stdout='full iteration history',
                              files={'xtbopt.xyz': 'last geometry'}),
            (np.zeros((2, 3)), -1.), (np.ones((2, 3)), -2.)])
        with self.assertRaisesRegex(RuntimeError, 'too few converged frames \\(2\\)'):
            ns['scan_1d']([6, 1], np.zeros((2, 3)), lambda s: [(1, 2, s)], [1.6, 1.4, 1.1])
        self.assertEqual(records[0]['stdout'], 'full iteration history')
        self.assertEqual(records[0]['engine_files']['xtbopt.xyz'], 'last geometry')
        self.assertEqual([r['converged'] for r in records], [False, True, True])

    def test_spliced_hydrogen_constraint_uses_length_not_spring(self):
        from openmm import System, HarmonicBondForce, unit
        from rdkit import Chem
        ns = functions(3, ['_splice_bonds'])
        mol = Chem.AddHs(Chem.MolFromSmiles('C'))
        force = HarmonicBondForce()
        force.addBond(0, 1, .109*unit.nanometer,
                      300000*unit.kilojoule_per_mole/unit.nanometer**2)
        for existing, expected in ((set(), 1), ({frozenset((0, 1))}, 0)):
            system = System()
            for _ in range(mol.GetNumAtoms()):
                system.addParticle(1.)
            combined = HarmonicBondForce()
            ns['_splice_bonds'](force, combined, system, mol, 0, existing)
            _, _, length, spring = combined.getBondParameters(0)
            self.assertAlmostEqual(length.value_in_unit(unit.nanometer), .109)
            self.assertAlmostEqual(spring.value_in_unit(unit.kilojoule_per_mole/unit.nanometer**2), 300000)
            self.assertEqual(system.getNumConstraints(), expected)
            if expected:
                self.assertAlmostEqual(system.getConstraintParameters(0)[2].value_in_unit(unit.nanometer), .109)

    def test_constrained_optimization_returns_unrestrained_energy(self):
        run = Mock(side_effect=[('GEOMETRY OPTIMIZATION CONVERGED\ntotal energy: -1.0 Eh',
                                {'xtbopt.xyz': 'saved geometry'}),
                               ('TOTAL ENERGY -2.0', {})])
        ns = functions(5, ['xtb_opt'], run_xtb=run,
                       _parse_xyz=lambda _: ([1, 1], np.zeros((2, 3))))
        _, energy = ns['xtb_opt']([1, 1], np.zeros((2, 3)), chrg=1,
                                  constraints=[(1, 2, 1.)])
        self.assertEqual(energy, -2.)
        self.assertEqual(run.call_args_list[1].args[2], ['--sp'])
        self.assertEqual(run.call_args_list[1].kwargs['chrg'], 1)
        self.assertNotIn('workdir', run.call_args_list[1].kwargs)

    def test_image_cache_invalidates_and_returns_force_copy(self):
        class Engine:
            name = 'fake xTB'
            def __init__(self):
                self.calls = 0
            def _run(self, numbers, positions):
                self.calls += 1
                return 1., np.ones_like(positions)

        ns = functions(4, ['_PerAtomCalc'], XTBWrap=Engine, ANIWrap=type('ANI', (), {}))
        positions = np.zeros((1, 3))
        atoms = types.SimpleNamespace(get_atomic_numbers=lambda: np.array([1]),
                                      get_positions=lambda: positions.copy())
        engine = Engine()
        calc = ns['_PerAtomCalc'](engine)
        forces = calc.get_forces(atoms)
        forces[:] = 999
        self.assertEqual(calc.get_potential_energy(atoms), 1.)
        np.testing.assert_array_equal(calc.get_forces(atoms), np.ones((1, 3)))
        self.assertEqual(engine.calls, 1)
        positions[0, 0] += .001
        calc.get_forces(atoms)
        self.assertEqual(engine.calls, 2)

    def test_warm_start_still_runs_requested_ci_steps(self):
        images = [Mock() for _ in range(3)]
        for i, im in enumerate(images):
            im.get_potential_energy.return_value = float(i == 1)
            im.get_atomic_numbers.return_value = np.array([1])
            im.get_positions.return_value = np.zeros((1, 3))
        neb = Mock(climb=False)
        neb.get_forces.return_value = np.zeros((1, 3))
        optimizers = []

        def optimizer(*args, **kwargs):
            opt = Mock()
            optimizers.append(opt)
            return opt

        class Warm:
            def _eval(self, atoms):
                return None

        modules = {'ase': types.SimpleNamespace(Atoms=Mock()),
                   'ase.mep.neb': types.SimpleNamespace(NEB=Mock(return_value=neb)),
                   'ase.optimize': types.SimpleNamespace(BFGS=optimizer, FIRE=optimizer),
                   'ase.io': types.SimpleNamespace(read=Mock(return_value=images), write=Mock())}
        ns = functions(4, ['stage2_neb', 'neb_fmax'], json=json, time=time,
                       ANIWrap=Warm, make_engine=lambda _: types.SimpleNamespace(name='xTB'),
                       _PerAtomCalc=lambda x: x, _log=Mock(), _warn=Mock(),
                       _images_to_xyz=lambda _: '', RESULTS={}, EV_TO_KCAL=23.06)
        args = types.SimpleNamespace(engine='auto', spring=.1, fmax=.05, max_neb_steps=500)
        with tempfile.TemporaryDirectory() as td, patch.dict(sys.modules, modules):
            self.assertTrue(ns['stage2_neb'](Path(td), args, force=True))
        self.assertEqual(len(optimizers), 3)
        optimizers[-1].run.assert_called_once_with(fmax=.05, steps=500)
        self.assertTrue(neb.climb)

    def test_ani_hessian_keeps_imaginary_mode_and_separates_cache(self):
        at = Mock()
        vibration = Mock()
        vibration.get_frequencies.return_value = np.array([100j, 200+0j])
        vibration.get_vibrations.return_value.get_mode.return_value = np.zeros((1, 3))
        factory = Mock(return_value=vibration)
        modules = {'ase': types.SimpleNamespace(Atoms=Mock(return_value=at)),
                   'ase.vibrations': types.SimpleNamespace(Vibrations=factory)}
        ns = functions(4, ['_ani_hessian'], ANIWrap=Mock(), _PerAtomCalc=Mock())
        with patch.dict(sys.modules, modules):
            first = ns['_ani_hessian']([1], [[0, 0, 0]])
            ns['_ani_hessian']([1], [[0, 0, 0]])
        self.assertEqual(first['frequencies'], [-100., 200.])
        paths = [c.kwargs['name'] for c in factory.call_args_list]
        self.assertNotEqual(paths[0], paths[1])
        self.assertFalse(Path(paths[0]).parent.exists())

    def test_failed_xtb_not_overridden_by_normal_text(self):
        proc = types.SimpleNamespace(returncode=1, stdout=b'normal termination', stderr=b'failure')
        ns = functions(5, ['run_xtb'], _xyz_text=lambda *a: '0\n\n',
                       xtb_env=lambda: {}, subprocess=types.SimpleNamespace(run=Mock(return_value=proc)))
        with self.assertRaisesRegex(RuntimeError, 'rc=1'):
            ns['run_xtb']([], [], ['--sp'])

    def test_optimization_requires_convergence(self):
        ns = functions(5, ['xtb_opt'], run_xtb=Mock(return_value=('not converged', {'xtbopt.xyz': 'partial'})))
        with self.assertRaisesRegex(RuntimeError, 'did not converge'):
            ns['xtb_opt']([1], [[0, 0, 0]])

    def test_esp_requires_real_charges(self):
        for charges in ({}, {'charges': 'nan'}):
            ns = functions(5, ['xtb_sp_charges'], run_xtb=Mock(return_value=('TOTAL ENERGY -1.0', charges)))
            with self.assertRaisesRegex(RuntimeError, 'ESP unavailable'):
                ns['xtb_sp_charges']([1], [[0, 0, 0]])

    def test_hessian_requires_spectrum(self):
        ns = functions(5, ['xtb_hess'], run_xtb=Mock(return_value=('TOTAL ENERGY -1.0', {'charges': '0'})))
        with self.assertRaisesRegex(RuntimeError, 'energy/frequencies'):
            ns['xtb_hess']([1], [[0, 0, 0]])

    def test_g98_partial_mode_block(self):
        ns = functions(4, ['_parse_g98_modes'])
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'g98.out'
            path.write_text(' Frequencies -- -100.0 200.0\n Atom AN\n 1 6 1 0 0 0 1 0\n 2 1 -1 0 0 0 -1 0\n')
            modes = ns['_parse_g98_modes'](path)
        self.assertEqual([m['freq'] for m in modes], [-100., 200.])
        np.testing.assert_array_equal(modes[0]['disp'], [[1, 0, 0], [-1, 0, 0]])

    def test_g98_missing_modes_rejected(self):
        ns = functions(4, ['_parse_g98_modes'])
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'g98.out'
            path.write_text(' Frequencies -- -100.0\n Atom AN\n malformed\n')
            with self.assertRaisesRegex(RuntimeError, 'no displacement rows'):
                ns['_parse_g98_modes'](path)


if __name__ == '__main__':
    unittest.main()
