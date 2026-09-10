"""Small regression tests; these are not substitutes for full scientific runs."""
import ast
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parent

class CampaignTests(unittest.TestCase):
    def test_public_evidence_redacts_paths_and_marks_unavailable(self):
        import publish_phase1_19_evidence as pub
        result=pub.portable({'path':str(pub.CAMPAIGN/'phase01'), 'nan':float('nan'), 'finite':1.25})
        self.assertTrue(result['path'].startswith('$CAMPAIGN'))
        self.assertIsNone(result['nan'])
        self.assertEqual(result['finite'],1.25)

    def test_public_evidence_excludes_mutable_running_results(self):
        import contextlib
        import io
        import json
        import time
        import publish_phase1_19_evidence as pub
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); campaign=root/'runs'; stamp=time.strftime('%Y%m%dT%H%M%S')
            attempt=campaign/'phase02'/stamp
            (attempt/'results_phase2').mkdir(parents=True)
            status={'phase':2,'status':'running','started':stamp,'scope':'default/full attempt',
                    'source_sha256':{},'scientific_acceptance':'pending output review'}
            (attempt/'rerun_status.json').write_text(json.dumps(status),encoding='utf-8')
            (attempt/'results_phase2/phase2_results.json').write_text('{"incomplete":true}',encoding='utf-8')
            with patch.object(pub,'ROOT',root),patch.object(pub,'CAMPAIGN',campaign),patch.object(pub,'OUT',root/'public'),contextlib.redirect_stdout(io.StringIO()):
                pub.main()
            manifest=json.loads((root/'public/manifest.json').read_text(encoding='utf-8'))
            self.assertEqual(len(manifest['published_files']),1)
            self.assertTrue(manifest['published_files'][0]['path'].endswith('run_status.json'))
            self.assertFalse((root/'public/phase02'/stamp/'results_phase2').exists())

    def test_phase9_mass_and_price_units(self):
        # AST-extracted pure model: no sklearn import or workflow execution.
        import math
        import numpy as np
        tree=ast.parse((ROOT/'run_phase9_self_driving_lab_compiler.py').read_text(encoding='utf-8'))
        names={'PHASE5_ANCHOR','SURROGATE','R_GAS_KCAL','solvent_mixture_properties',
               'bubble_point_c','antoine_vp_mmhg','eyring_rate','true_objectives'}
        nodes=[]
        for n in tree.body:
            if isinstance(n,ast.FunctionDef) and n.name in names: nodes.append(n)
            if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in names for t in n.targets): nodes.append(n)
        ns=dict(np=np,math=math)
        exec(compile(ast.Module(body=nodes,type_ignores=[]),'<model>','exec'),ns)
        for temp in [25.,50.,90.]:
            r=ns['true_objectives'](np.array([temp,5.,10.,.7]),np.random.default_rng(0),noisy=False)
            self.assertLessEqual(r['x_cat']+r['x_bg'],1.+1e-12)
            self.assertAlmostEqual(r['cost_usd_per_mol']*(r['yield_pct']/100),24000.)
            self.assertAlmostEqual(8.*5.*ns['SURROGATE']['catalyst_stock_mM']/1000./20.*100.,5.)
        self.assertAlmostEqual(ns['SURROGATE']['cp_dcm'],1.2045,places=3)

    def test_phase9_real_failure_not_overwritten_by_mock(self):
        tree=ast.parse((ROOT/'run_phase9_self_driving_lab_compiler.py').read_text(encoding='utf-8'))
        func=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='validate_ot2_protocol')
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            protocol=root/'protocol.py'
            protocol.write_text("metadata={'apiLevel':'2.15'}\ndef run(ctx):\n    pass\n",encoding='utf-8')
            sim=root/'sim.exe';sim.touch()
            import py_compile
            ns=dict(Path=Path,ROOT=root,ast=ast,os=os,sys=sys,py_compile=py_compile,
                    _install_opentrons_mock=lambda:None,_make_mock_ctx=lambda:None)
            exec(compile(ast.Module(body=[func],type_ignores=[]),'<extracted-source>','exec'),ns)
            failure=subprocess.CompletedProcess([],1,'','Error: invalid protocol')
            with patch.dict(os.environ,{'OPENTRONS_SIMULATE_EXE':str(sim)}),patch('subprocess.run',return_value=failure):
                result=ns['validate_ot2_protocol'](protocol)
            self.assertTrue(result['opentrons_simulate'].startswith('FAIL'))
            self.assertEqual(result['mock_simulation'],'PASS')
            self.assertFalse(result['hardware_validation_passed'])

    def test_all_sources_compile(self):
        for p in ROOT.glob('*.py'):
            compile(p.read_text(encoding='utf-8-sig'),str(p),'exec')

if __name__=='__main__': unittest.main()
