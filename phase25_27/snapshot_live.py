"""Freeze only hash-bound completed atomic gradients from an ongoing Hessian."""
from datetime import datetime,timezone
import hashlib
import json
import numpy as np
from .build_inputs import ROOT,OUT


def main():
    source=OUT/'phase25/stationary_pilot/pbe_E_freq_02/atomic_jobs'
    dest=OUT/'phase25/hessian_checkpoint_snapshot';dest.mkdir(exist_ok=True)
    records=[]
    for binding in sorted(source.glob('*.binding.json')):
        key=binding.name.split('.')[0];b=json.loads(binding.read_text(encoding='utf8'))
        inp=source/(key+'.input.json');res=source/(key+'.result.json')
        if hashlib.sha256(inp.read_bytes()).hexdigest()!=key or b['input_sha256']!=key:raise ValueError('Atomic input binding mismatch')
        if hashlib.sha256(res.read_bytes()).hexdigest()!=b['result_sha256']:raise ValueError('Atomic result binding mismatch')
        raw=json.loads(res.read_text(encoding='utf8'));plan=json.loads(inp.read_text(encoding='utf8'))
        spec=raw['input_data']['specification'];expected=plan['specification']
        if not raw['success'] or spec['driver']!='gradient' or spec['model']!=expected['model']:raise ValueError('Atomic model mismatch')
        if spec['model']!={'method':'pbe-d3bj','basis':'def2-svp'}:raise ValueError('Unexpected method')
        if spec['keywords']!=expected['keywords']:raise ValueError('Atomic option mismatch')
        for option,value in {'DFT_RADIAL_POINTS':75,'DFT_SPHERICAL_POINTS':302,'E_CONVERGENCE':1e-9,'D_CONVERGENCE':1e-8,'SCF_TYPE':'DF'}.items():
            if expected['keywords'].get(option)!=value:raise ValueError('Unexpected PES option: '+option)
        if plan['molecule']['molecular_charge']!=0 or plan['molecule']['molecular_multiplicity']!=1:raise ValueError('Unexpected charge/spin')
        np.testing.assert_allclose(raw['input_data']['molecule']['geometry'],plan['molecule']['geometry'],atol=0,rtol=0)
        gradient=np.asarray(raw['return_result'])
        if gradient.size!=84 or not np.all(np.isfinite(gradient)):raise ValueError('Invalid native gradient')
        gradient=gradient.reshape(28,3)  # QCSchema JSON flattens numeric arrays.
        for file in [inp,res,binding]:(dest/file.name).write_bytes(file.read_bytes())
        records.append(dict(input_sha256=key,result_sha256=b['result_sha256'],success=True,gradient_shape=list(gradient.shape)))
    meta=dict(snapshot_utc=datetime.now(timezone.utc).isoformat(),completed_native_gradient_jobs=len(records),expected_native_jobs=163,
      full_hessian_complete=False,frequencies=None,scientific_minimum_acceptance=False,records=records,
      units='QCSchema: geometry bohr; gradient Hartree/bohr',method='PBE-D3BJ/def2-SVP',grid=[75,302],
      source_geometry_sha256=hashlib.sha256((OUT/'phase25/stationary_pilot/pbe_E_opt_01/optimized.xyz').read_bytes()).hexdigest(),
      note='Frozen completed gradients from a still-running native finite-difference Hessian. Partial gradients are not a frequency result.')
    (dest/'manifest.json').write_bytes((json.dumps(meta,indent=2)+'\n').encode())
    status=json.loads((ROOT.parent/'local_campaign/status.json').read_text(encoding='utf8'))
    config=json.loads((ROOT.parent/'local_campaign/config.json').read_text(encoding='utf8'))
    snapshot=dict(snapshot_utc=meta['snapshot_utc'],worker_state=status['state'],current_jobs=status['jobs'],threads=config['threads'],memory_MB=config['memory_mb'],
      queue=config['jobs'],auto_publish=False,requires_local_computer_running=True,
      note='Process status at the snapshot time only; future success is not guaranteed. The finite queue stops after these jobs.')
    (OUT/'phase25/local_execution_snapshot.json').write_bytes((json.dumps(snapshot,indent=2)+'\n').encode())
    print('Frozen',len(records),'completed native gradients; full Hessian not yet accepted')


if __name__=='__main__':main()
