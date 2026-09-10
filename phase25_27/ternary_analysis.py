"""Finite-candidate structural/energy diagnostics. No rates, barriers, or ee."""
import json
import numpy as np
from .build_inputs import OUT
from .pilot_validation import xyz,HARTREE_KCAL


def rmsd(a,b):
    a=a-a.mean(axis=0);b=b-b.mean(axis=0)
    u,_,v=np.linalg.svd(a.T@b);d=np.eye(3);d[-1,-1]=np.linalg.det(u@v)
    return float(np.sqrt(np.mean(np.sum((a@(u@d@v)-b)**2,axis=1))))


def main():
    validation=json.loads((OUT/'phase25/xtb_pilot/validation.json').read_text())['records']
    rows=[r for r in validation if r['name'].startswith('C01_I01') and r['accepted_local_minimum']]
    if not rows:return
    candidates=[];geoms=[]
    ref=min(r['raw_engine_G_hartree'] for r in rows)
    for r in rows:
        symbols,pos=xyz(OUT/'phase25/xtb_pilot'/r['campaign']/r['name']/'xtbopt.xyz')
        geoms.append(pos[np.array(symbols)!='H'])
        candidates.append(dict(name=r['name'],initial_imine_E_Z=r['imine_E_Z_input'],final_imine_E_Z=r['imine_geometry_E_Z'],
          proton_location=r['acid_proton_location_by_distance'],relative_local_G50_kcal_mol=(r['raw_engine_G_hartree']-ref)*HARTREE_KCAL,
          minimum_frequency_cm1=r['lowest_frequency_cm1'],rotor_cutoff_G_shift_kcal_mol=r.get('entropy_only_G_shift_cutoff50_to100_kcal_mol'),
          physical_degeneracy=None,ensemble_weight=None,product_R_S=None))
    matrix=[[rmsd(a,b) for b in geoms] for a in geoms]
    controls=json.loads((OUT/'phase25/xtb_pilot/mirror_singlepoints_01/results.json').read_text())
    energies={r['name']:r['total_energy_hartree'] for r in controls}
    original_symbols,original=xyz(OUT/'phase25/xtb_pilot/mirror_singlepoints_01/xtbopt/input.xyz')
    mirror_symbols,mirror=xyz(OUT/'phase25/xtb_pilot/mirror_singlepoints_01/C01_ternary_mirror/input.xyz')
    reflection=original.copy();reflection[:,0]*=-1
    if original_symbols!=mirror_symbols or not np.allclose(reflection,mirror,atol=2e-12,rtol=0):raise ValueError('Control is not an exact reflection')
    result=dict(candidates=candidates,heavy_atom_fixed_mapping_RMSD_A=matrix,
      rmsd_note='Proper-rotation alignment only, no reflection or atom permutation; this diagnostic does not establish physical basin degeneracy.',
      mirror_control=dict(delta_energy_kcal_mol=(energies['C01_ternary_mirror']-energies['xtbopt'])*HARTREE_KCAL,
        interpretation='Parity/energy implementation check for an exactly reflected full complex. Does not validate sign inversion of product ee.'),
      signed_ee=None,activation_free_energy=None,conformer_search_converged=False,
      warning='Local-minimum ordering must not be inserted into an R/S Eyring selectivity formula.')
    (OUT/'phase25/ternary_diagnostics.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Diagnosed',len(rows),'complete-ternary local minima')


if __name__=='__main__':main()
