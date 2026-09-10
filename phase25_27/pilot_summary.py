"""Summarize only audited low-cost results; no selectivity or rate is inferred."""
import json
from datetime import datetime,timezone
from .build_inputs import OUT
from .pilot_validation import main as audit,HARTREE_KCAL,xyz,dihedral


def main():
    audit();val=json.loads((OUT/'phase25/xtb_pilot/validation.json').read_text())
    rows=val['records'];comparisons=[]
    for solvent in ['ALPB(toluene)','ALPB(ch2cl2)']:
        for substrate in ['I01','I02','I03','I04']:
            candidates=[r for r in rows if r['name'].startswith(substrate+'_') and r['solvation']==solvent and r['accepted_local_minimum']]
            best={}
            for ez in ['E','Z']:
                pool=[r for r in candidates if r.get('geometric_E_Z')==ez and r.get('raw_engine_G_hartree') is not None]
                if pool:best[ez]=min(pool,key=lambda r:r['raw_engine_G_hartree'])
            if len(best)!=2:continue
            diff=lambda field:best['Z'][field]-best['E'][field]
            comparisons.append(dict(substrate=substrate,solvation=solvent,
              E_source=best['E']['campaign']+'/'+best['E']['name'],Z_source=best['Z']['campaign']+'/'+best['Z']['name'],
              Z_minus_E_energy_including_ALPB_kcal_mol=diff('energy_including_ALPB_hartree')*HARTREE_KCAL,
              Z_minus_E_local_minimum_G_cutoff50_kcal_mol=diff('raw_engine_G_hartree')*HARTREE_KCAL,
              Z_minus_E_local_minimum_G_cutoff100_kcal_mol=(diff('raw_engine_G_hartree')*HARTREE_KCAL+diff('entropy_only_G_shift_cutoff50_to100_kcal_mol')) if all('entropy_only_G_shift_cutoff50_to100_kcal_mol' in r for r in best.values()) else None,
              interpretation='GFN2-xTB single-local-minimum E/Z difference; not activation energy, conformer ensemble, reaction ee, or validated experimental prediction',
              standard_state_note='Equal-composition 1-to-1 comparison: common standard-state shifts cancel. Absolute solution G not certified.',
              statistical_uncertainty=None))
    dft=[]
    for file in sorted((OUT/'phase25/stationary_pilot').glob('*/result.json')):
        r=json.loads(file.read_text());entry={k:r.get(k) for k in ['tag','method','stage','status','electronic_energy_hartree','wall_seconds','max_abs_gradient_hartree_per_bohr','rms_gradient_hartree_per_bohr']}
        if r['stage']=='opt' and r['status']=='OPTIMIZED_FREQUENCY_PENDING':
            elements,pos=xyz(file.parent/'optimized.xyz')
            if len(elements)==28:
                entry['optimized_aryl_C_N_aryl_dihedral_deg']=dihedral(pos[[9,1,2,3]])
                entry['final_E_Z_geometry']='E' if abs(entry['optimized_aryl_C_N_aryl_dihedral_deg'])>90 else 'Z'
        dft.append(entry)
    costs=[]
    for ledger in sorted((OUT/'phase25/xtb_pilot').glob('*/results.json')):
        entries=json.loads(ledger.read_text());costs.append(dict(campaign=ledger.parent.name,
          completed_jobs=sum('wall_seconds' in r for r in entries),pending_jobs=sum('wall_seconds' not in r for r in entries),
          finished_allocated_core_seconds=sum(r.get('allocated_core_seconds',0) for r in entries)))
    result=dict(snapshot_utc=datetime.now(timezone.utc).isoformat(),xTB_audited_jobs=val['total_jobs'],xTB_accepted_preliminary_minima=val['accepted_preliminary_minima'],
      E_Z_comparisons=comparisons,DFT_jobs=dft,xTB_costs=costs,
      full_ternary_starts=8,atoms_per_ternary=117,
      scientific_acceptance='NOT_PASSED',signed_ee=None,absolute_product='selectivity_undetermined',
      complete_catalytic_R_S_paths=0,explicit_solvent_free_energy_repeats=0,
      caveat='No reduction of model discrepancy to a statistical confidence interval. Search hits are not physical degeneracies.')
    (OUT/'phase25/pilot_progress.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Paired E/Z comparisons:',len(comparisons),'DFT jobs:',len(dft))


if __name__=='__main__':main()
