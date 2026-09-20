"""Compare matched fixed-geometry gas-phase state differences, never raw method energies."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ID_MAP={'Q01':'P01','Q02':'P04','Q03':'P09'}


def matching_gas_pairs(dft_rows,xtb_rows):
    rows=[]
    for dft in dft_rows:
        if dft['environment']!='gas':
            raise ValueError('DFT comparison is restricted to gas phase')
        matches=[r for r in xtb_rows if r['molecule_id']==ID_MAP[dft['molecule_id']] and r['environment']=='gas']
        if len(matches)!=2 or {int(r['gfn']) for r in matches}!={1,2}:
            raise ValueError('Expected complete GFN1/GFN2 gas pair')
        for xtb in matches:
            if dft['geometry_sha256']!=xtb['geometry_sha256']:
                raise ValueError('DFT/xTB geometry mismatch')
            dft_delta=float(dft['delta_E_radical_minus_cation_eV'])
            xtb_delta=float(xtb['delta_model_eV'])
            rows.append({'molecule_id':dft['molecule_id'],'matrix_id':ID_MAP[dft['molecule_id']],
                         'basis':dft['basis'],'gfn':int(xtb['gfn']),'environment':'gas',
                         'geometry_sha256':dft['geometry_sha256'],'PBE0_delta_E_eV':dft_delta,
                         'GFN_delta_E_eV':xtb_delta,'GFN_minus_PBE0_eV':xtb_delta-dft_delta,
                         'DFT_radical_spin_squared':float(dft['radical_spin_squared']),
                         'interpretation':'raw unaligned charge-state reference; not predictive error or potential'})
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dft',type=Path,required=True)
    parser.add_argument('--xtb-csv',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    dft_file=args.dft/'energy_differences.csv'
    dft_summary=json.loads((args.dft/'summary.json').read_text(encoding='utf8'))
    dft_rows=list(csv.DictReader(dft_file.open(encoding='utf8')))
    xtb_rows=list(csv.DictReader(args.xtb_csv.open(encoding='utf8')))
    rows=matching_gas_pairs(dft_rows,xtb_rows)
    if not rows:raise ValueError('No completed matched pairs')
    args.out.mkdir(parents=True,exist_ok=False)
    with (args.out/'comparison.csv').open('w',encoding='utf8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    relative=[]
    for row in rows:
        if row['basis']!='def2-svp':continue
        parent=next(r for r in rows if r['molecule_id']=='Q01' and r['basis']==row['basis'] and r['gfn']==row['gfn'])
        relative.append({'molecule_id':row['molecule_id'],'gfn':row['gfn'],
                         'PBE0_attachment_difference_from_Q01_eV':row['PBE0_delta_E_eV']-parent['PBE0_delta_E_eV'],
                         'GFN_attachment_difference_from_Q01_eV':row['GFN_delta_E_eV']-parent['GFN_delta_E_eV']})
    summary={'definition':'Delta E = E(neutral radical) - E(+1 cation), identical fixed nuclei within each pair',
             'evidence':'executed gas-phase model comparison; PBE0 is not experimental ground truth',
             'reference_warning':'Raw xTB charge-state differences have unaligned electron references. No empirical IP correction applied. Use Q01-centered double differences to cancel a method-wide constant offset.',
             'DFT_summary':dft_summary,'matched_comparisons':len(rows),'relative_to_parent':relative,
             'sources_sha256':{str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in (dft_file,args.xtb_csv)},
             'GFN_minus_PBE0_ranges_eV':{}}
    for gfn in (1,2):
        values=[r['GFN_minus_PBE0_eV'] for r in rows if r['gfn']==gfn and r['basis']=='def2-svp']
        summary['GFN_minus_PBE0_ranges_eV'][str(gfn)]={'min':min(values),'max':max(values)}
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    molecules=list(ID_MAP)
    for method in ('GFN1','GFN2','PBE0/def2-SVP'):
        values=[]
        for mol in molecules:
            if method.startswith('GFN'):
                row=next(r for r in rows if r['molecule_id']==mol and r['basis']=='def2-svp' and r['gfn']==int(method[-1]))
                values.append(row['GFN_delta_E_eV'])
            else:
                row=next(r for r in rows if r['molecule_id']==mol and r['basis']=='def2-svp')
                values.append(row['PBE0_delta_E_eV'])
        axes[0].plot(molecules,values,'o-',label=method)
        axes[1].plot(molecules,[v-values[0] for v in values],'o-',label=method)
    axes[0].set(ylabel='Raw charge-state difference / eV',title='Electron references unaligned; not an error metric')
    axes[1].set(ylabel='Difference from Q01 / eV',title='Centered contrasts cancel a constant reference')
    for ax in axes:
        ax.set_xlabel('Q01 parent | Q02 3-OMe | Q03 3-CN');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Fixed-nuclei gas-phase models: no redox potential or site prediction',fontsize=11)
    fig.tight_layout();fig.savefig(args.out/'matched_model_comparison.png',dpi=170);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
