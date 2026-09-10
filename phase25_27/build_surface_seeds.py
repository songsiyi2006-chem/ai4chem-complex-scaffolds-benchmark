"""Build dry unrelaxed Cu slabs only. These are NOT electrochemical interfaces."""
from ase.build import fcc100,fcc111,fcc211
from ase.io import write
from .build_inputs import OUT,dump


def main():
    out=OUT/'phase26'/'dry_slab_seeds';out.mkdir(parents=True,exist_ok=True)
    slabs={'100':fcc100('Cu',size=(4,4,4),a=3.61,vacuum=15),
           '111':fcc111('Cu',size=(4,4,4),a=3.61,vacuum=15,orthogonal=True),
           '211':fcc211('Cu',size=(6,4,4),a=3.61,vacuum=15)}
    records=[]
    for facet,slab in slabs.items():
        slab.info.pop('adsorbate_info',None)
        slab.info['evidence_kind']='dry_unrelaxed_seed'; slab.info['potential']='undefined'
        slab.arrays['atom_map']=__import__('numpy').arange(1,len(slab)+1)
        write(out/f'Cu{facet}.extxyz',slab)
        records.append(dict(facet=facet,atoms=len(slab),cell_A=slab.cell.tolist(),lattice_parameter_A=3.61,
          lattice_parameter_origin='provisional conventional Cu value; converge with chosen functional before production',
          water_molecules=0,ions=0,CO_molecules=0,electrode_potential_V=None,relaxed=False))
    dump(out/'manifest.json',records)


if __name__=='__main__':main()
