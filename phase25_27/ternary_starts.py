"""Construct intact C01/imine/HE starting complexes; face labels are NOT R/S.

Rigid placements sample imine-N/acid-OH binding and both imine face approaches.
Ranking uses only overlap penalties, never interpreted as physical energies.
"""
import json
import numpy as np
from rdkit import Chem
from scipy.spatial.transform import Rotation
from .build_inputs import ROOT,OUT


def unit(x):return x/np.linalg.norm(x)


def rotation(a,b):
    a,b=unit(a),unit(b);v=np.cross(a,b);c=np.dot(a,b)
    if np.linalg.norm(v)<1e-10:
        if c>0:return np.eye(3)
        q=unit(np.cross(a,[1,0,0] if abs(a[0])<.8 else [0,1,0]))
        return Rotation.from_rotvec(np.pi*q).as_matrix()
    return Rotation.from_rotvec(np.arccos(np.clip(c,-1,1))*unit(v)).as_matrix()


def load(name):
    mols=[m for m in Chem.SDMolSupplier(str(OUT/'phase25/starting_structures'/f'{name}.sdf'),removeHs=False) if m]
    mol=min(mols,key=lambda m:float(m.GetProp('MMFF_energy_kcal_mol')))
    maps={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()}
    return mol,mol.GetConformer().GetPositions(),maps


def write_xyz(path,symbols,xyz,title):
    path.write_bytes((f'{len(symbols)}\n{title}\n'+''.join(f'{s} {x:.12f} {y:.12f} {z:.12f}\n' for s,(x,y,z) in zip(symbols,xyz))).encode())


def main():
    dest=OUT/'phase25/ternary_starts';dest.mkdir(parents=True,exist_ok=False)
    cat,cxyz,cm=load('C01');don,dxyz,dm=load('HE_dimethyl')
    oh=cm[3];ah=next(a.GetIdx() for a in cat.GetAtomWithIdx(oh).GetNeighbors() if a.GetAtomicNum()==1)
    u=unit(cxyz[ah]-cxyz[oh]); dc=dm[416]
    dh=next(a.GetIdx() for a in don.GetAtomWithIdx(dc).GetNeighbors() if a.GetAtomicNum()==1)
    pt=Chem.GetPeriodicTable();records=[]
    for ez in ['E','Z']:
        im,ixyz,imaps=load('I01_'+ez);ni,ci=imaps[203],imaps[202]
        lp=-sum(unit(ixyz[a.GetIdx()]-ixyz[ni]) for a in im.GetAtomWithIdx(ni).GetNeighbors())
        initial=(ixyz-ixyz[ni])@rotation(lp,-u).T
        symbols=[a.GetSymbol() for m in [cat,im,don] for a in m.GetAtoms()]
        maps=[a.GetAtomMapNum() for m in [cat,im,don] for a in m.GetAtoms()]
        radii=np.array([pt.GetRcovalent(a.GetAtomicNum()) for m in [cat,im,don] for a in m.GetAtoms()])
        groups=np.repeat([0,1,2],[cat.GetNumAtoms(),im.GetNumAtoms(),don.GetNumAtoms()])
        pairmask=np.triu(groups[:,None]!=groups[None,:],1)
        for face in [-1,1]:
            candidates=[]
            for ia,angle in enumerate(np.arange(0,360,30)):
                ix=initial@Rotation.from_rotvec(u*np.radians(angle)).as_matrix().T+cxyz[ah]+1.8*u
                normal=face*unit(np.cross(ix[ni]-ix[ci],ix[imaps[210]]-ix[ci]))
                dr=(dxyz-dxyz[dc])@rotation(dxyz[dh]-dxyz[dc],-normal).T
                for phi in np.arange(0,360,30):
                    dx=dr@Rotation.from_rotvec(normal*np.radians(phi)).as_matrix().T+ix[ci]+3.7*normal
                    xyz=np.vstack([cxyz,ix,dx]);dist=np.linalg.norm(xyz[:,None]-xyz[None,:],axis=2)
                    ratio=dist/(radii[:,None]+radii[None,:])
                    minimum=float(ratio[pairmask].min())
                    score=float(np.maximum(1.25-ratio[pairmask],0).dot(np.maximum(1.25-ratio[pairmask],0)))
                    candidates.append((score,int(angle),int(phi),minimum,xyz))
            candidates.sort(key=lambda v:v[0]);selected=[]
            for c in candidates:
                if c[3]<.72:continue
                if any(min(abs(c[1]-p[1]),360-abs(c[1]-p[1]))<60 for p in selected):continue
                selected.append(c)
                if len(selected)==2:break
            if len(selected)<2:raise RuntimeError('Insufficient nonclashing starts')
            for i,(score,angle,phi,minimum,xyz) in enumerate(selected):
                name=f'C01_I01_{ez}_face{face:+d}_start{i}'
                write_xyz(dest/(name+'.xyz'),symbols,xyz,'Rigid ternary start; face is not product CIP; no energy')
                bonds=[];offset=0
                for mol in [cat,im,don]:
                    bonds += [[b.GetBeginAtomIdx()+offset,b.GetEndAtomIdx()+offset,float(b.GetBondTypeAsDouble())] for b in mol.GetBonds()]
                    offset+=mol.GetNumAtoms()
                row=dict(name=name,atom_count=len(symbols),atom_maps=maps,components=groups.tolist(),covalent_bonds=bonds,
                  imine_E_Z=ez,approach_face=face,product_CIP='UNASSIGNED',imine_rotation_deg=angle,donor_rotation_deg=phi,
                  overlap_penalty_NOT_energy=score,min_intercomponent_distance_over_covalent_radius_sum=minimum,
                  acid_H_index=ah,acid_O_index=oh,imine_C_index=cat.GetNumAtoms()+ci,imine_N_index=cat.GetNumAtoms()+ni,
                  donor_H_index=cat.GetNumAtoms()+im.GetNumAtoms()+dh,donor_C_index=cat.GetNumAtoms()+im.GetNumAtoms()+dc,
                  catalyst_axis='same input geometry as C01.xyz; absolute CIP pending',charge=0,multiplicity=1,
                  sampled_placements=144,physical_degeneracy=None,statistical_weight=None,
                  method='rigid placement; no quantum energy; complete catalyst and donor retained')
                records.append(row)
    (dest/'manifest.json').write_bytes((json.dumps(records,indent=2)+'\n').encode())
    print('Built',len(records),'intact ternary starts of',records[0]['atom_count'],'atoms')


if __name__=='__main__':main()
