"""Draw graph identities; does not depict certified axial stereochemistry."""
import json
from rdkit import Chem
from rdkit.Chem import Draw
from .build_inputs import OUT


def main():
    rows=json.loads((OUT/'phase25'/'catalysts.json').read_text(encoding='utf8'))[::2]
    mols=[];legends=[]
    names=['Phenyl','4-Methylphenyl','3,5-Bis(CF3)phenyl','2,4,6-Triisopropylphenyl','9-Anthracenyl','1-Naphthyl']
    for row,name in zip(rows,names):
        mol=Chem.MolFromSmiles(row['mapped_connectivity_smiles'])
        for atom in mol.GetAtoms():atom.SetAtomMapNum(0)
        mols.append(mol);legends.append(row['pair_id']+' | '+name+'\nConnectivity only; axial CIP pending')
    img=Draw.MolsToGridImage(mols,molsPerRow=3,subImgSize=(600,460),legends=legends)
    img.save(str(OUT/'phase25'/'catalyst_connectivities.png'))


if __name__=='__main__':main()
