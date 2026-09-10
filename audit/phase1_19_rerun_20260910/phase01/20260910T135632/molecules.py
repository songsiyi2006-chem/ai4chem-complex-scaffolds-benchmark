# -*- coding: utf-8 -*-
"""用 RDKit 计算咖啡因和阿司匹林的分子量、分子式，并保存结构图。"""

import os

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, rdMolDescriptors

MOLECULES = {
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
}


def main():
    mols = []
    legends = []

    for name, smiles in MOLECULES.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"[错误] 无法解析 {name} 的 SMILES: {smiles}")
            continue

        mw = Descriptors.MolWt(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        print(f"{name}: 分子量 = {mw:.2f} g/mol, 分子式 = {formula}")

        mols.append(mol)
        legends.append(name)

    if mols:
        img = Draw.MolsToImage(mols, legends=legends)
        os.makedirs("figures", exist_ok=True)
        img.save(os.path.join("figures", "molecules.png"))
        print(f"已保存 {len(mols)} 个分子的结构图到 figures/molecules.png")


if __name__ == "__main__":
    main()
