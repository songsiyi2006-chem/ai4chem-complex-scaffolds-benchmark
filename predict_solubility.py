# -*- coding: utf-8 -*-
"""极简水溶性预测：用 RDKit 计算 MolLogP / TPSA 特征，
训练决策树分类器，预测阿司匹林的水溶性倾向。"""

from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors
from sklearn.tree import DecisionTreeClassifier, export_text

# 名称: (SMILES, 水溶性标签) —— 1 = 易溶, 0 = 难溶
DATA = {
    "乙醇":     ("CCO", 1),
    "葡萄糖":   ("C(C1C(C(C(C(O1)O)O)O)O)O", 1),
    "丙酮":     ("CC(=O)C", 1),
    "苯":       ("c1ccccc1", 0),
    "己烷":     ("CCCCCC", 0),
    "四氯化碳": ("ClC(Cl)(Cl)Cl", 0),
}

ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
FEATURE_NAMES = ["MolLogP", "TPSA"]


def features(smiles):
    """由 SMILES 计算 [MolLogP, TPSA] 特征向量。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法解析 SMILES: {smiles}")
    return [Crippen.MolLogP(mol), rdMolDescriptors.CalcTPSA(mol)]


def main():
    # 1. 构建特征矩阵和标签
    X, y = [], []
    print("训练数据特征：")
    print(f"{'分子':　<6}\t{'SMILES':<28}\tMolLogP\tTPSA\t标签")
    for name, (smiles, label) in DATA.items():
        feats = features(smiles)
        X.append(feats)
        y.append(label)
        print(f"{name:　<6}\t{smiles:<28}\t{feats[0]:6.2f}\t{feats[1]:6.2f}\t{label}")

    # 2. 训练决策树（限制深度，保证可解释）
    clf = DecisionTreeClassifier(max_depth=2, random_state=0)
    clf.fit(X, y)

    # 3. 打印决策树的判断逻辑
    print("\n决策树判断逻辑：")
    print(export_text(clf, feature_names=FEATURE_NAMES))

    # 4. 预测阿司匹林
    aspirin_feats = features(ASPIRIN)
    pred = clf.predict([aspirin_feats])[0]
    proba = clf.predict_proba([aspirin_feats])[0]

    print(f"阿司匹林特征: MolLogP = {aspirin_feats[0]:.2f}, TPSA = {aspirin_feats[1]:.2f}")
    print(f"预测类别概率: 难溶(0) = {proba[0]:.2f}, 易溶(1) = {proba[1]:.2f}")
    verdict = "易溶" if pred == 1 else "难溶"
    print(f"预测结果: 阿司匹林在水里更偏向【{verdict}】(标签 = {pred})")


if __name__ == "__main__":
    main()
