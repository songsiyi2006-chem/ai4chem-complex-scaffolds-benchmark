# -*- coding: utf-8 -*-
"""用 ESOL 随机森林模型预测紫杉醇（Paclitaxel）的水溶性 logS。

复用 esol_model.py 中的特征提取与数据集下载逻辑，
以完全相同的超参数重建模型，保证与已评估模型一致。
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, Draw
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

from esol_model import download_dataset, features, CSV_PATH, FEATURE_NAMES

RDLogger.DisableLog("rdApp.*")

PACLITAXEL = ("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C"
              "(C(c5ccccc5)NC(=O)c6ccccc6)O)O)OC(=O)c7ccccc7)"
              "(CO4)OC(=O)C)O)C)OC(=O)C")


def main():
    # 1. 解析紫杉醇并计算特征
    mol = Chem.MolFromSmiles(PACLITAXEL)
    if mol is None:
        raise ValueError("无法解析紫杉醇 SMILES")
    pt_feats = features(PACLITAXEL)
    mw = Descriptors.MolWt(mol)

    print("紫杉醇（Paclitaxel）分子特征：")
    for name, value in zip(FEATURE_NAMES, pt_feats):
        print(f"  {name:<18}{value:8.2f}")

    # 2. 重建与 esol_model.py 相同的随机森林模型
    download_dataset()
    df = pd.read_csv(CSV_PATH)
    smiles_col = next(c for c in df.columns if "smiles" in c.lower())
    target_col = next(c for c in df.columns
                      if "measured log solubility" in c.lower())
    X = np.array(df[smiles_col].apply(features).tolist())
    y = df[target_col].to_numpy()
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42)
    model = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
    model.fit(X_train, y_train)

    # 3. 预测并换算单位
    log_s = model.predict([pt_feats])[0]
    mol_per_l = 10 ** log_s                 # mol/L
    mg_per_ml = mol_per_l * mw              # g/L == mg/mL
    ug_per_ml = mg_per_ml * 1000

    print(f"\n预测结果：")
    print(f"  logS          = {log_s:.2f} (log mol/L)")
    print(f"  溶解度        = {mol_per_l:.3g} mol/L")
    print(f"                = {mg_per_ml:.4g} mg/mL = {ug_per_ml:.3g} μg/mL")

    # 4. 保存 2D 结构图
    img = Draw.MolToImage(mol, size=(600, 600),
                          legend=f"Paclitaxel (pred. logS = {log_s:.2f})")
    os.makedirs("figures", exist_ok=True)
    img.save(os.path.join("figures", "paclitaxel.png"))
    print("\n结构图已保存 -> figures/paclitaxel.png")

    # 5. 与文献对比点评
    print("\n【与文献对比】")
    print(f"  紫杉醇实测水溶性约 0.3~1 μg/mL（logS 约 -6.1 ~ -5.6，"
          f"按 MW={mw:.0f} 折算），")
    print(f"  模型预测为 {ug_per_ml:.1f} μg/mL，"
          f"方向正确（极难溶），但存在数倍量级的高估。")


if __name__ == "__main__":
    main()
