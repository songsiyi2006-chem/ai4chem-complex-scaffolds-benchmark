# -*- coding: utf-8 -*-
"""科学级水溶性回归：Delaney ESOL 数据集 + RDKit 特征 + 随机森林。

流程：
1. 自动下载 ESOL 数据集（约 1128 个分子的 SMILES 与实测 logS）；
2. RDKit 计算 4 个物化特征：MolWt / MolLogP / NumRotatableBonds / TPSA；
3. 80/20 划分训练/测试集，训练 RandomForestRegressor；
4. 测试集评估 R² 与 RMSE，绘制实测值 vs 预测值散点图 esol_result.png；
5. 定量预测阿司匹林的 logS。
"""

import os
import urllib.request

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error

RDLogger.DisableLog("rdApp.*")

DATA_URL = ("https://raw.githubusercontent.com/deepchem/deepchem/"
            "master/datasets/delaney-processed.csv")
CSV_PATH = "esol.csv"
ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
FEATURE_NAMES = ["MolWt", "MolLogP", "NumRotatableBonds", "TPSA"]


def download_dataset():
    if os.path.exists(CSV_PATH):
        print(f"发现本地数据集 {CSV_PATH}，跳过下载。")
        return
    print(f"正在下载数据集：{DATA_URL}")
    urllib.request.urlretrieve(DATA_URL, CSV_PATH)
    print(f"下载完成 -> {CSV_PATH}")


def features(smiles):
    """由 SMILES 计算 [MolWt, MolLogP, NumRotatableBonds, TPSA]。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return [
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
        rdMolDescriptors.CalcTPSA(mol),
    ]


def main():
    # 1. 下载数据集
    download_dataset()
    df = pd.read_csv(CSV_PATH)

    smiles_col = next(c for c in df.columns if "smiles" in c.lower())
    target_col = next(c for c in df.columns
                      if "measured log solubility" in c.lower())
    print(f"数据集共 {len(df)} 个分子，特征列: {smiles_col}，标签列: {target_col}")

    # 2. RDKit 提取特征（丢弃解析失败的行）
    feats = df[smiles_col].apply(features)
    valid = feats.notna()
    X = np.array(feats[valid].tolist())
    y = df.loc[valid, target_col].to_numpy()
    print(f"成功解析 {len(X)} / {len(df)} 个分子")

    # 3. 划分数据集并训练随机森林
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42)
    model = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
    model.fit(X_train, y_train)
    print(f"训练集 {len(X_train)} 个，测试集 {len(X_test)} 个")

    # 4. 测试集评估
    y_pred = model.predict(X_test)
    r2 = r2_score(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    print(f"\n测试集评估：R² = {r2:.4f}，RMSE = {rmse:.4f} (log units)")

    importances = model.feature_importances_
    print("特征重要性：")
    for name, imp in sorted(zip(FEATURE_NAMES, importances),
                            key=lambda t: -t[1]):
        print(f"  {name:<18}{imp:.3f}")

    # 5. 实测值 vs 预测值 散点图
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y_test, y_pred, s=18, alpha=0.6, edgecolors="none")
    lims = [min(y_test.min(), y_pred.min()) - 0.5,
            max(y_test.max(), y_pred.max()) + 0.5]
    ax.plot(lims, lims, "r--", lw=1, label="y = x")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Measured logS (mol/L)")
    ax.set_ylabel("Predicted logS (mol/L)")
    ax.set_title(f"ESOL test set: R$^2$ = {r2:.3f}, RMSE = {rmse:.3f}")
    ax.legend()
    ax.set_aspect("equal")
    fig.tight_layout()
    os.makedirs("figures", exist_ok=True)
    fig.savefig(os.path.join("figures", "esol_result.png"), dpi=150)
    print("\n散点图已保存 -> figures/esol_result.png")

    # 6. 预测阿司匹林
    aspirin_feats = features(ASPIRIN)
    aspirin_logs = model.predict([aspirin_feats])[0]
    print(f"\n阿司匹林特征: "
          + ", ".join(f"{n} = {v:.2f}"
                      for n, v in zip(FEATURE_NAMES, aspirin_feats)))
    print(f"预测 logS = {aspirin_logs:.2f} (log mol/L)，"
          f"对应溶解度约 {10 ** aspirin_logs:.3g} mol/L")


if __name__ == "__main__":
    main()
