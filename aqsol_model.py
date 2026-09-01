# -*- coding: utf-8 -*-
"""工业级水溶性预测系统：AqSolDB + ECFP4 指纹 + 梯度提升。

1. 自动下载 AqSolDB（Harvard Dataverse, doi:10.7910/DVN/OVHAW8）；
2. RDKit 提取物化描述符 + 2048 位 Morgan 指纹（ECFP4, 半径 2）；
3. HistGradientBoostingRegressor + 5 折交叉验证；
4. 与旧 ESOL 极简模型公平对比（同一 ESOL 测试集，
   训练前剔除与该测试集重复的分子，避免数据泄漏）；
5. 挑战紫杉醇 logS 预测。
"""

import os
import time
import urllib.request

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors, rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import train_test_split, KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_squared_error

RDLogger.DisableLog("rdApp.*")

AQSOL_URL = ("https://dataverse.harvard.edu/api/access/datafile/3407241")
AQSOL_PATH = "aqsol.csv"
ESOL_PATH = "esol.csv"
PACLITAXEL = ("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C"
              "(C(c5ccccc5)NC(=O)c6ccccc6)O)O)OC(=O)c7ccccc7)"
              "(CO4)OC(=O)C)O)C)OC(=O)C")

DESC_FUNCS = [
    ("MolWt", Descriptors.MolWt),
    ("MolLogP", Crippen.MolLogP),
    ("TPSA", rdMolDescriptors.CalcTPSA),
    ("NumRotatableBonds", rdMolDescriptors.CalcNumRotatableBonds),
    ("NumHDonors", Lipinski.NumHDonors),
    ("NumHAcceptors", Lipinski.NumHAcceptors),
    ("FractionCSP3", rdMolDescriptors.CalcFractionCSP3),
    ("RingCount", rdMolDescriptors.CalcNumRings),
    ("NumAromaticRings", rdMolDescriptors.CalcNumAromaticRings),
    ("HeavyAtomCount", Descriptors.HeavyAtomCount),
    ("FormalCharge", Chem.GetFormalCharge),
]

_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def canonical(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else Chem.MolToSmiles(mol)


def featurize(mol):
    """[11 个描述符 | 2048 位 ECFP4] -> 一维 float 向量。"""
    desc = [f(mol) for _, f in DESC_FUNCS]
    fp = np.zeros(2048, dtype=np.float32)
    DataStructs.ConvertToNumpyArray(_morgan_gen.GetFingerprint(mol), fp)
    return np.concatenate([np.array(desc, dtype=np.float32), fp])


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_aqsol():
    if not os.path.exists(AQSOL_PATH):
        log(f"下载 AqSolDB: {AQSOL_URL}")
        req = urllib.request.Request(AQSOL_URL,
                                     headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=300) as r, \
                open(AQSOL_PATH, "wb") as f:
            f.write(r.read())
        log(f"下载完成 -> {AQSOL_PATH}")
    df = pd.read_csv(AQSOL_PATH, sep="\t")
    df["canon"] = df["SMILES"].map(canonical)
    df = df.dropna(subset=["canon"])
    # 同一分子可能出现在多个来源子库中：按规范 SMILES 去重并取均值
    n0 = len(df)
    df = df.groupby("canon", as_index=False)["Solubility"].mean()
    log(f"AqSolDB: 解析成功 {n0} 条, 去重后 {len(df)} 个分子 "
        f"(标签 Solubility = logS/molL, 范围 {df['Solubility'].min():.2f} ~ "
        f"{df['Solubility'].max():.2f})")
    return df


def build_matrix(smiles_list):
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    return np.array([featurize(m) for m in mols], dtype=np.float32)


def main():
    # 1. 数据与特征
    aqsol = load_aqsol()
    log(f"计算 {len(aqsol)} 个分子 x {len(DESC_FUNCS)} 描述符 + 2048 位 ECFP4 ...")
    X_all = build_matrix(aqsol["canon"].tolist())
    y_all = aqsol["Solubility"].to_numpy()
    log(f"特征矩阵: {X_all.shape}")

    # 2. 5 折交叉验证
    model = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.1,
                                          random_state=42)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    log("5 折交叉验证训练中 ...")
    y_oof = cross_val_predict(model, X_all, y_all, cv=kf, n_jobs=1)
    r2 = r2_score(y_all, y_oof)
    rmse = np.sqrt(mean_squared_error(y_all, y_oof))
    print(f"\n=== AqSolDB 5 折交叉验证 ===")
    print(f"R²   = {r2:.4f}")
    print(f"RMSE = {rmse:.4f} (log units)")

    # 3. 与旧 ESOL 极简模型公平对比：同一测试集 + 泄漏防护
    log("加载 ESOL 测试集并与旧模型对比 ...")
    esol = pd.read_csv(ESOL_PATH)
    esol_smiles = next(c for c in esol.columns if "smiles" in c.lower())
    esol_target = next(c for c in esol.columns
                       if "measured log solubility" in c.lower())
    esol["canon"] = esol[esol_smiles].map(canonical)
    esol = esol.dropna(subset=["canon"])
    idx = np.arange(len(esol))
    _, test_idx = train_test_split(idx, test_size=0.2, random_state=42)
    esol_test = esol.iloc[test_idx]
    y_test = esol_test[esol_target].to_numpy()
    X_test = build_matrix(esol_test["canon"].tolist())

    # 泄漏防护：从 AqSolDB 训练集中剔除与 ESOL 测试集相同的分子
    test_set = set(esol_test["canon"])
    mask = ~aqsol["canon"].isin(test_set).to_numpy()
    log(f"泄漏防护: 从 AqSolDB 剔除 {int((~mask).sum())} 个与 "
        f"ESOL 测试集重复的分子")

    # 旧模型：4 特征随机森林（与 esol_model.py 完全同参）
    old_feats_names = [n for n, _ in DESC_FUNCS[:4]]
    def old_feats(mol):
        return [f(mol) for _, f in DESC_FUNCS[:4]]
    X_old_train = np.array([old_feats(Chem.MolFromSmiles(s))
                            for s in aqsol.loc[mask, "canon"]])
    y_old_train = aqsol.loc[mask, "Solubility"].to_numpy()
    X_old_test = np.array([old_feats(Chem.MolFromSmiles(s))
                           for s in esol_test["canon"]])
    old_model = RandomForestRegressor(n_estimators=500, random_state=42,
                                      n_jobs=-1)
    old_model.fit(X_old_train, y_old_train)
    old_pred = old_model.predict(X_old_test)

    # 新模型：全量特征梯度提升（在泄漏防护后的 AqSolDB 上训练）
    new_model = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.1,
                                              random_state=42)
    new_model.fit(X_all[mask], y_all[mask])
    new_pred = new_model.predict(X_test)

    old_r2, old_rmse = r2_score(y_test, old_pred), np.sqrt(mean_squared_error(y_test, old_pred))
    new_r2, new_rmse = r2_score(y_test, new_pred), np.sqrt(mean_squared_error(y_test, new_pred))
    print(f"\n=== 公平对比（同一 ESOL 测试集, n={len(y_test)}）===")
    print(f"旧模型 (4描述符 + RF, AqSolDB训练): R² = {old_r2:.4f}, RMSE = {old_rmse:.4f}")
    print(f"新模型 (描述符+ECFP4 + HistGBR):     R² = {new_r2:.4f}, RMSE = {new_rmse:.4f}")
    print(f"RMSE 相对下降: {(1 - new_rmse / old_rmse) * 100:.1f}%")

    # 4. 对比散点图
    fig, ax = plt.subplots(figsize=(6.5, 6))
    lims = [min(y_test.min(), new_pred.min()) - 0.5,
            max(y_test.max(), new_pred.max()) + 0.5]
    ax.plot(lims, lims, "k--", lw=1, label="y = x")
    ax.scatter(y_test, old_pred, s=16, alpha=0.5, label=f"Old RF (RMSE={old_rmse:.2f})")
    ax.scatter(y_test, new_pred, s=16, alpha=0.5, label=f"New HistGBR+ECFP4 (RMSE={new_rmse:.2f})")
    ax.set_xlim(lims); ax.set_ylim(lims); ax.set_aspect("equal")
    ax.set_xlabel("Measured logS (mol/L)")
    ax.set_ylabel("Predicted logS (mol/L)")
    ax.set_title("Same ESOL test set: old vs new model")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    os.makedirs("figures", exist_ok=True)
    fig.savefig(os.path.join("figures", "aqsol_result.png"), dpi=150)
    log("对比图已保存 -> figures/aqsol_result.png")

    # 5. 紫杉醇挑战
    pt_mol = Chem.MolFromSmiles(PACLITAXEL)
    pt_X = featurize(pt_mol).reshape(1, -1)
    pt_desc = [f(pt_mol) for _, f in DESC_FUNCS]
    mw = pt_desc[0]
    pt_logS = float(new_model.predict(pt_X)[0])
    ug_ml = 10 ** pt_logS * mw * 1000
    print(f"\n=== 紫杉醇 (Paclitaxel, MW={mw:.1f}) ===")
    print("描述符:", ", ".join(f"{n}={v:.2f}" for (n, _), v in zip(DESC_FUNCS, pt_desc)))
    print(f"预测 logS = {pt_logS:.2f} (log mol/L)")
    print(f"预测溶解度 = {10 ** pt_logS:.3g} mol/L = {ug_ml:.3g} μg/mL")
    print("文献实测: 约 0.3~1 μg/mL (logS 约 -6.1 ~ -5.6)")
    print("旧极简模型 (ESOL 4特征RF): logS = -4.82, 12.8 μg/mL")


if __name__ == "__main__":
    main()
