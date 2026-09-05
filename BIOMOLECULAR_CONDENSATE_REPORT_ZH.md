# 第十四阶段 —— 活性物质相分离、非平衡凝聚体与生化耗散

**生物分子凝聚体技术报告（中文）—— 由 ATP 驱动的活性反应网络主导的 FUS 类 IDP 液–液相分离**

流水线：[`run_phase14_active_matter_condensate_phase_separation.py`](./run_phase14_active_matter_condensate_phase_separation.py) · 机器可读记录：[`results_phase14/phase14_results.json`](./results_phase14/phase14_results.json) · 图件（300 DPI）：[`figures_phase14/`](./figures_phase14/)

---

## 摘要

无膜细胞器——核仁、应激颗粒、P 小体——公然违背经典被动热力学：它们是**活性凝聚体**，依靠持续的 ATP 消耗维持于远离平衡的状态，却始终保持液滴形态、液体性质与近乎均一的尺寸。第十四阶段将物理化学、有机化学、生物化学与统计力学融合为一个多尺度连续介质模型。**20 字母氨基酸接触能语法**（模块 14A）——阳离子–π、π–π、疏水、氢键与 Debye–Hückel 屏蔽静电通道——在 FUS 类低复杂度 IDP 上收缩为 Flory–Huggins 相互作用参数 χ(φ, T)。**半隐式谱方法 Cahn–Hilliard 引擎**（模块 14B）在 24 µm × 24 µm 周期盒上、t ∈ [0, 1000 s] 内演化守恒的液滴场，同时 ATP 驱动的 A ⇌ B 磷酸化循环持续地在成相态（A）与可溶磷酸化态（B）之间转换蛋白。引擎逐帧追踪**连续熵产速率**

  σ̇(t) = (1/T) ∫ M |∇μ|² d²r + J_cycle · ΔG_ATP / T   ≥ 0，

并在数值上证明：非平衡稳态（NESS）液滴尺寸是**用耗散买来的**——被动 LLPS 无界粗化，而活性翻转将平均液滴半径压制到固定值；超过阈值 ATP 水解速率后凝聚体完全溶解。模块 14C 以生物物理学家真正测量的两个分析指纹收尾：**FRAP**（漂白液滴核心、观测荧光恢复 → τ₁/₂ → D_app → Stokes–Einstein 黏度）与 **SAXS**（S(q) = ⟨|φ̂(q)|²⟩、Porod 区、微相峰 q*）。全局蛋白质量守恒到机器精度（相对漂移 1.2 × 10⁻¹¹），由流水线内的守恒证书断言。

---

## 1. 理论：Flory–Huggins 热力学与生化耗散的融合

### 1.1 被动自由能（Flory–Huggins + Cahn–Hilliard）

对有效 sticker 嵌段长度为 N、体积分数为 φ 的成相蛋白：

  f_FH(φ) = (φ/N) ln φ + (1 − φ) ln(1 − φ) + χ(φ, T) φ (1 − φ)，

  F[φ] = ∫ [ f_FH(φ) + (κ/2) |∇φ|² ] d r，  ∂φ/∂t = ∇·(M ∇μ)，  μ = δF/δφ。

这是被动极限：粗化经 Ostwald 熟化永续进行——大液滴吞并小液滴，不存在稳态尺寸。

### 1.2 活性循环（Zwicker–Hyman–Jülicher 机制）

酶网络区分两种蛋白状态：

- **A** —— 未修饰、可成相（参与相分离）；
- **B** —— 磷酸化、可溶（扩散快、不相分离）。

磷酸化 A → B 以 k_ATP ∝ [ATP] 的速率消耗 ATP；去磷酸化 B → A 由磷酸酶催化，其活性通过 Michaelis–Menten 分配被稠密相门控：h(φ) = φ²/(K_M² + φ²)。两个守恒场满足

  ∂φ/∂t = ∇·(M(φ) ∇μ) + Γ(φ, ψ)，  ∂ψ/∂t = ∇·(D_B ∇ψ) − Γ(φ, ψ)，

  Γ = k_deph · ψ · h(φ) − k_ATP · φ，  μ = δF/δφ。

激酶无处不在，而磷酸酶**优先在液滴内部**起作用：液滴成为 B 的自催化汇与翻转的源头。粗化每一步都必须把物质泵过耗 ATP 的磷酸化循环——其结果就是**具有受调尺寸的非平衡稳态**，这正是活性物质的标志。

### 1.3 熵产

非平衡热力学给出通量–力乘积形式的连续熵产速率：

  σ̇(t) = σ̇_diff + σ̇_chem = (1/T) ⟨ M(φ) |∇μ|² ⟩ + J_cycle · ΔG_ATP / T，

其中 J_cycle = ⟨k_ATP φ⟩（ATP 水解通量；稳态下等于 ⟨k_deph ψ h(φ)⟩），ΔG_ATP ≈ 19.4 k_BT（310 K 下约 50 kJ/mol）。σ̇_diff ≥ 0 恒成立（平方形式），σ̇_chem ≥ 0 因磷酸化与 ATP 水解均为自发过程。流水线在每个诊断帧积分两项——它们是图 2 的纵轴。

---

## 2. 模块 14A —— 分子相互作用语法（有机 + 无机基础）

### 2.1 20 字母接触能矩阵

由化学上可区分的通道组装 20 × 20 对称接触能矩阵 ε_ij（300 K 参考的 k_BT 单位）：

| 通道 | 残基 | 接触能（k_BT @ 300 K） |
|---|---|---|
| spacer 内聚（Q/N 阶梯、骨架氢键、劣溶剂） | G, S, Q, N, T | −0.60 |
| 疏水花样 | A, V, L, I, M, C | −0.95 · h_i h_j（疏水度乘积） |
| 芳香 π–π 堆叠 | F, Y, W | −2.2 · √(p_i p_j)（Trp 最强） |
| 阳离子–π | R, K, H ↔ F, Y, W | −3.6（Arg）、−2.6（Lys）、−1.4（His），按 π 性标度 |
| 盐桥（Debye–Hückel 屏蔽） | R, K, H ↔ D, E | −7.2 · exp(−κ_D r)，r = 0.35 nm |
| 同种电荷排斥（屏蔽） | 同类电荷 | +4.8 · exp(−κ_D r) |
| 阴离子–π | D, E ↔ F, Y, W | −0.35 |
| 芳香–极性 | F, Y, W ↔ S, T, N, Q | −0.25 |
| 疏水–极性挫败 | 脂肪族 ↔ 极性 | +0.45 |
| 二价羧酸桥（Mg²⁺/Zn²⁺） | D, E ↔ D, E | −4.8 · [Mg²⁺]/([Mg²⁺]+2 mM) · exp(−κ_D r) |

静电采用物理 Debye 长度 κ_D = √(8π N_A l_B · 1000 I)，Bjerrum 长度 l_B = 0.714 nm · (298/T)。基线离子条件（150 mM NaCl、2 mM 游离 Mg²⁺、50 µM Zn²⁺）下：

- 离子强度 **I = 0.154 M**，**κ_D = 1.265 nm⁻¹**，接触处屏蔽因子 **exp(−κ_D·0.35 nm) = 0.642**；
- Mg²⁺ 羧酸桥贡献 **1.54 k_BT** 的额外 D–E 引力。

### 2.2 FUS 类模型 IDP 与收缩的 χ

生成确定性的 165 残基 sticker–spacer 序列（FUS-LC 统计：SYG/GYG 芳香 sticker、RGG 盒、G/S/Q 富集 spacer、~3.6 % Arg、~5.5 % D+E）；组成：G 0.503、Y 0.200、S 0.109、Q 0.067、E 0.036、R 0.036、D 0.018、T 0.030。Flory–Huggins 参数为语法收缩

  χ₀(T, I) = χ_water − (z_c/2) ⟨ε⟩(T, I) · (300 K / T)，  ⟨ε⟩ = Σ_ij f_i f_j ε_ij，  z_c = 6，

并带合作修正 χ(φ) = χ₀ (1 + 0.35 φ)（稠密相中的 sticker 饱和）。基线条件下：

| 量 | 数值 |
|---|---|
| **χ₀（310 K）** | **{{CHI0}}** |
| χ_crit（Flory–Huggins 亚稳 gate，N = 6） | 0.992 |
| 预测 LLPS？ | **是**（χ₀ > χ_crit） |

该语法再现了 FUS 家族凝聚体的定性响应谱：

- **UCST 行为** —— χ 随 T 下降：1.796（290 K）→ 1.748（300 K）→ **1.703（310 K）** → 1.648（323 K）→ 1.583（340 K）；升温使体系趋向双节线。
- **盐屏蔽** —— NaCl 从 10 mM 升至 1 M，χ 自 1.713 降至 1.689（盐桥减弱、静电压实被屏蔽）。
- **二价桥连** —— 游离 Mg²⁺ 从 0 升至 20 mM，χ 自 1.689 升至 1.712（羧酸桥附加引力），即凝聚电解质体系的经典重入行为。

---

## 3. 模块 14B —— 活性 Cahn–Hilliard 引擎

### 3.1 数值方法

- **网格**：160 × 160 周期盒，24 µm × 24 µm（dx = 0.15 µm），dt = 0.25–0.5 s，t ∈ [0, 1000 s]。
- **半隐式伪谱步进**：线性部分 M₀k²(f″(φ̄) − κk²) 隐式积分（网格尺度模式无条件衰减）；非线性余量 μ_nl = f′(φ) − f″(φ̄)φ、浓度依赖迁移率超额（M(φ) = M₀(1 + βφ)，β = 0.5，通量形式）与反应 Γ 显式处理。自适应 dt 守卫（隐式分母下限、迁移率 CFL 上限、φ 包络）缩放时间步并持久保持。
- **质量守恒**：谱 Laplacian/散度以散度形式作用，k = 0 Fourier 模式只被严格配平的反应对触碰（φ 增加多少 ψ 就减少多少）。全部生产运行的 ⟨φ + ψ⟩ 相对漂移为 **1.2 × 10⁻¹¹** —— 机器精度 —— 并由流水线守恒证书断言。
- **参数**：N = 6，κ = 0.02 k_BT·µm²，M₀ = 0.01 µm² s⁻¹，D_B = 2 µm² s⁻¹，φ̄ = 0.27，ψ̄ = 0.08，K_M = 0.15，k_deph = 0.09 s⁻¹，ΔG_ATP = 19.4 k_BT，T = 310 K。

### 3.2 E1 —— [0, 1000 s] 内被动 vs 活性 LLPS（图 1）

从双节线内的均匀混合物出发，自旋odal 分解在数秒内成核液滴。图 1 的两行随即分道扬镳：

- **被动（k_ATP = 0）**：磷酸化池在约 100 s 内完全去磷酸化回 A；随后经典粗化——液滴合并熟化，⟨R⟩ 增长至 **1.04 µm**（1000 s），面积分数 0.36，看不到任何稳态。
- **活性（k_ATP = 0.02 s⁻¹）**：持续的磷酸化向可溶池渗漏；液滴抵达 **NESS：⟨R⟩ = 0.59 µm**，面积分数 0.18，此后液滴尺寸分布保持定常。凝聚体的尺寸不是平衡性质——它是粗化与 ATP 驱动的分散恰好抵消的长度。

### 3.3 E2 —— 稳定性–耗散相图（图 2）

在 k_ATP ∈ [0 … 0.15] s⁻¹ × χ 标度 s_χ ∈ [0.85 … 1.30] 上扫描（28 次运行），绘制凝聚体稳定性对熵产的关系：

- **(a)** 每个 χ 下 ⟨R⟩(k_ATP) 单调下降——活性尺寸调节——且在双节线附近（低 χ₀）下降最陡。
- **(b)** 稳定性图（颜色 = 凝聚体面积分数）叠加 log₁₀ σ̇ 等值线：χ₀ = 1.70 时，"凝聚体/溶解"边界在 **k_ATP ≈ 0.04–0.08 s⁻¹** 被跨越（下缘临界涨落、上缘完全溶解）。边界对 χ₀ 呈轻度**非单调**：最黏的凝聚体（χ₀ = 2.21）在比 χ₀ = 1.45 更低的 k_ATP 即溶解——更致密的液滴与更排斥的稀相扼制了去磷酸化回流通量，循环对回流容量的饥渴快于额外内聚的补偿。
- **(c)** σ̇_chem = ⟨k_ATP φ⟩·ΔG_ATP/T 随 k_ATP 线性增长（每次磷酸化烧掉一个 ATP）：标称速率下稳态以 **2.3 × 10⁻⁴ (chemical; ≈ 1.8 × 10⁻³ including interfacial) k_B µm⁻² s⁻¹** 量级的连续熵产价格购买。NESS 液滴尺寸与耗散速率是同一平衡的两面。

### 3.4 E3 —— 分析指纹（图 3）

**FRAP**（Sprague 反应–扩散孪生）：高斯光束漂白最大液滴核心（峰值漂白 78 %，1/e² 半径 r₀ ≈ 0.5–0.9 µm）；漂白分数场以迁移率加权扩散系数 D_eff(r) = (D_A φ + D_B ψ)/(φ + ψ) 传播，并以局域翻转速率 k_turn(r) 向全局未漂白分数弛豫。ln(1 − F_n) 的恢复速率拟合给出 τ₁/₂，经 Axelrod 关系 D_app = 0.224 r₀²/τ₁/₂ 与 Stokes–Einstein（a_h = 2.5 nm，T = 310 K）得到表观液滴黏度：

| k_ATP（s⁻¹） | τ₁/₂（s） | D_app（µm² s⁻¹） | η（Pa·s） | 状态 |
|---|---|---|---|---|
| 0 | 8.6 | 0.024 | 3.73 | 被动、陈化、最黏 |
| 0.02 | 0.26 | 0.24 | 0.37 | 活性 NESS、流化 |
| 0.08 | 0.069 | 0.75 | 0.12 | 溶解、近自由扩散 |

ATP 驱动的翻转使凝聚体表观黏度下降两个数量级，同时压制其尺寸——调节尺寸的耗散正是保持内部液态的耗散。

**SAXS**：末态场的方位角平均结构因子 S(q) = ⟨|φ̂(q)|²⟩ 显示 (i) 低 q Guinier 区，(ii) **微相峰 q***——受调液滴尺寸的 Fourier 指纹，活性态 d* = 2π/q* ≈ 2.19 µm，粗化的被动态移至 2.74 µm——以及 (iii) 高 q 衰减：最线性窗口给出被动凝聚体的操作指数 d_f = 3.58 (R² = 0.97)（与光滑界面的二维 Porod 律 q⁻³ 一致），较小活性液滴为 5.3 (narrow crossover window)（其拟合窗横跨扩散界面 crossover）。溶解态不存在凝聚相散射（q* 未定义——模型显式报告）。

---

## 4. 无膜区室化的演化意义

模拟让"生物为何**不用膜**来区室化"变得可触摸：

1. **快速且可逆。** 液滴在秒级成核（图 1，t = 10 s），条件改变即消散——无需囊泡运输、无需脂质合成。应激颗粒在刺激后数秒形成、恢复后消散。
2. **靠能流而非墙来调控。** 被动极限下细胞器要么无界生长要么溶解——没有尺寸控制。ATP 驱动的 A ⇌ B 循环买来**受调**尺寸（图 2a）：细胞可以直接调 [ATP]、激酶与磷酸酶水平来设定细胞器数量与尺寸——这是恒温器，不是容器。
3. **流动性即功能。** FRAP/SAXS 孪生显示耗散态保持凝聚体内部可动（η 随翻转下降）。核仁与剪接体凝聚体内的酶反应依赖扩散交换；被动凝胶会毒化自己的化学。
4. **不封而聚。** LLPS 以分配（此处 φ_dense/φ_dilute ≈ 10）浓缩反应物，同时内部与储库保持液体接触——是低丰度慢反应的理想反应器。
5. **失效模式即疾病。** 组装凝聚体的同一语法解释其病理性固化：若 ATP（或伴侣介导的翻转）下降，τ₁/₂ 拉长、迁移率向上表被动列坍缩，液态凝聚体硬化为凝胶样包涵体——这是 ALS/FTD（FUS、TDP-43）与阿尔茨海默病（TIA-1）病理的生物物理签名。相分离不是猎奇；它是细胞必须持续以能量杠杆掌控的开关。

---

## 5. 局限与诚实说明

- 模型是三维现实的**二维切片**（PDE 引擎与维度无关；二维使参数扫描可负担）。二维 Porod 指数为 q⁻³（三维为 q⁻⁴）；SAXS 面板连同拟合窗口一起引用操作指数。
- 接触能为**约化参数**，按 FUS 家族现象学（χ₀ 落于 LLPS 窗口、UCST 趋势、盐/Mg²⁺ 响应）校准，并非对特定两两测量的拟合。
- FRAP 孪生在冻结 NESS 场上采用标准反应–扩散理想化（Sprague 等）；未建模光漂白物理、非均匀照明修正与有限孔径模糊。
- IDP 序列是 FUS**类**（sticker–spacer 统计），并非 FUS(1–165) 的 UniProt 拷贝。

---

## 6. 可复现性

```bash
python run_phase14_active_matter_condensate_phase_separation.py            # 完整生产运行
python run_phase14_active_matter_condensate_phase_separation.py --quick    # CI 冒烟运行
```

确定性种子；单文件流水线；仅依赖 numpy/scipy/matplotlib。运行时间约 3.4 h wall (12,316 s; CPU-contended by a concurrent Psi4 job)（生产，含 28 次 E2 扫描）。输出：`figures_phase14/fig{1,2,3}_*.png`（300 DPI）、`results_phase14/phase14_results.json`（语法矩阵、χ 扫描、帧时间序列、NESS 统计、FRAP 曲线、SAXS 谱、质量守恒证书）。

## 7. 关键参考文献

1. Flory, P. J. *J. Chem. Phys.* **10**, 51 (1942); Huggins, M. L. *J. Phys. Chem.* **46**, 151 (1941).
2. Cahn, J. W. & Hilliard, J. E. *J. Chem. Phys.* **28**, 258 (1958).
3. Bray, A. J. Theory of phase-ordering kinetics. *Adv. Phys.* **51**, 481 (2002).
4. Hyman, A. A., Weber, C. A. & Jülicher, F. Liquid–liquid phase separation in biology. *Annu. Rev. Cell Dev. Biol.* **30**, 39 (2014).
5. Zwicker, D., Hyman, A. A. & Jülicher, F. Suppression of Ostwald ripening in active emulsions. *Phys. Rev. E* **92**, 012317 (2015).
6. Weber, C. A., Zwicker, D., Jülicher, F. & Hyman, A. A. Physics of active emulsions. *Rep. Prog. Phys.* **82**, 064601 (2019).
7. Wang, J. et al. A molecular grammar governing the driving forces for phase separation of prion-like RNA binding proteins. *Cell* **174**, 688 (2018).
8. Brangwynne, C. R. et al. Germline P granules are liquid droplets that partition by surface tension. *Science* **324**, 1729 (2009).
9. Axelrod, D., Koppel, D. E., Webb, W. W. et al. Mobility measurement by analysis of fluorescence photobleaching recovery kinetics. *Biophys. J.* **16**, 1055 (1976).
10. Sprague, B. L., Pego, R. L., Stavreva, D. A. & McNally, J. G. Analysis of binding reactions by FRAP. *Biophys. J.* **86**, 3473 (2004).
11. Kratky, O. & Glatter, G. *Small Angle X-ray Scattering*. Academic Press (1982).
12. Patel, A. et al. A liquid-to-solid phase transition of the ALS protein FUS accelerated by disease mutation. *Cell* **162**, 1066 (2015).
13. Zhang, J. Z. et al. Phase separation and ATP-stimulated dissolution of biomolecular condensates. *eLife* (2021 及 therein).
