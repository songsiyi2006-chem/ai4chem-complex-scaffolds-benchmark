# 第23阶段：材料内计算与分子忆阻器——可复现的有效模型基准

## 定位与实际结果

本阶段是经典数值模拟，不是已制备的双三联吡啶/多金属氧酸盐器件，不是量子输运，
也不是实现 AGI 的证据。配位化学可以提供可调的氧化还原和离子弛豫性质，但本程序
的参数是指定值，未由某一真实配合物、第一性原理或实验拟合得到。
材料内计算有望减少部分任务的数据搬运；“硅路线不可逾越”及全面能效优势不能据此推出。

| Quantity / 指标 | Computed value / 计算值 |
|---|---:|
| One-step test NMSE (mean x,y,z) | 0.25545867 |
| One-step test NMSE x / y / z | [0.10287173447306878, 0.26713650512751624, 0.39636777835630116] |
| Ten-step test NMSE (mean x,y,z) | 0.98501813 |
| One-step persistence baseline, x only | 0.010697335 |
| Device dissipation J / input sample | 2.9845378e-09 |
| Device + wire + source-resistor J / sample | 1.3640139e-08 |
| Full system J/op | Not calculated / 未计算 |
| High-frequency log(area)-log(f) slope | -0.999997 |
| Hysteresis area refinement relative error | 1.71153e-05 |
| Positive / negative STDP sign tests | True / True |
| KCL max residual (A) | 9.98575e-14 |
| Halved network-step max state difference | 0.000129068 |
| Halved plasticity-step max normalized-window difference | 1.60583e-07 |

## 23A：为什么必须修改原公式

原式 I=G(w)V sinh(alpha V) 是电压的偶函数，负电压时 I V<0，不能用作这里的被动耗散器件。
采用 I=G(w)sinh(alpha V)/alpha，确保小电压电导为 G，且任意电压 I V>=0。
状态方程为：

    q = sign(V) max(|V|-V_threshold,0)
    dw/dt = eta [1-(2w-1)^(2p)] sinh(q/V0) - (w-w_eq)/tau

保留 Joglekar 窗口，p=1、w_eq=0.5、tau=10 ms。正反指数反应速率相减产生 sinh，
所以这是受 Butler-Volmer 启发的模型，而非定量推导的界面动力学；电流和状态变化
没有通过法拉第电荷守恒关联，V0 是有效参数。没有计算热涨落、反应自由能或化学热。
原模型在 w=0 会锁死；向内部平衡态恢复后，两个边界的导数都指向区间内部。
RK4 根据速率上界细分步长，越界时报错，不用裁剪状态掩盖数值不稳定。

频率覆盖 1 Hz—100 kHz，通过周期轨道射击法排除任意初态瞬变。
报告两瓣面积绝对值之和，避免有向积分相互抵消。高频时状态摆幅为 O(1/f)，
趋于冻结状态的单值 I-V 曲线；实际输出有限频率斜率和加密网格误差。
过原点的滞回与忆阻动力学相容，但不是鉴定真实忆阻器的唯一实验证据。

## 23B：脉冲塑性与长期保持的边界

每个脉冲由 +0.35 V、1 ms 的头部，以及初始 -0.15 V、8 ms 衰减常数的尾部构成，
尾部在 51 ms 截断。器件电压定义为 V_post-V_pre，交换电极会交换学习窗口方向。
扫描 -50 至 +50 ms 共101个延迟，统一积分同一个器件方程，不直接植入正负指数学习规则。

Delta W=int[G_pair(t)-G0]dt 的量纲为 S·s，不是电导；图上无量纲量为 Delta W/(G0 tau)，
不能误标为 Delta G/G0。统一 230 ms 窗口包含所有脉冲及至少69 ms的最终恢复。
3—40 ms尾部进行对数空间指数拟合，并保存拟合质量；这不是“实验 Hebbian 保真度”。
另外导出最终状态和最大状态偏移。PPF 定义为第二脉冲后电导/第一脉冲后电导，
不是扣除基线后的增量之比。零偏压下状态恢复到0.5，因此只能称挥发性 STP/STDP-like，
不能声称长期 LTP/LTD；需要额外慢变量、非易失态和保持实验才能研究长期学习。

## 23C：25节点物理网络与预测协议

5×5节点通过40条电阻边连接，每个节点有接地非线性忆阻器、有限输入电阻和2 pF寄生电容。
固定随机输入增益、偏置、电导和速率，以及3—30 ms弛豫时间代表指定的器件差异。
输入是仅用训练段标准化的 Lorenz x，经过固定电压掩码分配到节点，不加入延迟嵌入、
多项式或手工特征。拓扑是带类交叉阵列寄生项的电阻网格，不是完整行列选通交叉阵列。

每个子步用 Newton 法和正定线性求解器求解：

    (C/h)(v-v_previous)+G_source(v-u)+G_wire L v+I(v,w)=0

然后在保持的节点电压下推进分子状态。整体是后向Euler电容加算子分裂，
全局精度不是四阶。导出KCL残差和子步减半差异。两个初始状态在相同直流输入下
收敛用于展示衰退记忆，但不是任意输入下的全局回声状态稳定性证明。

Lorenz参数为10、28、8/3，每0.02无量纲时间映射为0.5 ms物理采样，是指定时钟，
不是实测计算速度。丢弃最初1500个生成样本后保留3600个输入；洗去200个样本，
训练截止2200，验证截止2800，最后800个样本完全留作测试。
边界预留10个样本，防止未来标签跨分区泄漏。标准化、读出中心化和缩放只用训练集。
只在验证集比较9个岭回归系数；测试集不参与器件参数、读出或超参数选择。

通过SVD训练一个25×3线性矩阵及偏置，预测未来1步和10步的x、y、z。
NMSE分别按测试段每个坐标的方差归一化，同时报告三坐标平均值。
另提供只使用当前输入的线性基线和x坐标的持续性基线。
这是持续输入真实x的驱动式预测，不是自主闭环长期生成混沌。
单步x预测容易因信号平滑而获得低误差，不等于三维或较长步长的预测性能。
本次为确定性单种子结果，未计算跨器件噪声或多个独立混沌样本的置信区间。

## 能耗：必须说明每次操作是什么

一次操作定义为整个25节点网络接收一个输入样本，不随意拆成内部事件以缩小数值。
分别积分器件VI耗散、线阻损耗和输入电阻损耗，并保存电源提供能量、电容储能变化、
后向Euler人工数值损耗，验证同一离散电路的能量平衡。
右端点求积与KCL中保持状态一致；离散守恒残差不是连续时间能耗精确性的证明。
尚未计入ADC、DAC、掩码分配、读出、训练、复位、噪声裕度、制造及化学反应热，
所以报告的是电路子总量，系统总能耗为空；不得标成已实现亚飞焦耳系统。
精度和能耗是否达到原定阈值均由程序计算并保存，未通过也保留。

## 23D：阻抗分析模型

分别求 -0.25、0、+0.25 V下的稳定平衡态，对状态方程F和电流线性化：

    Y_j = I_V + I_w F_V/(i omega-F_w)
    Z_W = R_W tanh(sqrt(i omega tau_D))/sqrt(i omega tau_D)
    Z = R_s + [i omega C_geom + (1/Y_j + Z_W)^(-1)]^(-1)

明确扫描的是f=0.01—1e6 Hz，角频率为2πf。保存名义10 mV交流扰动的电流相量。
Warburg采用有限长度、透射边界，低频电阻有限，并非阻塞边界下发散的赝电容。
R_s=500欧姆、C_geom=2 nF、R_W=30千欧、tau_D=0.5 s都是指定的薄膜/电解质参数，
与网络的单节点寄生电容不同。单一红氧状态无法产生连续扩散谱，故Warburg是额外
等效电路假设，不是从前述状态ODE直接推导出来的扩散。
冻结状态R_ct=1/I_V与完整直流微分电阻不同。输出稳定性检查和10 mV静态割线误差，
并未完成非线性交流谐波仿真、参数反演或实验谱峰解卷积。图谱也不能验证配位环境。
偏压下的状态响应项可能产生感性符号的小环，不能把所有弧线都解释为纯电容半圆。

## 实验路线与复现

下一步应选择明确组成的氧化还原薄膜，测量脉冲、保持、温度依赖和电化学状态，
独立测定线阻/电容，再开展含噪声的化学传感边缘计算测试。
分子电子传感接口是值得验证的方向；亚飞焦耳AGI协处理器仍是展望，不是本阶段成果。

安装 requirements_phase23.txt，运行脚本的 --self-test，再正常运行。
支持 --output-dir 指定输出目录；脚本本身不联网、不执行Git、不自动改写README。
交付包含全部数据、岭回归验证记录、读出矩阵、数据划分、摘要、阈值判断、SHA-256校验
及四张300 DPI图。提交和发布由独立的人工可审查流程完成。

## References / 参考文献

1. Joglekar & Wolf, *The elusive memristor: properties of basic electrical circuits*:
   https://arxiv.org/abs/0807.3994 — window-model background, not parameter calibration.
2. Du et al., *Reservoir computing using dynamic memristors for temporal information processing*:
   https://www.nature.com/articles/s41467-017-02337-y — physical short-term memory and readout-only learning.
3. Zhong et al., *Dynamic memristor-based reservoir computing for high-efficiency temporal signal processing*:
   https://www.nature.com/articles/s41467-020-20692-1 — experimental reservoir work; its measured accuracy is not ours.
4. Goswami et al., *Decision trees within a molecular memristor*:
   https://www.nature.com/articles/s41586-021-03748-0 — molecular redox-state computing; a different device/task.
5. Gamry Instruments, *Basics of electrochemical impedance spectroscopy*:
   https://www.gamry.com/application-notes/EIS/basics-of-electrochemical-impedance-spectroscopy/
   — small-signal equivalent circuits and diffusion impedance.
