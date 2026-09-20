# 产业需求边界 / Industrial requirements

主用户限定为**氯碱电极供应商**。本次没有企业委托、正式客户需求、合作装置数据或报价；下表中的工厂字段保留 `not_obtained`。

| 字段 | 当前可交付物/假设 | 实际企业输入状态 |
|---|---|---|
| 部署位置 | 电极供应商研发阶段的证据审计与设计软件 | not_obtained |
| 交付物 | 来源可追溯的稳定性再分析、失活识别诊断、同产出核算代码 | 本地已交付；不是获批准控制器 |
| 允许停机 | 场景中 48 h 连续生产；不允许用瞬时更换掩盖失效 | not_obtained |
| 最低负荷/爬坡 | 合成演示 0.75–1.25 倍负荷，0.25 倍/h | not_obtained |
| 杂质容忍度 | 需要供应商对膜、电解液和分子层联立定义 | not_obtained |
| 涂层更换与补充 | 按消耗寿命比例分摊涂层成本；未假定自修复 | not_obtained |
| 产品连续性 | 每小时需求等于额定产量；库存相对共同初始缓冲量的偏差≤2额定小时；初始/最终实际库存2小时，容量4小时 | not_obtained |
| 质量 | FE=0.985、纯度=0.998 为场景假设；阈值=0.995 | not_obtained |
| 温度/膜工况 | 简化热状态上限 35°C 仅是软件场景；不同于文献 90°C | not_obtained |
| 功能单位 | 1 t 合格 Cl₂ 对应的总装置输入，不给 NaOH/H₂分配抵扣 | not_obtained |
| 核算边界 | 槽电积分+辅助用电+分子层消耗+载体/膜+维护+处置+不合格预留 | 未含真实工厂固定资本、融资、盐/水、分离物流 |
| 验收负责人 | 应由实际电极供应商及装置方指定 | not_obtained |

## 能耗、联产品与成本解释

公开论文电耗通常以 kWh/t NaOH 表示；本地场景以 kWh/t Cl₂ 表示，不在同一排名中直接比较。场景给出总模块能耗、产量、总成本与逐项费用，未对 NaOH/H₂ 分摊或给副产品收入。其环境项只给电力相关排放敏感性，不是完整 LCA。

Faraday 量纲检查：生成 1 mol Cl₂ 需 2F C。FE=1 时每 1 V 槽电压相当于约 756 kWh/t Cl₂；只有全槽电压变化、同 FE 且同产量时，才能按这一关系算用电变化。30 mV 的单阳极差只能在其他条件确实相同的假设下推导约 22.7 kWh/t Cl₂，不能直接等同于论文的 65 kWh/t NaOH 系统差。

本地电价为合成 0.11±周期项 USD/kWh，无真实地区/年份报价。电力排放系数 0.15–0.65 kgCO₂e/kWh 和材料价格全部为敏感性区间。Nature Synthesis 2026 SI 的 0.07 USD/kWh、100 t Cl₂/day 等同样是其经济场景，不能称为本项目采购报价。

预设“有产业价值”阈值为同产出成本降低≥5%。本次合成敏感性结果未达到此门槛；且缺乏实测参数，因此不能确认商业收益。简单规则优于当前 MPC 的负结果保留，不重新调参把复杂策略包装成优胜者。

English: The target customer is an electrode supplier; all factory-specific requirements remain not_obtained. The output is an R&D audit/analysis package, not a deployed process controller. Scenario cost accounting includes listed non-energy categories but is not a plant-wide TEA or validated LCA.
