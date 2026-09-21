# Phase31 执行记录

本次是本地确定性角色管线，不是已连接 LLM 的蜂群。
参考体系：BRD4 BD2–VHL / macroPROTAC-1 (6SIS)。
预算上限：100,000,000；本管线实际提供方 Token：0。
该计数不包含编写本项目的 Codex 对话消耗。

| 角色 | 状态 | 证据 |
|---|---|---|
| Alpha | SUCCEEDED | EXECUTED_MMFF_BENCHMARK |
| Beta | BLOCKED | Requested original omegaB97X-D not wired; no reviewed QM partition |
| Gamma | SUCCEEDED | dependency/method blockers are not SCF oscillation; no retries |
| Delta | BLOCKED | No validated forcefield/topology/partition or PLUMED/REST2 sampling adapter |
| Epsilon | BLOCKED | No verified route supplied |
| Zeta | BLOCKED | No matched binary/conditional binding free energies |
| Omega | SUCCEEDED | No candidate passes scientific evidence gates |

没有执行 DFT/CASSCF、QM/MM、FEP、REST2、CI-NEB 或 IRC；没有合格候选。
MMFF94s 图只比较同一分子的收敛构象相对能量，既非绝对张力能，也非结合自由能。
无真实采样，不生成三维自由能曲面；无 SCF 执行，自愈日志为空。
SQLite 仅适用于同机共享磁盘；不能部署到 NFS 充当多节点队列。
