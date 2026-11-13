# 结果目录

统一入口：[模型、架构与历次结果 HTML](index.html)；机器可读数据：[catalog.json](catalog.json)。页面数字来自现存聚合表，不从 README 的四舍五入数字反推。

**已更正：** `escherichia_coli_ast_amr_external_eval_20261112/` 的冻结分数读取缺少 drug 键，同基因组的两药分数发生覆盖，该目录的点估计/配对比较不能用于模型排名，原表原址保留。更正版 `escherichia_coli_ast_amr_external_eval_20261113/` 已按 drug/model/seed/genome_id 修复并统一逐种子指标平均（commit `dc48648`，含同基因组两药回归测试）。修复前用正确 ciprofloxacin 行、分数先平均口径复算的 620 例值 Bacformer AUROC/AP 为 0.9164205/0.8257622、ESM 为 0.6265888/0.5200288，与更正版 `verification.json` 并存记录。

| 用途 | 位置 | 阅读规则 |
|---|---|---|
| 最新两药表示比较 | `escherichia_coli_ast_amr_frozen_eval_20261112/` | 主干冻结；逐种子指标平均；相同队列、划分与轻量头 |
| 环丙沙星专项模型参照 | `escherichia_coli_ast_amr_external_eval_20261112/` | 发布模型；先排除已知训练 ID 重叠；非共同训练预算 |
| 环丙沙星更正结果 | `escherichia_coli_ast_amr_external_eval_20261113/` | 修正跨药索引并统一逐种子指标平均；旧外部比较失效，更正版仍非受控重训 |
| 历史小队列与敏感性 | `escherichia_coli_amr_frozen_eval_20261108/`、`…_20261109/` | primary 复用预测；direct_refit 改变样本并重拟合，不视为新增独立测试 |
| 单菌株 operon | `escherichia_coli_bacformer_frozen_eval_20261106/operon/` | 评价单位为基因对；官方与本地对齐输入分开 |
| 部署 / 检索诊断 | 含 `smoke` 的目录 | 小样本诊断，不能替代功能任务 |
| 数据、标签与划分证据 | 其他审计、来源、质量、split 目录 | 是模型结果的支撑，不计作模型性能 |
| 第三方快照 | `bacbench_pinned_20261106/` 等本地忽略目录 | 不再分发，也不是第二套正式实现 |

新运行沿用 `物种_任务或对象_条件_日期/` 命名。指标表、预测、配置/contract、摘要、completion 和简要 README 按运行归属；不把新结果覆盖写入旧运行。大规模输入和权重保留在冷数据目录。

历史目录已有哈希契约，本次不移动。`catalog.json` 列出每个目录的用途、已读指标数量和指标输入 SHA-256；无指标的运行仍保留为记录。报告会收录本地 smoke22 的聚合表，但不会自动把其大型或未跟踪产物纳入 Git。

重建：`PYTHONPATH=src python3 -m strainfm report`。HTML 内嵌聚合数据，无 CDN 或运行时网络依赖。逐种子平均与预测集成分别标记；配对区间不冒充预测集成区间。不同队列的差值不解释为训练数据规模的因果收益，AP 不解释为概率校准。
