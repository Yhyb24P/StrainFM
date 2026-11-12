# 扩大队列联合评测（2026-11-12）：Bacformer vs ESM++ 均值

本目录是在**冻结联合划分**（D0-30，`../escherichia_coli_ast_joint_split_20261112/`，units SHA-256 `387f9a4e…`）上，对扩大后的 E. coli AMR 队列做两药冻结主干 + 轻头的评测。StrainAMR 第三方对照见 `../escherichia_coli_ast_amr_external_eval_20261112/`（D0-31）；本目录只包含 ESM++ 均值与 Bacformer-L 两个表示模型。

## 数据与冻结条件

| 项 | 值 |
|---|---|
| 身份代表基因组 | 3,456（train/validation/test = 2,418/344/694 基因组） |
| 注释 | 统一 pyrodigal 3.7.1；16,915,861 蛋白；0 失败；manifest SHA-256 `0f96d033a35d66008b3fa71308f8ae63e7b5c46627a122cb56d3aa74abab49c2` |
| 蛋白截断 | 前 6,000 蛋白（contig 顺序）；1 个基因组触发截断；>1,024aa 蛋白 218,820 条 |
| Bacformer 快照 | `84f85c37a05e55559142b842c1569e265183b554`（bf16，cuda:0） |
| ESM++ 快照 | `1e40c1b8ef46c33f93a9b817eb0bd81279ab4088`（bf16，cuda:0） |
| `genome_embeddings.parquet` | SHA-256 `7b496a6fb5ddf7083bc97e01b60840bceb7929315b3287a2d9cf22708034b002` |
| `esm_mean_embeddings.parquet` | SHA-256 `d92cecf21c96ea7f6138a1b204a5b7e06b9f398d5f76c272edccc78b7baf7293` |
| 主分析规模 | cipro train/val/test = 1,931/286/620；gent = 1,912/272/609（冲突组排除后） |
| 划分隔离 | 近克隆组整组隔离；两药共用分区；来源组件跨分区（D0-30，来源独立主张不支持） |

嵌入目录：`/mnt/e/StrainFM/processed/escherichia_coli_ast_embeddings_20261112/`（`summary.json`/`completion.json`，3,456 全成功、0 失败）。注释目录：`/mnt/e/StrainFM/processed/escherichia_coli_ast_inputs_20261112/`。

## 两模型结果（primary；3 seeds 均值，区间 = 近克隆组 bootstrap 95%）

| 药物 | 模型 | test N (S/R) | AUROC [95% CI] | AP [95% CI] |
|---|---|---:|---|---|
| ciprofloxacin | Bacformer-L | 620 (392/228) | 0.902 [0.856, 0.936] | 0.810 [0.746, 0.877] |
| ciprofloxacin | ESM++ 均值 | 620 (392/228) | 0.636 [0.580, 0.704] | 0.519 [0.430, 0.624] |
| gentamicin | Bacformer-L | 609 (505/104) | 0.893 [0.853, 0.925] | 0.569 [0.476, 0.682] |
| gentamicin | ESM++ 均值 | 609 (505/104) | 0.722 [0.665, 0.777] | 0.382 [0.301, 0.492] |

`with_conflicts`（含混合标签近克隆组）：cipro Bacformer 0.869 / ESM 0.680；gent Bacformer 0.882 / ESM 0.674，配对结论方向一致。

**配对差异（Bacformer − ESM++ 均值；同一重采样组、seed 平均后取区间）**

| 药物 | 分析集 | N | ΔAUROC [95% CI] | ΔAP [95% CI] |
|---|---|---:|---|---|
| ciprofloxacin | primary | 620 | **+0.266 [0.198, 0.328]** | **+0.291 [0.198, 0.376]** |
| ciprofloxacin | with_conflicts | 643 | +0.189 [0.130, 0.245] | +0.197 [0.131, 0.274] |
| gentamicin | primary | 609 | **+0.171 [0.114, 0.224]** | **+0.188 [0.098, 0.277]** |
| gentamicin | with_conflicts | 655 | +0.208 [0.135, 0.280] | +0.206 [0.107, 0.298] |

所有配对区间在两种分析集和两药上均排除 0；无需单类别跳过，有效重复均 1000/1000。cipro Bacformer 的 seed 点估计为 0.887/0.931/0.887（优化随机性，划分固定）。

## 与 AMR-02 的关系

AMR-02（`../escherichia_coli_amr_frozen_eval_20261109/`，546/545 单位）在相同方法下**未观察到稳定差异**（配对 CI 跨 0）。本次队列扩大到约 6 倍（2,837/2,793 主分析单位），测试集 620/609 例、R 例数 228/104（AMR-02 为 46/13）。两者不矛盾：AMR-02 的区间宽度覆盖了本次观测到的差值；本次结果提高的是估计精度，而非改变方法或划分规则。AMR-02 的 `direct_refit` 探索性线索（gent，N=95、8 个 R）仍未被本次主分析直接检验——本次 gent 主分析测试 R 为 104 例，但标签口径是队列主口径而非 direct-only。

## 解释与边界

1. **Bacformer 的种内上下文表征在两药上稳定优于均值池化 ESM**：AUROC 与 AP 同时提高，且在 `primary` 与 `with_conflicts` 两口径一致。该差异是表示层面的，不能归因于头架构（两模型同架构、同调参预算、同划分）。
2. 该结果支持“单菌株基因组内蛋白上下文（操纵子/邻近基因）携带 AMR 相关信号”，但**不区分**是基因共现、水平转移元件邻近性还是其他上下文结构。
3. 结论限定于 pyrodigal 注释、前 6,000 蛋白均值池化、BV-BRC 二元 AST 标签和本划分；未做 MIC、未做物种专属继续预训练。
4. 与 AMR-02 的 3 个 BioSample 重叠（D0-30）是研究流程暴露，不是本评测的训练泄漏：本目录的两个头均在本次 train 划分上重新训练，旧头未参与。

## 产物与核验

- `amr_metrics.tsv`：逐药、逐分析集、逐模型、逐 seed 的指标与选择超参（`selected_lr`/`selected_epoch`）。
- `amr_predictions.tsv`：逐 seed 测试预测（比较与独立复算的输入）。
- `amr_model_comparison.tsv`：配对差异与区间。
- `summary.json` / `contract.json` / `completion.json`：输入哈希、代码哈希、分析集定义与产物哈希。
- `verification.json`：产物哈希复核、从 `amr_predictions.tsv` 独立重算的 AUROC/AP、缺失嵌入计数。

复现：

```bash
PYTHONPATH=src python3 -m strainfm amr-eval \
  --embeddings-dir /mnt/e/StrainFM/processed/escherichia_coli_ast_embeddings_20261112 \
  --units data/results/escherichia_coli_ast_joint_split_20261112/amr_task_units.tsv \
  --output data/results/escherichia_coli_ast_amr_frozen_eval_20261112
```

## 限制

- 标签为 BV-BRC 已发表二元 AST，单一来源；无 MIC/折点复核；混合标签近克隆组只做敏感性。
- 来源组件跨分区（最大单来源占比 45–48%），测试集非来源独立；不能据此主张跨来源泛化。
- 统一 pyrodigal 注释；circular 基因组不旋转；蛋白顺序 = contig 顺序 + 起点。
- Bacformer 预训练语料对这些基因组的暴露未知；本结果度量的是种内结构恢复能力，不是未见语料泛化。
- 头训练为论文式轻头（LayerNorm+Dropout+Linear，Adam），与 pinned BacBench 的 AdamW 有已记录偏差；不声称精确复现。
