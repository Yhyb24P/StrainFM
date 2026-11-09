# AMR-02 补充运行（2026-11-09）：direct_refit 敏感性与跨种子配对 CI

本目录是对 2026-11-08 首轮结果（`../escherichia_coli_amr_frozen_eval_20261108/`，保留不改写）的**统计补全**：①新增 D0-27 计划的 direct-only 证据受限敏感性（`direct_refit`）；②配对模型差异改为“同一近克隆组重采样索引、3 seeds 配对差先平均再取区间”；③修复首轮运行中配对 bootstrap 的排序 bug；④单类别重采样时 AUROC/AP 同时跳过并报告有效重复次数。首轮点估计不变，二者差异仅在 CI 统计量与新增分析集。

## 分析集定义

| 分析集 | 训练/测试数据 | 问题 |
|---|---|---|
| `primary` | paired-accession 标签；排除近克隆组内 S/R 冲突组 | 主分析 |
| `with_conflicts` | 含冲突组 | 冲突组敏感性 |
| `direct_refit` | 仅 label_source=file_accession（各 496），沿用冻结组分区，**重新拟合**轻头 | 证据口径敏感性；不声称达到 500 门槛 |

## 结果（3 seeds 均值 ± sd；区间 = 组级 bootstrap，seed 平均后取百分位）

| 药物 | 分析集 | 模型 | test N (S/R) | AUROC | AP |
|---|---|---|---:|---:|---:|
| ciprofloxacin | primary | Bacformer-L | 93 (47/46) | 0.636 ± 0.039 | 0.645 ± 0.037 |
| ciprofloxacin | primary | ESM++ 均值 | 93 (47/46) | 0.621 ± 0.014 | 0.653 ± 0.017 |
| ciprofloxacin | direct_refit | Bacformer-L | 83 (44/39) | 0.590 ± 0.067 | 0.607 ± 0.102 |
| ciprofloxacin | direct_refit | ESM++ 均值 | 83 (44/39) | 0.626 ± 0.027 | 0.614 ± 0.039 |
| gentamicin | primary | Bacformer-L | 105 (92/13) | 0.571 ± 0.075 | 0.221 ± 0.063 |
| gentamicin | primary | ESM++ 均值 | 105 (92/13) | 0.645 ± 0.003 | 0.196 ± 0.036 |
| gentamicin | direct_refit | Bacformer-L | 95 (87/8) | **0.495 ± 0.164** | 0.102 ± 0.029 |
| gentamicin | direct_refit | ESM++ 均值 | 95 (87/8) | **0.755 ± 0.100** | 0.190 ± 0.067 |

**配对差异（Bacformer − ESM；CI 与 gpt 独立复算逐项一致）**

| 药物 | 分析集 | N | ΔAUROC [95% CI] | ΔAP [95% CI] |
|---|---|---:|---|---|
| ciprofloxacin | primary | 93 | +0.015 [−0.087, 0.121] | −0.008 [−0.095, 0.088] |
| ciprofloxacin | with_conflicts | 96 | +0.020 [−0.076, 0.127] | −0.008 [−0.103, 0.094] |
| ciprofloxacin | direct_refit | 83 | −0.036 [−0.163, 0.088] | −0.007 [−0.141, 0.114] |
| gentamicin | primary | 105 | −0.073 [−0.317, 0.209] | +0.026 [−0.165, 0.233] |
| gentamicin | with_conflicts | 115 | −0.021 [−0.245, 0.177] | +0.019 [−0.128, 0.159] |
| gentamicin | direct_refit | 95 | **−0.260 [−0.479, 0.009]** | −0.088 [−0.258, 0.066] |

有效重复：除 gentamicin `direct_refit` 为 **999/1000**（1 次重采样单类别，已同时跳过 AUROC 与 AP）外，其余均为 1000/1000。

## 解释

1. **主分析维持首轮结论**：两药均未见 Bacformer 优于匹配 ESM；配对 CI 跨 0。
2. **direct_refit 出现方向性信号（探索性）**：gentamicin 上 Bacformer 均值 AUROC 0.495（近似随机、seed 间极不稳定 ±0.164），ESM 0.755；ΔAUROC −0.260，CI 上界 0.009 接近但不排除 0。测试集仅 8 个 R、LR 选择在验证集完成，**不能据此断言模型缺陷**；这是下一步（更多标签或池化敏感性）最值得检查的方向。
3. ciprofloxacin 在两套证据口径下均无模型差异；direct_refit 的 N=83 是冲突组与 paired-only 单位共同减少的结果。

## 限制

- 标签为已发表二元 AST；paired 是数据库配对非独立序列核验；direct_refit 各 496 未达 500 门槛，仅作敏感性。
- 来源组件跨分区（D0-28）；跨药共享 BioSample 的分区不同（542 共享，243 跨区、160 跨药 train/test），当前逐药独立头无泄漏，但未来共享训练编码器必须统一划分。
- 统一 pyrodigal 注释；>1,024aa 蛋白截断；基因组均值池化。
- 优化器：本实现 `torch.optim.Adam`（无显式 weight_decay）；pinned BacBench AMR 脚本使用 AdamW（同样无显式 weight_decay）。不声称精确复现论文。
- 预训练暴露 unknown；近克隆内不同 BioSample 的表型差异不视为标签错误。

## 文件

`contract.json`、`summary.json`（含逐模型 bootstrap 与有效重复）、`amr_metrics.tsv`（逐 seed）、`amr_predictions.tsv`（逐样本）、`amr_model_comparison.tsv`（配对差异）、`verification.json`（与 gpt 独立复算的逐项对照）。

**verification 性质**：`verification.json` 中的对照是在**同一批预测文件**上对统计量实现的复算，用于核验代码与口径一致性；它不是独立确认性试验，也不构成新的科研证据。
