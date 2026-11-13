# 扩大队列三模型对照（2026-11-13 更正版，环丙沙星）：Bacformer / ESM++ / StrainAMR

本目录是 D0-31 预注册第三方对照的**更正重跑**：修复 `_frozen_scores` 未按药物索引的 bug（旧目录 `../escherichia_coli_ast_amr_external_eval_20261112/` 保留为无效证据）。冻结分数来自 `../escherichia_coli_ast_amr_frozen_eval_20261112/` 的 `primary` 分析，现按 **drug → model → seed → genome_id** 索引；测试行筛选不重训任何头。外部分数为 StrainAMR 官方发布模型 `Ecoli_ciprofloxacin_3fold` Fold1 的 `predict --mode model-only` 输出。**仅环丙沙星**：官方清单无庆大霉素模型，庆大霉素保持两模型。

## bug 与更正

- 旧实现把两药预测合并到 `model → seed → genome_id`，同一基因组在两药下评分/标签不同，后写入的 gentamicin 行覆盖 ciprofloxacin 行；输入中跨药同键冲突 3,306 条。
- 修正后 cipro 冻结分数与 `amr-eval` 分药结果一致（Bacformer primary test AUROC 0.9018，ESM 0.6363）；旧外部表为 0.9272/0.7126，属跨药混入而非集成口径差异。
- 回归测试 `tests/test_amr_external_eval.py::test_same_genome_in_two_drugs_does_not_mix_scores`：同一基因组写入两药、gent 分数反向，断言 cipro 指标仍 >0.9。
- 口径：冻结模型点估计 = **先逐 seed 计算指标再平均**（与 `amr-eval` 一致）；StrainAMR 单一检查点。若对 seed 分数先平均再算指标，Bacformer AUROC 为 0.9164（与 gpt 独立复算一致）——两种聚合都在记录中，报告采用前者。

## 协议不变项

| 项 | 值 |
|---|---|
| 外部模型包 | SHA-256 `6f240090453432bbeda7ed1352308c3393bcb9b848d80be1022c60df68608169` |
| 测试基因组 | 环丙沙星 test 全部 643 例；特征构建 54 分钟、CPU model-only 预测 |
| 训练重叠 | 发布模型 1,012 个训练基因组中 **194/643（30.2%）** 与 test 的 BV-BRC ID 重合 |
| 公平集 | `primary_disjoint` = 无冲突且无重叠 = **433 例（228S/205R）**；`primary_all` = 620 例（含重叠，混杂已知） |
| 统计 | 近克隆组重采样 1000 次；StrainAMR 分数跨 seed 不变 |

## 结果（环丙沙星 test）

| 分析集 | 模型 | N (S/R) | AUROC | AP |
|---|---|---:|---|---|
| `primary_all` | Bacformer-L | 620 (392/228) | 0.902 [0.855, 0.938] | 0.810 [0.742, 0.877] |
| `primary_all` | ESM++ 均值 | 620 (392/228) | 0.636 [0.572, 0.696] | 0.519 [0.435, 0.625] |
| `primary_all` | StrainAMR | 620 (392/228) | 0.936 [0.881, 0.975] | 0.943 [0.898, 0.975] |
| `primary_disjoint` | Bacformer-L | 433 (228/205) | 0.885 [0.839, 0.923] | 0.836 [0.761, 0.905] |
| `primary_disjoint` | ESM++ 均值 | 433 (228/205) | 0.614 [0.544, 0.683] | 0.591 [0.502, 0.689] |
| `primary_disjoint` | StrainAMR | 433 (228/205) | 0.930 [0.875, 0.974] | 0.948 [0.908, 0.980] |

**配对差异（同组重采样）**

| 分析集 | 比较 | N | ΔAUROC [95% CI] | ΔAP [95% CI] |
|---|---|---:|---|---|
| `primary_all` | Bacformer − StrainAMR | 620 | −0.034 [−0.100, 0.030] | **−0.133 [−0.209, −0.057]** |
| `primary_all` | ESM − StrainAMR | 620 | **−0.300 [−0.373, −0.220]** | **−0.424 [−0.516, −0.316]** |
| `primary_all` | Bacformer − ESM | 620 | **+0.266 [0.199, 0.330]** | **+0.291 [0.203, 0.372]** |
| `primary_disjoint` | Bacformer − StrainAMR | 433 | −0.045 [−0.107, 0.014] | **−0.112 [−0.183, −0.045]** |
| `primary_disjoint` | ESM − StrainAMR | 433 | **−0.316 [−0.403, −0.228]** | **−0.357 [−0.446, −0.262]** |
| `primary_disjoint` | Bacformer − ESM | 433 | **+0.271 [0.205, 0.337]** | **+0.245 [0.170, 0.319]** |

（`external_paired.tsv` 为精确值；上表四舍五入。）

## 解释

1. **Bacformer 与 StrainAMR 的 AUROC 无可分辨差异**（公平子集 Δ −0.045，CI 跨 0）；**StrainAMR 的 AP 稳定更高**（Δ −0.114，CI 全负），即高分段富集更好。top-50 的 R 精确率 StrainAMR 1.00 vs Bacformer 0.90（诊断性，不参与预注册比较）。
2. **ESM++ 均值在两个参照下都明显落后**（对 StrainAMR ΔAUROC −0.316；对 Bacformer −0.271）。
3. 训练重叠签名：194 个重叠测试单位上 StrainAMR AUROC = 1.000（187 个主分析重叠例）；公平集 433 例为 0.930。`primary_disjoint` 是预注册的公平比较集。
4. 即便无基因组重叠，StrainAMR 仍来自同一 BV-BRC 队列分布；其 AP 优势不能全部归因于特征工程，也不能排除队列邻近性。model-only 模式（机制规则层关闭）与无规则先验的表示模型可比。
5. 对下一步模型工作：环丙沙星上排序能力已接近 StrainAMR，缺口在**高精度端（AP）**；本评测未测概率校准（无 Brier/可靠性曲线），不能把 AP 差异解释为校准差异。庆大霉素无发布模型对照，Bacformer − ESM 的 +0.171 AUROC 是当前最干净的表示层证据。

## 产物与核验

- `external_metrics.tsv` / `external_paired.tsv` / `external_predictions.tsv` / `summary.json` / `contract.json` / `completion.json` / `verification.json`。
- `verification.json` 记录：产物哈希；从 `external_predictions.tsv` 独立重算的 StrainAMR 指标；重叠子集诊断；并列分数统计；以及“逐 seed 平均指标”与“分数先平均”两种口径的 Bacformer/ESM 值，便于与独立复算对照。
- 运行命令：

```bash
PYTHONPATH=src python3 -m strainfm amr-external-eval \
  --eval-dir data/results/escherichia_coli_ast_amr_frozen_eval_20261112 \
  --units data/results/escherichia_coli_ast_joint_split_20261112/amr_task_units.tsv \
  --external-scores /mnt/e/StrainFM/processed/strainamr_ecoli_cipro_20261112/output_sample_prob_predict.txt \
  --external-name strainamr-cipro-frozen --drugs ciprofloxacin \
  --overlap-ids /mnt/e/StrainFM/tools/planA/strainamr_cipro_test_train_overlap_ids.txt \
  --output data/results/escherichia_coli_ast_amr_external_eval_20261113
```

## 限制

- 非受控重训比较；StrainAMR 预训练数据与算力未知；重叠仅按 BV-BRC genome ID 判定（保守口径）。
- 发布模型来自同一 BV-BRC 队列；无重叠子集仍可能存在队列级分布邻近。
- 外部模型单一检查点（Fold1，自报 val AUC 0.98）；冻结头 3 seeds。
- StrainAMR 许可未定（无 LICENSE），仅内部研究用途，不再分发。
- 其余数据/注释限制同 `../escherichia_coli_ast_amr_frozen_eval_20261112/README.md`。
