# 扩大队列三模型对照（2026-11-12，环丙沙星）：Bacformer / ESM++ / StrainAMR

本目录是 D0-31 预注册的第三方对照评测。冻结分数来自 `../escherichia_coli_ast_amr_frozen_eval_20261112/` 的 `primary` 分析（测试行筛选不重训任何头）；外部分数来自 StrainAMR 官方发布模型 `Ecoli_ciprofloxacin_3fold` Fold1 的 `predict --mode model-only` 输出。**仅环丙沙星**：StrainAMR 官方清单的 24 个物种-药物模型中没有庆大霉素模型（E. coli 仅 amoxicillin/ampicillin/ciprofloxacin/levofloxacin），庆大霉素保持两模型比较。

## 协议（看结果前写入 D0-31）

| 项 | 值 |
|---|---|
| 外部模型包 | `Ecoli_ciprofloxacin_3fold__Fold1.zip`，SHA-256 `6f240090453432bbeda7ed1352308c3393bcb9b848d80be1022c60df68608169`（含 2.7GB 训练库与检查点；下载器校验） |
| 运行方式 | 官方 `build-test`（pyrodigal + RGI/CARD + PC/k-mer/graph 通道）→ `predict --mode model-only --skip-analysis`；CPU 推理（cu102 不支持 RTX 5090 的 sm_120） |
| 测试基因组 | 环丙沙星 test 全部 643 例（decompressed FASTA）；特征构建 54 分钟、预测数分钟 |
| 训练重叠 | 发布模型训练集 1,012 个 BV-BRC 基因组中，**194/643（30.2%）**与我们 test 的 BV-BRC ID 重合（逐代表 genome 取 `source_bvbrc_genome_ids` 交集） |
| 公平集 | `primary_disjoint` = 无标签冲突且不在其训练集中的测试单位 = **433 例（228S/205R）**；`primary_all` = 620 例（含重叠，混杂已知） |
| 统计 | 同一近克隆组重采样，1000 次；冻结头为 3 seeds 平均，外部分数单一检查点 |

重叠 ID 清单：`/mnt/e/StrainFM/tools/planA/strainamr_cipro_test_train_overlap_ids.txt`；外部原始分数：`/mnt/e/StrainFM/processed/strainamr_ecoli_cipro_20261112/output_sample_prob_predict.txt`。

## 结果（环丙沙星 test）

**点估计（组级 bootstrap 95% CI）**

| 分析集 | 模型 | N (S/R) | AUROC [95% CI] | AP [95% CI] |
|---|---|---:|---|---|
| `primary_all` | Bacformer-L | 620 (392/228) | 0.927 [0.890, 0.952] | 0.845 [0.773, 0.915] |
| `primary_all` | ESM++ 均值 | 620 (392/228) | 0.713 [0.634, 0.768] | 0.608 [0.492, 0.732] |
| `primary_all` | StrainAMR（冻结） | 620 (392/228) | 0.936 [0.881, 0.975] | 0.943 [0.898, 0.975] |
| `primary_disjoint` | Bacformer-L | 433 (228/205) | 0.916 [0.875, 0.946] | 0.863 [0.780, 0.937] |
| `primary_disjoint` | ESM++ 均值 | 433 (228/205) | 0.709 [0.631, 0.774] | 0.682 [0.557, 0.802] |
| `primary_disjoint` | StrainAMR（冻结） | 433 (228/205) | 0.930 [0.875, 0.974] | 0.948 [0.908, 0.980] |

**配对差异（同组重采样；StrainAMR 分数跨 seed 不变）**

| 分析集 | 比较 | N | ΔAUROC [95% CI] | ΔAP [95% CI] |
|---|---|---:|---|---|
| `primary_all` | Bacformer − StrainAMR | 620 | −0.010 [−0.065, 0.050] | **−0.101 [−0.178, −0.021]** |
| `primary_all` | ESM − StrainAMR | 620 | **−0.234 [−0.309, −0.161]** | **−0.338 [−0.451, −0.214]** |
| `primary_all` | Bacformer − ESM | 620 | **+0.224 [0.166, 0.280]** | **+0.237 [0.143, 0.330]** |
| `primary_disjoint` | Bacformer − StrainAMR | 433 | −0.016 [−0.073, 0.043] | **−0.087 [−0.163, −0.011]** |
| `primary_disjoint` | ESM − StrainAMR | 433 | **−0.226 [−0.306, −0.150]** | **−0.263 [−0.386, −0.147]** |
| `primary_disjoint` | Bacformer − ESM | 433 | **+0.209 [0.152, 0.269]** | **+0.176 [0.099, 0.264]** |

## 解释

1. **Bacformer 与 StrainAMR 的 AUROC 无可分辨差异**（公平子集 Δ −0.016，CI 跨 0；全测试集同）。但 **StrainAMR 的 AP 稳定更高**（Δ −0.087，CI 全负），即高分段富集更好；后验检查：公平子集 top-50 的 R 精确率 StrainAMR 1.00 vs Bacformer 0.90（诊断性，不参与预注册比较）。
2. **ESM++ 均值在两个参照下都明显落后**（ΔAUROC −0.23 对 StrainAMR、−0.21 对 Bacformer）。
3. **训练重叠的签名清晰**：194 个重叠测试单位上 StrainAMR AUROC = 1.000（187 个主分析重叠例），非重叠 433 例为 0.930。因此全测试集数字不能作为“泛化”证据；`primary_disjoint` 是预注册的公平比较集。
4. 即便在无基因组重叠子集，StrainAMR 仍来自同一 BV-BRC 队列分布（同标签来源、同物种、时间段相近），其优势可能部分来自队列级邻近性；本评测不能进一步拆分“机制规则/特征工程”与“数据邻近性”的贡献。StrainAMR 以 model-only 模式运行（机制规则层关闭），与无规则先验的表示模型可比。
5. 对下一步模型工作的含义：环丙沙星上，**达到 StrainAMR 的排序能力不是主要缺口，缺口在高精度端（AP）与校准**；庆大霉素没有发布模型可对照，两个表示模型的差距（ΔAUROC +0.171）是当前最干净的表示层证据。

## 产物与核验

- `external_metrics.tsv` / `external_paired.tsv` / `external_predictions.tsv` / `summary.json` / `contract.json` / `completion.json`。
- `verification.json`：产物哈希复核；从 `external_predictions.tsv` 独立重算的 StrainAMR 指标；重叠子集诊断；分数并列统计。
- 运行命令：

```bash
PYTHONPATH=src python3 -m strainfm amr-external-eval \
  --eval-dir data/results/escherichia_coli_ast_amr_frozen_eval_20261112 \
  --units data/results/escherichia_coli_ast_joint_split_20261112/amr_task_units.tsv \
  --external-scores /mnt/e/StrainFM/processed/strainamr_ecoli_cipro_20261112/output_sample_prob_predict.txt \
  --external-name strainamr-cipro-frozen --drugs ciprofloxacin \
  --overlap-ids /mnt/e/StrainFM/tools/planA/strainamr_cipro_test_train_overlap_ids.txt \
  --output data/results/escherichia_coli_ast_amr_external_eval_20261112
```

## 限制

- 非受控重训比较：StrainAMR 的预训练数据与算力预算未知；冻结头在本次 train 划分训练、验证选择。
- 训练重叠只按 BV-BRC genome ID 判定（代表基因组的任一来源 ID 命中即算重叠），是保守口径；未做序列级近重复复核。
- StrainAMR 发布模型训练集来自同一 BV-BRC 队列；无重叠子集仍可能存在队列级分布邻近。
- 外部模型单一检查点（Fold1，其自报 val AUC 0.98）；冻结头为 3 seeds 平均。
- StrainAMR 许可未定（仓库无 LICENSE），仅内部研究用途，不对外再分发模型或数据。
- 其余数据与注释限制同 `../escherichia_coli_ast_amr_frozen_eval_20261112/README.md`。
