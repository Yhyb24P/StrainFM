# AMR-02：冻结表示 + 论文式轻头（gentamicin / ciprofloxacin）

状态：**探索性结果（2026-11-08）**。预注册近克隆组隔离划分（D0-27/D0-28），冻结主干仅训练轻头；测试集每个配置只评分一次。样本量与组数有限，所有区间较宽，**不构成强结论**。

## 协议

- 单位：每药 × 每 BioSample × 版本化组装；标签 = 已发表二元 AST（0=S，1=R），主分析用 paired-accession 证据（gent 546、cipro 545，达 500/50）；direct-only（各 496）未达到 500，仅作敏感性口径保留
- 分组：near-clone 组整组隔离；来源组件允许跨分区（D0-28，写入 limitation）；冲突组（组内 S/R 混合）主分析排除（gent 6 组 64 例、cipro 3 组 13 例），`with_conflicts` 为敏感性分析
- 划分（seed 20261108，冻结于看模型结果之前）：gent train/val/test = 329/48/105（285/41/82 组），cipro = 396/43/93（284/40/81 组）；均 ≥20 组/分区
- 模型：Bacformer-L `84f85c37` 与 ESM++ `1e40c1b8` 的基因组均值向量（同蛋白集合、共享前 6,000、prodigal 注释）；轻头 = LayerNorm→Dropout(0.1)→Linear(960,1)，BCE，Adam，LR∈{1e-4,1e-3,1e-2} 与早停均在验证集选择
- 指标：AUROC 与 average precision（`sklearn.average_precision_score`）；3 seeds（20261108–10）的均值±标准差；区间为测试集 near-clone 组级 bootstrap（1000 次）；模型差异用同一重采样索引的配对 bootstrap

## 结果（主分析，3 seeds 均值 ± sd）

| 药物 | 模型 | test N (S/R) | AUROC | AP | 验证 AUROC |
|---|---|---:|---:|---:|---:|
| ciprofloxacin | Bacformer-L | 93 (47/46) | 0.636 ± 0.039 | 0.645 ± 0.037 | 0.693 |
| ciprofloxacin | ESM++ 均值 | 93 (47/46) | 0.621 ± 0.014 | 0.653 ± 0.017 | 0.726 |
| gentamicin | Bacformer-L | 105 (92/13) | 0.571 ± 0.075 | 0.221 ± 0.063 | 0.615 |
| gentamicin | ESM++ 均值 | 105 (92/13) | 0.645 ± 0.003 | 0.196 ± 0.036 | 0.653 |

组级 bootstrap（95% 区间，seed 20261108 的预测）：

| 药物 | 模型 | AUROC 区间 | AP 区间 |
|---|---|---|---|
| ciprofloxacin | Bacformer-L | 0.547–0.799 | 0.526–0.843 |
| ciprofloxacin | ESM++ 均值 | 0.482–0.764 | 0.490–0.830 |
| gentamicin | Bacformer-L | 0.491–0.820 | 0.111–0.487 |
| gentamicin | ESM++ 均值 | 0.440–0.809 | 0.088–0.362 |

**配对模型差异（Bacformer − ESM，同一组重采样索引）**

| 药物 | ΔAUROC（3 seeds 均值±sd；boot CI） | ΔAP（均值±sd；boot CI） |
|---|---|---|
| ciprofloxacin | +0.015 ± 0.028；[−0.057, 0.166] | −0.008 ± 0.027；[−0.080, 0.131] |
| gentamicin | −0.073 ± 0.060；[−0.273, 0.337] | +0.026 ± 0.079；[−0.137, 0.308] |

含冲突敏感性（`with_conflicts`）与主分析方向一致，见 `amr_metrics.tsv` / `amr_model_comparison.tsv`。

## 解释与限制

1. **未观察到 Bacformer 上下文表示优于匹配 ESM 均值的证据**：两药配对差的 95% 区间均跨 0；gentamicin 上 ESM 的 AUROC 名义更高、Bacformer 的 AP 名义更高，但区间很宽。该结果只适用于本协议、本数据版本与本轻头。
2. 测试集小（93/105 例；gent 仅 13 R），组级 CI 宽；不能据此断言“模型无能力”。
3. 标签为已发表二元 AST；paired-accession 是数据库配对而非独立序列核验；折叠前的实验条件与冲突不可观测；direct-only 口径（496）未达 500 门槛。
4. 来源组件跨分区（D0-28）；near-clone 组隔离已满足，但共享项目不等于谱系独立。
5. 注释统一为 pyrodigal（D0-23）；蛋白 >1,024aa 截断 33,895 个；基因组均值池化与 Bacformer 论文的下游用法可能存在差异（未做逐蛋白池化敏感性）。
6. 预训练暴露 unknown；本结果不涉及“未见菌株”或 “未见谱系”泛化声明。

## 文件

`contract.json`（协议与哈希）、`summary.json`（含逐模型 bootstrap）、`amr_metrics.tsv`（逐 seed）、`amr_predictions.tsv`（测试集逐样本分数）、`amr_model_comparison.tsv`（配对差异）。

## 下一步（待用户决定）

- 采 NCBI AST 扩大独立标签；或
- 逐蛋白/注意力池化等表示层敏感性；或
- 停止 AMR 主线，回到其他真实标签任务。
