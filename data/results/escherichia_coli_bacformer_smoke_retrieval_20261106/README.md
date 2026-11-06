# Phase 2A smoke：22 基因组冻结表示与检索（pipeline verification only）

状态：**smoke / pipeline verification**（2026-11-06）。用于验证 adapter → 冻结向量 → 检索评测的端到端链路与产物契约，**不构成科研结论**：gallery 只有 12 个训练基因组，指标基于 3 个有正例的验证查询。

## 输入与链路

| 环节 | 产物 |
|---|---|
| 蛋白输入 | `data/processed/escherichia_coli_bacformer_inputs_smoke22_20261106/`（22 基因组：19 GBFF + 3 个 0-CDS GBFF 回退 pyrodigal；89,724 蛋白；SHA-256 见其 summary） |
| 冻结向量 | `data/processed/escherichia_coli_bacformer_embeddings_smoke22_20261106/`（22 × 960，两种模型同蛋白集合；无 NaN；0 截断；契约含设备变更与代码哈希历史） |
| 检索 | 本目录 `bacformer/` 与 `esm_mean/`（协议 A，split=validation，k={1,5,10}） |

基因组：10 个验证集（高完整度、1 contig）+ 12 个同邻域训练 gallery；3 个查询（`GCA_048839225.1`、`GCA_048841005.1`、`GCA_048842605.1`）各有 4 个已嵌入正例。测试集未使用。

## 指标（3 个有正例查询；4 正例/查询）

| 模型 | Hit@1 | Hit@5 | Recall@1 | Recall@5 | Recall@10 | MRR | mAP@10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Bacformer-L 冻结 | 0.667 | 1.000 | 0.167 | 0.917 | 1.000 | 0.833 | 0.856 |
| ESM++ small 均值 | 0.333 | 1.000 | 0.083 | 0.333 | 0.833 | 0.500 | 0.397 |

开放集：本目录 `open_set_n = 7` 是**小 gallery 造成的伪开放集**（这 7 个验证查询在完整训练集中有正例，见 `retrieval_query_support.tsv`），不是数据交付口径的 val 160 / test 153。正式运行需嵌入完整训练 gallery 后重算。

## 运行条件

- Bacformer-L `84f85c37`（CPU, bf16）；ESM++ `1e40c1b8`（前 16 个基因组 cuda:0，后 6 个按用户指示 CPU 续跑；契约 `resume_runs` 记录）。
- 检索实现：`src/strainfm/eval_retrieval.py`（qw 实现，协议 A/B、余弦距离、ties 按 genome_id、开放集单列）。
- 单查询结果：各自目录下 `retrieval_results.parquet`；指标表 `retrieval_metrics.tsv`；契约 `contract.json`。

## 未包含

- 全量开发集（train gallery 26,462 / validation 2,113）运行；正式开放集计数；评分稳定性与置信区间；第二/三层评价（更细遗传距离、实验标签）。测试集继续封存。
