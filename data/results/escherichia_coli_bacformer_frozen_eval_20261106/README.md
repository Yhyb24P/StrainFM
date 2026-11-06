# E. coli Bacformer 冻结部署验证与闭环 A（官方输入）

状态：**deployment verification**（2026-11-06）。本目录数字只证明“本地部署 + BacBench 官方评测代码可运行、可复算”，**不构成菌株级表征能力或多菌株泛化的科研证据**。测试/验证划分未开封，未运行任何检索协议。

## 内容

```text
config.json                         方法/模型/策略
verification.json                   输入哈希、环境、embedding 检查、资源、命令、指标
operon/official_input/bacformer/     官方 run_embed_prot_seqs.py 输出（parquet，不入 Git）
operon/official_input/esm/           ESM++ 对照输出（parquet，不入 Git）
operon/metrics_bacformer_official_input.csv   官方 run_evaluation_long_read_rna_seq.py 原始输出
operon/metrics_esm_official_input.csv         同上（ESM 对照）
operon/metrics_compact.csv                    标量指标（正式报告引用此文件）
```

## 输入

- 作者原始 E. coli（DH5α，`GCF_002899475.1`，contig `NZ_CP026085.1`）蛋白顺序 + `operon_prot_indices`：
  `data/processed/escherichia_coli_operon_reference_20261106/ecoli_operon_reference.parquet`
  （SHA-256 `b54cdb7c1dbd003ab51a01175d42af810156bec128a5f24eabbb9d055134c487`；4,440 蛋白 / 369 operon / 1 contig）
- 数据集 revision `d77df90f19544475a4da28c0cb80cda2e4d00652`；BacBench commit `64b7d45823b2ff2ef88458e3aa235caa461a12aa`。

## 结果（BacBench 官方协议：(cos+1)/2、异链置零；CPU）

| 方法 | AUROC | average precision | 对 | 正例对 |
|---|---:|---:|---:|---:|
| Bacformer-L 冻结（rev `84f85c3`） | 0.7125 | 0.2194 | 4,439 | 574 |
| ESM++ small 均值（rev `1e40c1b8`） | 0.6724 | 0.1777 | 4,439 | 574 |

代码中的 “AUPRC” 实为 `sklearn.average_precision_score`。单菌株，无 std/CI；**不可外推为多菌株泛化**。

## 验收记录

- 设备：CPU（GPU 运行当时仅 1.6GB 空闲、97% 利用率，按部署策略回退；`device` 已记录）。
- 确定性：同设备复跑 embedding 逐位一致（max abs diff = 0.0）；容差标准 1e-6（见 `docs/BACFORMER_LOCAL_DEPLOYMENT.md` §2.1）。
- 资源：Bacformer 运行 37.13s / 峰值 RSS 9.83GB；ESM 运行 73.21s / 峰值 RSS 5.27GB（OMP_NUM_THREADS=16）。
- embedding：4,440 × 960 float32，无 NaN/Inf；`split` 字段由官方脚本覆盖（已知坑），本评估不依赖该字段。
- 工件说明：`operon/official_input/bacformer/` 同时有首跑与 `_rerun` 文件，两者逐位一致（diff 0.0）；`embeddings` 列是官方单行扁平嵌套列表，消费方需 `reshape(4440, 960)`。
- 版本：torch 2.11.0 / transformers 4.48.1 / bacformer 0.2.0；模型从固定本地快照离线加载（`HF_HUB_OFFLINE=1` + `HF_HUB_CACHE=/mnt/e/StrainFM/models/hf`）。

## 本地 DH5α 输入复算（验收 #3，已完成）

同版本 GBFF 重建输入（4,444 蛋白）经显式映射（4,429/4,440 精确、0 重复）后，在 4,408 个对齐邻接对上复算：Bacformer-L 0.71262/0.22060，官方同对集合 0.71186/0.22041；ESM++ 0.67177/0.17868 vs 0.67176/0.17868。注释漂移无可观测影响。详见 `operon/local_dh5a/README.md`。

## 未完成/未包含

- 官方输入之外的任何多菌株评价；AMR 标签交集不足（D0-18），不在 BacBench 标签文件上建模。
- 检索诊断协议（REQ-002 契约）冻结待用，未生成产物。
