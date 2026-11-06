# 本地 DH5α 输入 operon 复算（闭环 A 验收 #3）

状态：**deployment verification**（2026-11-06）。比较对象是作者官方输入与“同版本 GBFF 重建的本地输入”，用于分离注释漂移与流水线差异；不构成多菌株证据。

## 输入

| 工件 | 值 |
|---|---|
| GBFF（NCBI FTP，同版本） | `GCF_002899475.1_ASM289947v1_genomic.gbff.gz`，SHA-256 `0dacb3a72ac2698079b0eb02298b4ec4efd88dbd580e87ac27434d6bdfbd72b0` |
| 本地输入 parquet | `data/processed/escherichia_coli_dh5a_local_input_20261106/ecoli_dh5a_local_input.parquet`，SHA-256 `c8962229bd768b983fd716c1d0ae5d0b501f57d40c7694e69f25ebb91891dfd0`（4,444 蛋白） |
| 显式映射 | `data/processed/escherichia_coli_dh5a_local_mapping_20261106/mapping.tsv`，SHA-256 `214912c761c36b8f125b64f10d47b20e285a09880c8500749f062fec23c35183` |
| 作者参考 | `ecoli_operon_reference.parquet`，SHA-256 `b54cdb7c…c487`（4,440 蛋白 / 369 operon） |

映射：4,429/4,440 按 (protein_id, 坐标) 精确映射、序列一致、0 重复；11 个作者蛋白未映射（author_index 123/1330/2379/3230/3646/3828/3856/3865/3923/3924/4118）；15 个本地蛋白无作者对应；369 operon 全部完整（成员均在映射中）。

## 协议

- 只评“对齐邻接对”：本地相邻 (l, l+1) 且映射到作者连续索引 (a, a+1)；标签取官方作者成对标签；分数用官方 `predict_pairwise_operon_boundaries`（(cos+1)/2、异链置零）。
- 排除分类（不静默裁剪）：`local_extra_protein_pair` = 29（含 15 个本地新增蛋白的邻接对）；`author_unmapped_between` = 6（跨未映射作者蛋白）；对齐对 4,408 / 作者 4,439 对（99.30%）。排除项几乎全部是负例对，基率与全对集合略有差异。
- 对照：把官方输入分数投影到**同一对齐对集合**（P1 口径），本地 vs 官方的差异只来自 embedding 输入本身。

## 结果（CPU）

| 方法 | 本地输入 AUROC / AP | 官方同对集合 AUROC / AP | ΔAUROC / ΔAP |
|---|---:|---:|---:|
| Bacformer-L（rev `84f85c3`） | 0.71262 / 0.22060 | 0.71186 / 0.22041 | +0.0008 / +0.0002 |
| ESM++ small 均值（rev `1e40c1b8`） | 0.67177 / 0.17868 | 0.67176 / 0.17868 | ~0 / ~0 |

结论：在本组装上，同版本 GBFF 重建输入与作者蛋白列表在 99.3% 对齐对上给出几乎相同的 operon 边界恢复表现；注释漂移对结论无可观测影响（差值与单菌株噪声同量级）。Bacformer 相对匹配 ESM 基线的增益在本地输入上保持（+0.041 AUROC）。

## 复现

```bash
PYTHONPATH=src:data/results/bacbench_pinned_20261106/BacBench \
python -m strainfm operon-local-eval \
  --embeddings <local embed parquet> --mapping <mapping.tsv> \
  --reference data/processed/escherichia_coli_operon_reference_20261106/ecoli_operon_reference.parquet \
  --official-embeddings <official embed parquet> \
  --method-name <name> --model-revision <rev> --output <dir>
```

其余工件：`eval_bacformer/`、`eval_esm/`（summary/metrics/aligned_pairs/contract/completion）；embedding parquet 不进 Git。
