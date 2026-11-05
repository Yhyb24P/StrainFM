# BacBench 药敏标签可得性审计

负责人：Codex 当前会话。日期：2026-11-04。状态：已执行精确 accession 连接；不是独立评价数据集。源码/输入 SHA-256 见 `summary.json`，运行源码见 `code_snapshot/`。无 Git 提交。

```bash
PYTHONPATH=src python3 -m strainfm labels \
  --labels data/results/bacformer_source_audit_20261104/hf_amr_binary_labels.csv \
  --download-manifest /mnt/e/StrainFM/data/raw/download_manifest_all.tsv \
  --output data/results/bacbench_label_audit_rerun
```

公开标签来自 `macwiatrak/bacbench-antibiotic-resistance-protein-sequences`，revision `fa5f826a5a6ae28126483a751c6e8acffaca0eec`。文件 SHA-256：`644b91687d53c36e8645f0a5eb560b33cd807276ae41d81ae9a2f411050b0f63`。来源卡、URL 与抓取日期见相邻 source audit 的 `sources.json`。

- `matched_accessions.tsv`：26,052 个标签 accession 中精确匹配的 2,606 个。162 个同时匹配 genome_id/file_accession，2,444 个仅匹配 file_accession；不改前缀或版本。该表连接物种信息，不确认底层 FASTA 一致性。
- `species_label_coverage.tsv`：分母为每种精确匹配的唯一标签 accession；“至少一个标签”允许不同 accession 具有不同药物标签。
- `species_drug_counts.tsv`：每个物种×药物的非缺失标签计数、0/1 计数及 `n_label_1/n_observed`。没有使用模型指标选择药物；同一基因组的多药物观测不是独立基因组。
- `summary.json`：26,052 行、37 列药物、23,428 行至少一个标签、292,996 个已观测标签单元格。空值不填 0。

本轮没有拉取整套蛋白或 DNA 数据集，没有运行预测。只用于判断可得性，未完成 BioSample 去重、真实采样来源核验、实验药敏条件整理、预训练重叠或数据使用许可审查。精确连接漏掉的配对/历史版本可在官方身份核验后补充，不能直接进行字符串归一化。
