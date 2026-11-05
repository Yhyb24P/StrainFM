# NCBI 抽样组装元数据审计

负责人：Codex 当前会话。日期：2026-11-04。状态：100/100 个请求 accession 返回报告。无 Git 提交，执行源码、输入 SHA-256 见 `summary.json` 和 `code_snapshot/`。

```bash
PYTHONPATH=src python3 -m strainfm metadata \
  --sample data/results/manifest_local_audit_20261104/fasta_sample.tsv \
  --output data/results/ncbi_sample_metadata_rerun
```

仅发送样本表里的公开 versioned assembly accession，使用 NCBI Datasets v2 官方 API，每批 20 个，顺序查询，未上传序列或个人元数据。API 内容会随时间变化，复核本次结论应使用保存的响应而非假定最新响应相同。

- `ncbi_batch_*.json`：官方 JSON 原始响应；仅本地保存，其中存在公开 BioSample 来源描述，不默认对外发布。
- `requests.json`：精确 URL、accession 列表、UTC 抓取时间及响应 SHA-256。
- `sample_metadata.tsv`：本地样本与官方报告的连接、组装级别、BioProject/BioSample、来源字段有无、真实 contig/scaffold 统计、原注释蛋白编码基因数。
- `summary.json`：按物种计数。各类 known 字段的分母为该物种返回的报告数，本轮均为 25；常见缺失占位值被排除，未核验字段语义或不常见属性名称。

100/100 个本地 FASTA 总长度与官方 accession 报告一致；100/100 个原 genome_id 等于该 accession 或 paired accession。长度相同和配对成立都不等于序列完全一致。一个 BioProject 也不自动成为独立研究。

蛋白数来自当前公开注释，不是统一 Bakta 结果。P. aeruginosa 的 25 个报告有 16 个蛋白数已知，10/16 超过 6000，9 个未知；不能用 10/25 报告已知样本中的超限比例。E. coli、K. pneumoniae 各 20 个已知，B. subtilis 22 个已知，已知者均不超过 6000。

FASTA 记录可能为 scaffold，因此不能把上一审计的记录数/记录 N50 直接写成 contig 数/N50。本目录保留 NCBI 的两组统计供区分。
