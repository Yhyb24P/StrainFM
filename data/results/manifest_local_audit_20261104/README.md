# 本地清单与 FASTA 抽样审计

负责人：Codex 当前会话。日期：2026-11-04。状态：已执行并核验；不是训练就绪或生物质量认证。无 Git 提交，输入、配置、源码哈希在 `summary.json`，实际运行源码见 `code_snapshot/`。

正式入口（项目根目录，Python 3.12.3，无第三方依赖）：

```bash
PYTHONPATH=src python3 -m strainfm audit \
  --config config/manifest_audit.json \
  --output data/results/manifest_local_audit_rerun
```

输出目录必须尚不存在；复核时使用不同目录。当前 CLI 后续增加了 labels/metadata 子命令，audit 的计算实现哈希与本次运行一致。

输入为 `/mnt/e/StrainFM/data/raw/manifest_all.tsv`、`download_manifest_all.tsv` 和本地 `genomes/`；仅只读访问，没有重命名、修复、删除或写回原始数据。

- `manifest_species_counts.tsv`：按原始 kingdom/species_id 逐行计数；原清单 554,459 行，不含表头。
- `species_inventory.tsv`：四个物种的全部映射行和目录非空文件覆盖，共 125,942 个映射记录、125,936 个非空文件。文件存在不是格式、来源身份或许可确认。
- `fasta_sample.tsv`：按 seed=20261104 的 SHA-256 排名，每种取 25 条；选样早于存在性检查。100/100 完整解压且格式有效。逐个记录压缩文件 SHA-256；失败不替换。

指标口径：`sequence_count` 为 FASTA 记录数；`total_bases` 含 N 与其他 IUPAC 模糊字符；`n50_bases` 为记录长度排序后的 N50，分母是总碱基数；`gc_fraction_acgt` 分母仅为 A/C/G/T；`ambiguous_fraction` 分母为全部碱基。`max_contig_bases` 和 `contig_multiset_sha256` 在此次输出中也按 FASTA 记录计算，记录可能是 scaffold，不能视为 NCBI contig 统计。真正 contig 数和 contig N50 见相邻 NCBI 元数据审计。

序列多重集合哈希忽略大小写、换行、header 和记录顺序，但保留重复次数，不做 reverse-complement 或环状起点归一。样本内重复超额数为 0，不代表总体无重复或无近克隆。

首次严格解析因第 453,053 行的额外空尾列停止，未生成完整结果。修订后对空尾列计数并允许读取，对非空多余列继续拒绝。原清单有 101,408 行带额外空尾列，四个正式字段均非空。失败记录保留在这里，原清单没有改动。

未完成：上游 checksum 比对、每个文件 accession 身份验证、近克隆、谱系、统一注释、污染/完整度与数据权限核查。抽样规模是每种 25，不能外推为全体质量保证。
