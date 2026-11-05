# StrainFM 数据治理（Phase 0 草案）

状态：2026-10-29 草案。数据治理服务于可复现的菌株级研究；未完成来源、权限和质量确认前，不把清单记录数当作可用训练样本数。

## 1. 当前清单事实

- 路径（用户提供）：`/mnt/c/Users/12087/Downloads/manifest_all.tsv`
- 大小：82,082,573 B（79 MB）；数据 554,459 行 + 1 行表头；`genome_id` 无重复。
- SHA-256：`cb4d06fe58c957e22d3c4685247a866f8adb03304f24541c0c4f0ad95d6961e7`。
- 权限：777，位于 Windows Downloads，非项目受控路径；建议复制到项目只读路径并登记哈希，下载目录不作为长期数据源。
- 表头：`genome_id`、`species_id`、`fasta`、`kingdom`。`kingdom` 取值为小写 `bacteria`/`virus`/`archaea`/`fungi`；`species_id` 是名称而非 TaxID，存在名称变体与复合群标签噪声。
- `fasta` 有两种路径模式：约 437,516 条为 `/fast/genomes/`，约 15,526 条为 `ncbi_dataset/data/`；101,408 条病毒 ID 非 GCF/GCA 且无版本号。
- 分域统计：
  - 细菌：443,876 条，1,864 个物种，7 个物种 ≥10,000 条；
  - 病毒：109,891 条，2,157 个物种；
  - 古菌：428 条，15 个物种；
  - 真菌：264 条，8 个物种。
- 细菌组装来源：GCF 381,837 条（RefSeq），GCA 62,039 条（GenBank），合计 443,876 条。
- 物种计数（与材料一致）：E. coli 62,560；K. pneumoniae 44,401；P. aeruginosa 17,579。
- `fasta` 字段指向 `/data2/herui/StrainTrace/Ref_download/info/genomes/ncbi_dataset/...`；该路径在本机不可读，尚未核验文件存在性、压缩完整性和校验和。
- 清单**缺少**：完整度/污染度、N50、contig 数、基因组长度、组装级别、BioSample/BioProject、样本来源、采集时间、测序平台、注释版本、物种 TaxID/GTDB 标签。不能仅凭该表判断质量或独立样本数。

## 2. 数据来源、权限与发布

- 清单及其 FASTA 为用户提供的项目输入；仅用于本项目研究，不对外重分发，不提交到公开仓库。
- 使用前确认：数据获取权限、机构数据保护要求、临床样本相关限制。
- 若使用 NCBI/AllTheBacteria 补充数据，记录 accession、下载日期、数据库版本和许可；去重后保留“来源数据库”字段。
- 镜像源：GWH 的 `/assembly/ncbi/` 是观察到的静态文件布局；公开 API 只覆盖 `GWH*` 记录，GCA/GCF 查询依赖浏览页内部接口，不作为稳定 API 或主下载源。GWH 仅用于人工/一次性补充，使用前核验 MD5；NCBI 为主队列。
- 发布代码/结果时只发布可复现脚本、聚合统计和必要元数据；不发布受限原始序列、临床敏感信息和权重。

## 3. 质量分层

Bacformer 的 MAG 阈值（完整度 >50%、污染 ≤10%、长度 250 kb–20 Mb、contig ≤1,000）服务于跨物种语料，不足以支撑菌株级 synteny/质粒研究。该阈值的原文页码待核验（v1 Methods 或 v3 Methods；引用前核对）；本项目采用分层策略：

| 层级 | 用途 | 操作阈值（草案，试点后校准） |
|---|---|---|
| A 完整/近完整基因组 | 主训练、synteny、质粒、基因邻接 | NCBI assembly level = Complete Genome 或 Chromosome；contig ≤ 50；N50 ≥ 1 Mb；长度 0.5–10 Mb；有完整注释 |
| B 高质量 draft | 主训练/稳健性 | assembly level = Scaffold/Contig；contig ≤ 500；N50 ≥ 50 kb；长度 0.5–10 Mb；若为 MAG，完整度 ≥90%、污染 ≤5% |
| C 低质量 MAG/draft | 仅稳健性/敏感性 | MAG 完整度 50–90% 或污染 5–10%；明确标记；不进入主结论；不与 A/B 混合调参 |
| 排除 | 不用于训练/评价 | 完整度 <50%、污染 >10%、嵌合、物种标签错误、无法确认来源 |

- 质量来源：NCBI assembly stats（assembly_level、contig_count、N50、长度）为主；CheckM/CheckM2 只对 MAG 子集抽样，不对 44 万 NCBI 组装全量运行；BUSCO/注释完整性按需抽样。阈值需在试点 100–200 个基因组后校准并版本化。
- 统一注释：Bakta 生成 CDS/蛋白/坐标；保留原始注释做敏感性对照。注释失败或碎片化严重的基因组单独标记，不把“注释缺失”当作真实基因缺失。
- 质量分层必须在 train/validation/test 划分之前完成，并在结果中分层报告。

## 4. 去重、谱系与独立划分

1. **精确去重**：`genome_id` 已无重复；仍需检查同一 BioSample/BioProject 的重复组装、不同版本（GCF/GCA）和跨库重复。
2. **近重复**：以 Mash/sourmash、core-SNP 距离或 PopPUNK 簇定义近重复；FastANI 仅作辅助（在 ≥99.9% 区间分辨力不足）。预注册阈值敏感性（99.9/99.95/99.99 或 SNP 距离），报告各阈值下保留数量；每个近克隆簇最多保留一个代表进入训练，测试集按簇留出。
3. **谱系**：PopPUNK/PopPIPE、cgMLST（chewBBACA）或 core-SNP 树；划分按谱系/克隆簇不相交，避免同一暴发克隆同时出现在训练和测试。
4. **来源留出**：按国家/研究/BioProject/测序平台/时间留出独立来源，作为 H2/H4 的主检验。清单缺这些字段，需从 NCBI Datasets 补齐并先审计缺失率。
5. **预训练重叠**：48,932 个完整基因组可按 accession 匹配；1.3M MAG 的逐基因组清单未必公开，需向作者申请，否则用 Mash screen 近似并报告检测限。重叠样本从**全部评价集**（不只 test）排除并报告比例；ESM-2/UniRef 预训练污染无法排除，作为解释限制。
6. **划分记录**：每个基因组记录 `split`、`dedup_cluster`、`lineage`、`source_group`、`quality_tier`，固定随机种子与划分脚本版本。

## 5. 元数据字段（最低集）

| 类别 | 字段 |
|---|---|
| 身份 | genome_id、assembly accession、BioSample、BioProject、物种原始名、统一名、TaxID/GTDB |
| 质量 | 完整度、污染度、N50、contig 数、基因组长度、组装级别、BUSCO/CheckM 版本 |
| 来源 | 数据库（NCBI/AllTheBacteria/自建）、下载日期、许可、国家/研究/平台/时间 |
| 注释 | 注释工具与版本、CDS/蛋白数、假基因数、注释日期 |
| 生物学 | 核心/accessory 基因数、质粒/岛、MLST/cgMLST/PopPUNK 谱系 |
| 训练 | 去重簇、quality_tier、split、是否与预训练重叠、是否进入训练 |

字段缺失要显式记录为 unknown，不用默认值填充；任何填充规则写入数据处理日志。

## 6. 版本、校验与审计

- 清单：记录文件 SHA-256、行数、生成日期、生成脚本/来源；后续更新保留旧版本。
- FASTA：记录文件路径、大小、MD5/SHA-256、下载日期、来源 accession；抽样校验 gzip 完整性与序列数。
- 注释/模型/权重：记录版本、哈希、配置文件、运行命令、随机种子。
- 结果：每个结果目录包含 `README`（数据版本、代码提交、参数、指标定义、负责人、日期）。

## 7. 当前缺口与下一步

**缺口**：
- FASTA 在本机不可读；远程路径权限/存在性/校验和未知；
- 无质量、来源、组装级别、样本元数据；
- 无去重/近克隆/谱系划分；
- 无预训练重叠检查；
- 目标物种与用途未定。

**下一步（Phase 0）**：
1. 确认 FASTA 访问方式与数据权限；
2. 生成清单 SHA-256 与数据字典；
3. 对候选物种抽样 100–200 个基因组做质量/注释/去重试点；
4. 建立 `data/raw`（只读）、`data/processed`、`data/results` 目录；
5. 把质量分层、去重、划分脚本写入正式实现；
6. 在导师确认物种与用途后，记录数据权限、来源与可用范围。
