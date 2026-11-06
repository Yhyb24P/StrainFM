# E. coli 候选数据准备与划分

日期：2026-11-05。负责人：当前 Codex 执行会话。状态：候选数据的全量文件快照、元数据、身份合并、Mash 多样性分组及来源隔离划分均已完成并核验。**E. coli 正式目标物种选择仍待用户确认**，因此保留候选状态；没有开展注释、蛋白编码、模型预处理、推理或训练。

## 第一阶段结果

| 对象 | 数量 | 定义 |
|---|---:|---|
| 原始映射记录 | 62,560 | 下载清单内 E. coli 记录 |
| 完整读取且格式有效文件 | 62,555 | 完整解压与 DNA IUPAC 检查通过；另 5 个缺失 |
| 精确 FASTA 记录多重集合 | 62,060 | 对全部有效文件忽略 header、换行、大小写和记录顺序后计数 |
| 满足当前候选纳入条件的文件 | 60,984 | 使用下述组装统计及状态条件 |
| 身份代表 | 60,416 | 合并同一 BioSample 或精确记录多重集合副本 |
| 近克隆操作组与保留代表 | 30,687 | Mash 距离 ≤0.001 的边形成单连通分量，每组保留一个质量优先代表 |
| 较粗序列邻域 | 2,125 | Mash 距离 ≤0.005 的边形成单连通分量，定义于 60,416 个身份代表 |

以上数量逐层口径不同，不能把文件数、BioSample 数或操作组数直接称为统计独立菌株数量。近克隆最大组有 3,335 个身份代表；较粗邻域最大组有 9,003 个，均未把组内所有两点距离解释为阈值以内。

全量压缩输入为 94,859,610,772 字节。原始文件保持原位只读使用，没有复制整套序列或修改共享权限。每个纳入/排除记录和内容哈希保存在 `data/processed/escherichia_coli_download_20261105/assemblies.tsv`。

## 候选开发集与独立测试集

| 集合 | 基因组代表 | 占全部代表 | 完整成员涉及的不同 BioProject | 较粗序列邻域 | 有训练参考库正例的查询 |
|---|---:|---:|---:|---:|---:|
| 训练 | 26,462 | 86.23% | 3,480 | 1,833 | — |
| 验证 | 2,113 | 6.89% | 1,469 | 452 | 1,953 |
| 测试 | 2,112 | 6.88% | 1,496 | 480 | 1,959 |

开发集为训练+验证，共 28,575 个代表；测试为 2,112 个代表。来源项目数包含被合并身份/近克隆成员的全部已知来源，不只统计最终代表自己的项目。

全体成员来源关系形成 2,520 个来源连接分量，最大分量含 26,462 个保留代表，必须整体进入一个集合，超过 80% 的目标比例。因此实际划分约 86.2/6.9/6.9，优先保证完整来源隔离；未通过拆组来凑比例。固定 seed=20261105。

验证有 160 个查询、测试有 153 个查询在训练参考库中没有同邻域正例，显式保留为无正例/开放集合查询，不填零混入正例支持检索平均。测试覆盖 480 个较粗邻域，其中 331 个在训练中可见、149 个不在训练中可见。这里的“邻域”不是经 MLST 或系统发育验证的谱系；来源留出也不自动等于所有谱系留出。

检索正例定义为同一 Mash 较粗序列邻域，用于种内遗传结构诊断，不代表独立表型或传播标签。基础模型预训练成员关系仍为 unknown；独立性已验证的范围是当前候选适配/开发数据的身份、近克隆和已知来源隔离，不能写成基础模型确定从未见过这些基因组。

### 实际划分的质量覆盖

下表 contig 统计来自对应 accession 的 NCBI 元数据，FASTA 长度包含模糊字符。

| 集合 | Complete Genome / Chromosome | Scaffold / Contig | 长度中位数 bp | contig 数中位数 | contig N50 中位数 bp |
|---|---:|---:|---:|---:|---:|
| 训练 | 3,125 / 244 | 7,239 / 15,854 | 5,055,711.5 | 127 | 158,764 |
| 验证 | 447 / 25 | 435 / 1,206 | 5,084,721 | 115 | 178,671 |
| 测试 | 460 / 42 | 451 / 1,159 | 5,052,048 | 104.5 | 181,906 |

测试中的完整/染色体级组装占比高于训练，已保留这一质量差异；第二阶段应分层解释，不将未匹配的质量组成当作模型能力差异。这里没有重新估计污染或完整度。

## 元数据与纳入规则

62,560 个实际带版本 accession 以每批不超过 100 个、2 个并发请求查询 NCBI Datasets v2。626 批无批次失败，62,554 个取得精确版本报告，6 个标记 unknown，不用当前/配对版本静默替代。元数据表 SHA-256：`5bf3b6f4c619226b586676aa6d7924d5b1c7018be5b69d419a88067a22abc33e`。

以全部请求 accession 为分母，BioSample 和 BioProject 各覆盖 62,554 条，采集日期 56,850 条，地理字段 59,120 条，分离来源 45,174 条，contig N50 和组装级别各 62,554 条。字段存在不保证语义正确。全部返回元数据涉及 62,456 个不同 BioSample、6,845 个 BioProject；同一 BioSample 多个 accession 的超额记录为 98 条。最大项目有 2,052 条，前十个项目合计约占返回报告的 11.98%。见 [metadata_coverage.json](metadata_coverage.json)。

候选纳入使用：有效文件、精确版本报告、匹配物种名、非 suppressed、未标为 taxonomy check failed、长度 0.5–10 Mb、NCBI contig 数 ≤500、contig N50 ≥50 kb、模糊字符比例 ≤1%。参数见 `config/escherichia_coli_cohort.json`，这是操作性筛选，不是生物质量认证。

共排除 1,576 条记录。原因可重叠：N50 不足 697、suppressed 463、模糊字符过高 415、长度超范围 68、taxonomy check failed 62、contig 数过高 39、精确元数据未知 6、文件缺失 5。不能直接相加原因计数。所有纳入样本有已知 BioSample 和 BioProject。

## 版本与本地交付物

运行基于 Git 根提交 `340e1cb96576782d1fe742143bc8a2fbee35c464` 加本轮未提交实现，不能仅用该提交复现实验；实际代码哈希和源码快照分别保存在各输出的 contract/summary 与 `code_snapshot/`。正式实现统一在 `src/strainfm/`。原始序列、原样元数据和大规模索引保留本地，不自动对外发布。

| 交付物 | 本地路径或摘要 |
|---|---|
| 文件版本与哈希索引 | `data/processed/escherichia_coli_download_20261105/assemblies.tsv`；[摘要](file_snapshot_summary.json) |
| 精确版本元数据与请求来源 | `data/processed/escherichia_coli_ncbi_metadata_20261105/` |
| 全体纳入/排除与身份组成员 | `data/processed/escherichia_coli_identity_cohort_20261105/cohort_manifest.tsv`；[摘要](identity_cohort_summary.json) |
| 序列草图、分组与生成森林 | `data/processed/escherichia_coli_mash_groups_20261105/`；[摘要](sequence_grouping_summary.json) |
| 固定候选划分 | `data/processed/escherichia_coli_source_split_20261105/split_manifest.tsv`；[摘要](candidate_split_summary.json) |
| 每个验证/测试查询的正例支持 | `data/processed/escherichia_coli_source_split_20261105/retrieval_query_support.tsv` |
| 覆盖全体输出的最终一致性检查 | [candidate_delivery_verification.json](candidate_delivery_verification.json) |
| 一次针对性实现复核与修复 | [implementation_review.json](implementation_review.json) |

下载清单 SHA-256：`a96528c903f7c8422ecb6a2e3b399b24edd1aaf8d0555b1313e4aa55c8637bb7`。

文件索引 SHA-256：`1babfba1bb0c89a36bcddc1a1b19796cda330909ffa9191b5f48be509ad27451`。

分组表 SHA-256：`b509a575ccd2a08832b67f881570c13a537117e19dc092e9c8d44c74b91bef05`。

候选划分 SHA-256：`34ba67134a67b416d342a6449c2019cd09c353dfcfb759bf7c29a8b02b9b2a59`。

## 复现入口

在项目根目录执行。大输出保留在被 Git 忽略的 `data/processed/`。快照和元数据可按相同契约断点恢复；Mash 草图批次可复用，但距离阶段再次执行会重新计算，因此不要重跑已完成的 diversity 命令。正式复核请使用新输出目录；`cohort`、`split` 会拒绝已有目录，不能静默重写固定划分。

```bash
PYTHONPATH=src python3 -m strainfm snapshot \
  --config config/escherichia_coli_snapshot.json \
  --output data/processed/escherichia_coli_download_20261105 --workers 8

PYTHONPATH=src python3 -m strainfm metadata-bulk \
  --inventory data/processed/escherichia_coli_download_20261105/inventory.json \
  --output data/processed/escherichia_coli_ncbi_metadata_20261105

PYTHONPATH=src python3 -m strainfm cohort \
  --snapshot data/processed/escherichia_coli_download_20261105/assemblies.tsv \
  --metadata data/processed/escherichia_coli_ncbi_metadata_20261105/assembly_metadata.tsv \
  --config config/escherichia_coli_cohort.json \
  --output data/processed/escherichia_coli_identity_cohort_20261105

PYTHONPATH=src python3 -m strainfm diversity \
  --representatives data/processed/escherichia_coli_identity_cohort_20261105/identity_representatives.tsv \
  --config config/escherichia_coli_mash.json \
  --output data/processed/escherichia_coli_mash_groups_20261105 --threads 8

PYTHONPATH=src python3 -m strainfm split \
  --representatives data/processed/escherichia_coli_identity_cohort_20261105/identity_representatives.tsv \
  --groups data/processed/escherichia_coli_mash_groups_20261105/sequence_groups.tsv \
  --cohort-manifest data/processed/escherichia_coli_identity_cohort_20261105/cohort_manifest.tsv \
  --output data/processed/escherichia_coli_source_split_20261105
```

Mash 2.3 二进制 URL 与归档/二进制哈希在本机 `~/.local/share/strainfm/tools/mash-2.3/source.json`、运行 contract 和吞吐记录中。配置为 k=21、s=1000、seed=42，分批草图后计算全部下三角比较，共 1,825,016,320 对；筛出较粗阈值内关系 19,460,122 条，近克隆阈值内 962,684 条。只存储能重建相同分组的生成森林：近克隆 29,729 条、较粗邻域 58,291 条。Mash 是 MinHash 距离估计，不是 core-SNP 或精确 ANI，单连通分组也不保证组内任意两点均满足阈值。

## 核验和保留的失败

12 项测试通过，包含真实 Mash 的重复/无关序列、损坏与缺失文件、断点恢复契约、精确版本连接、身份合并和来源隔离。

正式实现的一次只读复核发现：只保留代表自己的项目，会丢失被合并成员的来源约束。已修复为划分时读取完整 cohort manifest，传播身份组及近克隆组所有来源并按来源连接分量分配；回归同时覆盖被移除身份成员和近克隆成员来自不同项目的情况，修复经同伴核对通过，无需重跑 Mash。

自动排队的旧 split 调用缺少修复后新增的必需 `--cohort-manifest` 参数，以 argparse 退出码 2 停止，未写出划分；补齐参数后正式命令成功完成。没有因观察超时而重启数据处理，也没有覆写旧划分。

最终检查从保存的两个生成森林重建了全部 60,416 个分组，逐一核对了原始合格成员的来源和 BioSample 无跨集合、全部查询的正例数和实际划分的质量/邻域覆盖。它验证数据交付完整性，不是模型科研结果。

预训练来源仅做一次有界核查：官方仓库固定版本只提供训练样本文件；公开的 `bacformer-genome-embeddings-corpus` 数据卡说明这些基因组由已预训练模型编码，不能据此判定逐 checkpoint 的预训练成员关系，数据 split 名称 test 也不能证明未见。因此保留 exposure unknown；来源、哈希和原始响应在本地 `data/processed/escherichia_coli_upstream_exposure_20261105/`，未下载整套嵌入或联系作者。

**剩余决定**：由用户/导师确认是否正式采用 E. coli。本阶段已把候选数据版本、数量、多样性和划分做成可审阅产物，不以更多审计或模型预处理替代该决定。
