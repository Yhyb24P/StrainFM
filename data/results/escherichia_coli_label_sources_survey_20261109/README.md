# 大肠杆菌标签扩展与三方法共同评价：资源实测

日期：2026-11-09。负责人：gpt。状态：**已完成资源查询和候选计数，尚未形成正式训练集合或模型结果**。代码基线 `0d93927`；用户附件为 `/home/yhshy/.codex/attachments/e903cf80-479f-4779-8993-cf8565780398/pasted-text-1.txt`。

本轮回答：能否扩大“同一菌株的基因组—实验标签”集合，以冻结 Bacformer、ESM++ 均值及 StrainAMR 做共同评价？已有 AMR-02 结果保留；不把数据源记录数、模型预测或物种汇总当成独立菌株实验真值。

## 结论与下一步

**AMR 扩样有具体数据支撑，主要增量来自带标签但本地缺序列的菌株。非 AMR 任务有候选，但基因组链接覆盖明显较低。**

1. BV-BRC 两药共得到 30,961 条 `Laboratory Method` 记录，12,217 个 genome ID。完成 genome 元数据查询和 1,905 个 ENA 样本别名转换后，有 7,814 个不同标准 BioSample 标识；其中两药至少一药有非冲突 S/R 的并集为 **3,945**，两药交集为 **3,444**。这是标识去重后的候选规模，不是近克隆独立样本数。
2. 这 3,945 个候选中 **3,351** 个没有可直接复用的本地质量合格组装匹配。先补这些候选的序列/身份，不扩展无标签全量编码。两条 accession 与 BioSample 对应冲突已阻止自动本地匹配。
3. MetaTraits 的物种汇总只用于找性状；进一步从上游 BacDive 实际取得 **1,975** 个 E. coli 条目，只有 **222** 个含基因组链接。可优先核查鸟氨酸脱羧酶或 α-葡萄糖苷酶的菌株身份、实验条件和有效分组支持。
4. 三方法比较可沿已有冻结编码流程扩展，但必须先形成共同样本/划分；StrainAMR 的特征构建需严格区分训练、验证、测试。另用固定验证/测试集的嵌套训练子集实验判断增加训练样本是否有效。

## BV-BRC：已下载的实际证据

本次查询范围固定为 `taxon_id=562`、`public=true`、`evidence=Laboratory Method`、药物名精确匹配 `ciprofloxacin` 或 `gentamicin`。**不声称覆盖 E. coli 所有下级 TaxID、所有证据类型或所有药物。** 官方 schema 的示例值不能替代实际字段检查；本次响应的实验标记实际为 `Laboratory Method`。[AMR 数据接口](https://www.bv-brc.org/api/doc/genome_amr)

| 药物 | AST 记录 | 标准 BioSample | S | R | I | 无 S/R/I 判定 | 非空判定冲突 | 有本地合格候选的 S/R |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ciprofloxacin | 16,741 | 7,378 | 2,084 | 1,596 | 59 | 3,611 | 28 | 261 / 318 |
| gentamicin | 14,220 | 7,622 | 2,787 | 922 | 99 | 3,811 | 3 | 375 / 103 |

分母与规则：

- 每药按标准 BioSample 合并所有命中记录；缺失判定不投票，有 S 或 R 而另有未知记录不自动构成冲突。两个以上不同的非空判定记冲突（包括 S/I、R/I，不仅 S/R）。I 单列，不并入 S 或 R；只给 MIC/抑菌圈而没有判定的记录不自行套断点。候选记录保留原实验方法、标准、年份、测量及引用。
- 两药分别还有 4,400 / 4,366 个缺标准 BioSample 的 genome ID，暂按未解析身份保留，**不计入表中标准 BioSample 数**。可用 SRA 标识进一步连到 ENA，但尚未批量完成。
- genome 元数据中 7,816 个非空 `biosample_accession` 字段含 1,905 个不同 `ERS…` 别名；实际 ENA XML 的 `EXTERNAL_ID namespace="BioSample"` 将其连到 `SAMEA…`。按原字符串计数或直接拼接 ATB URL 都不可靠。转换后合并重复标识。[ENA XML 文档](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html)
- 本地比较使用现有 `cohort_manifest.tsv` 的精确版本 accession（file/current/paired 分开保留）及 BioSample。若 accession 与 BioSample 指向不同本地样本，停止自动链接。匹配到本地质量合格组装的 **781 是 BV-BRC genome ID 数**，不能当作新增独立样本数。本轮复用既有质量结论，未重新全量哈希本地 FASTA。
- 实验标记仍需解释条件：两药缺实验方法 62 / 47 条，缺 testing standard 7,548 / 4,767 条，缺标准年份 11,641 / 8,761 条；这些是**记录级**缺失数。候选不等于已逐条验证实验可信度，也不能按记录数增加统计样本量。

来源下载按 `sort(+id)&limit(5000,start)` 分页，校验每页物种、药物、证据、public、总数及全体 record ID 唯一性。12,217 个 genome ID 分 62 批查询指定身份/质量字段，全部返回。每页 URL、时间和 SHA-256 均本地保存；本次是查询响应快照，不是上游提供的不可变 release。

**同伴复核后补充的下载前分解**：3,351 个无本地合格匹配候选中，3,329 个在 BV-BRC 根本未报告 assembly accession，19 个报告的 accession 不在本地冻结 cohort，2 个身份冲突，1 个 BioSample 在本地但无合格匹配；仅版本主体命中的情形为 0。类别按身份冲突优先、BioSample、精确 accession、accession 主体依次判断，互斥且和为 3,351。这不是磁盘全量搜索，也不意味着 3,329 个序列不存在；应通过 BioSample/SRA 和 ATB 找对应组装，而非批量把空 accession 当下载失败。见 `missing_local_breakdown.json`。

另将已发布 BacBench 二元标签经本地精确版本 file/current/paired accession 连到 BioSample，与本次 BV-BRC 非冲突 S/R 比较：交集环丙沙星 58、庆大霉素 57，当前交集判定全部一致。仅对交集成立，不代表两个数据库整体等价；未自动合并标签或改写旧队列。未来不一致记录须保留两来源并单独裁决。见 `cross_source_labels_summary.json`。

## MetaTraits / BacDive：非 AMR 性状可行性

[MetaTraits 文档](https://metatraits.embl.de/documentation)区分 culture-based 与 prediction-based 数据。物种 API 请求本次返回 HTML 404；不能解释为 E. coli 没有标签。网站及 culture-only 物种汇总文件可读，后者 E. coli 有 133 行。**汇总缺菌株 ID，不能展开成每个基因组的标签**；汇总观察数与独立条目数不同。

沿汇总的来源追到 [BacDive v2 API](https://api.bacdive.dsmz.de/)，实际分页取回全部 1,975 个条目。本次只解析 `Physiology and metabolism` 中 enzymes 的 `activity` 及 metabolite utilization 的 fermentation；完全不读取 `Genome-based predictions` 作为真值。`+/-` 分别为阳/阴，同一 BacDive ID 的同名性状冲突单列；保留 `@ref` 回溯来源，包括由实验 API 生化面板转录的注释。还需核对面板条件及保藏号重复。

| 候选性状 | 条目阳/阴 | 带基因组链接的阳/阴 | 本地 accession 主体命中候选阳/阴 |
|---|---:|---:|---:|
| 鸟氨酸脱羧酶 | 306 / 847 | 74 / 61 | 22 / 37 |
| α-葡萄糖苷酶 | 195 / 210 | 61 / 38 | 14 / 16 |
| α-半乳糖苷酶 | 270 / 138 | 61 / 41 | 22 / 10 |
| β-葡萄糖醛酸苷酶 | 269 / 141 | 56 / 47 | 25 / 8 |
| 蜜二糖发酵 | 321 / 428 | 46 / 5 | 39 / 4 |

这里的单位是 **BacDive 条目**，尚未按 BioSample/保藏号/序列去重；表中阳/阴不含冲突条目。完整冲突和覆盖计数见 `bacdive_trait_summary.tsv`。基因组链接常为无版本 GCA，末列仅是 accession 主体候选匹配，**不能直接转标签**。BacDive 的 genome score 是组装质量排序，不是身份匹配置信度。[基因组分数定义](https://hub.dsmz.de/wiki/bacdive/genomes/)

因此优先核实 135 个鸟氨酸脱羧酶二元条目的准确组装、实验来源和独立分组，再判断是否足以评价；不因其条目数最多就承诺正式任务。StrainAMR 仅用于 AMR 对照，不强行用于这些酶/代谢任务。

## AllTheBacteria：补序列，不能补实验真值

官网当前文档指向 r0.2 + 2024-08 + 2025-05 的聚合数据；按 `sample_accession` 可从 AWS 单样本下载，并取得组装、物种和 QC 元数据。[元数据](https://allthebacteria.org/docs/sample_metadata/)、[组装下载](https://allthebacteria.org/docs/assemblies/)。补样本时记录 release、下载日期、样本/测序 run、文件哈希；ATB 同 BioSample 的重组装不当作另一独立菌株。

实测：从“有标准 BioSample、有至少一药非冲突 S/R、没有本地合格匹配”的候选中，按 BioSample 字典序取前 10 个，AWS HEAD **10/10 返回 200**，记录长度与 ETag；不是随机抽样、不是全候选可用率，也没有下载或验证 FASTA。此前用 ERS 别名试探的 10 次 404 保留为失败证据，不能解释成 ATB 缺样本。候选里也包含仅因身份冲突而暂不能复用本地文件者，应先解冲突再下载。

ATB 的 AMR 输出是 AMRFinderPlus 检测，不当作实验药敏标签。[AMR 注释来源](https://allthebacteria.org/docs/amr/)。ATB HQ 标记也不自动满足本项目原有基因组上下文质量要求。

## StrainAMR：固定版本后的接入边界

本次 `v1.0` 解析为 commit **`c57619eed0fd7121f2b475f5f6401bf679f629f9`**；已保存 README、四个入口和 requirements 的本地源码快照及哈希。[固定版 README](https://github.com/liaoherui/StrainAMR/blob/c57619eed0fd7121f2b475f5f6401bf679f629f9/README.md)

- 方法：变异图 SNV、蛋白簇与筛选 k-mer 三通道，经多模态 Transformer 得分；当前版本还提供 AMRFinderPlus/CARD 机制层与 hybrid 输出。
- 共同评测使用相同菌株、DNA、标签和外部冻结划分；ESM++/Bacformer 另须相同蛋白输入。StrainAMR 可以采用其规定的内部预处理，作为专项系统比较，不把差异全归为表示能力。
- `build-train` **仅看训练集**，其特征筛选、蛋白簇、词表、图、机制校准均不可先看全队列。验证与测试分别 `build-test --db train_db`；训练明确使用 `--val-dir validation_features`，最终测试不出现在 `train` 的输入中。README 示例把验证目录称作 test_features，不能照名称把真正测试集交给选模。[训练入口](https://github.com/liaoherui/StrainAMR/blob/c57619eed0fd7121f2b475f5f6401bf679f629f9/StrainAMR_model_train.py)
- 默认包含机制层；`Model_prob` 的 AUROC/AP 与 hybrid 的分类结果分开报告。规则阈值只用验证集选择，不能把二元 hybrid call 当作连续分数进行同等排名比较。
- 代码含 `setup_seed(10)`，并未提供已核实的三种子 CLI 契约；正式重复运行前需显式处理并记录，不能声称重复命令就是三个独立随机种子。
- 需要 RGI/CARD、MMseqs/CD-HIT、minimap2、seqwish、Mash，机制层另需 AMRFinderPlus 及数据库；CPU fallback 在代码中存在，但本轮未安装或跑通，不承诺 CPU 耗时/内存。资源探测时本机为 RTX 5090 D（约 32 GB），E: 可用约 905 GB；不是长任务资源预约或预算。
- README 明写许可待定、未附 LICENSE；当前仅本地研究源码，不重分发代码/权重。离线执行需先固定所有数据库和依赖，本轮未验证离线运行。

## 下一轮共同实验的具体顺序（方案，未执行）

1. 把本轮候选表导入现有单一 CLI：保留原记录、身份证据与冲突，解决必要的 SRA→BioSample、accession 版本、实验来源链接；已有文件复用，仅补确有实验标签的缺失组装。先完成上述两个药物，不同时铺开所有药物。
2. 在模型测试表现不可见时，固定新增队列和三方法共同划分。BioSample 与近克隆整组隔离；报告 BioProject/研究来源支持，若来源隔离可行则留出来源，不自动沿用旧队列的来源共享例外。不同生物样本的近克隆 S/R 差异不是记录错误，与同一样本实验冲突分开处理。
3. 旧 AMR-02 测试已用于探索解释，应作为已见历史证据；新确认性测试优先来自新样本/组。训练和验证也扩充，不能只换测试集后把差异归为训练数据量。统一跨药样本分区便于以后多任务扩展。
4. 在同一固定验证/测试上预先设定嵌套的训练组子集（例如可用训练组的 25%、50%、100%，按组取整并检查两类支持），对三方法分别重新训练各自下游组件。ESM/Bacformer 主干冻结，固定编码/预处理；StrainAMR 每个数据规模重新构建训练特征，不能复用全训练规模选出的特征。训练组数量之外，还报告菌株数与两类数。
5. ESM/Bacformer 轻头协议一致；StrainAMR 按固定官方流程训练并记录不同的特征资源、调参预算、CPU/GPU 时长。共同报告每样本分数、AUROC/AP、验证集选择阈值下的分类指标及配对组级不确定性。先验证数据规模是否改善固定测试上的泛化，再决定池化、物种适配或新架构；非 AMR 任务独立判定可行性。

## 文件、复查与协作

- `bvbrc_counts.tsv`、`summary.json`：已查询范围的聚合计数。
- `bacdive_trait_summary.tsv`：本次解析范围内正负各至少 20 条的性状候选；阈值仅用于展示，不是正式任务纳入标准。
- `artifacts_sha256.json`：本地响应、映射与候选表的哈希索引。记录级数据位于忽略目录 `data/processed/ecoli_external_sources_20261109/`，包括 `sample_drug_candidates.tsv`、`identity_matches.tsv`、`bacdive_trait_candidates.tsv`，未发布。
- `verification.json`：计算一致性复核；`collaboration.json`：注册、协调与一次只读同伴复核记录。ds 在 `UPD-ds-20261109-label-survey-review-044` 核对原 275 个文件哈希、计数及方法边界，未独立访问上游或运行 StrainAMR。上面下载前分解及跨源交集为收到反馈后的补充，已本地校验，未再要求第二轮全面复核。
- 这是探索性资源调研，查询与计数规则已写明；正式生产入口尚未整合进 CLI，不能把本目录称为已冻结可复现训练数据集。既有源码、规范、划分、模型预测未改动。
- AIA：用户要求恢复 `gpt`；join 按 CLI 提示显式 reconnect，返回 `ready:true`、沿用 agent ID `agent-fc213dcd0cf14b568019df325f3da966`。实际 scope 为仓库 `.state/aia`，与 AGENTS 描述的主 scope 不一致；没有替换会话或擅改规范。已发 `UPDATE-gpt-20261109-label-expansion-start`，收到 ds `UPD-ds-20261109-label-survey-sync-043` 确认没有重复工作并记录读取回执。恢复接收不等于恢复所有旧角色权限。

本轮未做全量新序列下载、新划分、编码或训练；未得到新的模型优劣结论。候选规模仍会因身份核实、实验条件、质量和分组支持而减少。
