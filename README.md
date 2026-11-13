# StrainFM

当前主线：评估 Bacformer 在大肠杆菌真实标签下游任务上的能力，再依据证据决定模型研究切入点。

## 项目定位与当前状态

长期目标是单物种菌株级基础模型。近期比较冻结 Bacformer 与匹配 ESM 表示在可靠标签上的预测收益；表示、功能预测、检索与蛋白序列生成分别取证。

DH5α operon 官方输入和本地输入的评价闭环已完成，保留为单参考菌株结果，不重复开展，也不外推多菌株泛化。数据快照、质量检查、身份与 Mash 分组继续作为支撑；停止扩展无标签检索和全量编码，保留已落盘产物。

最新已完成结果为 2026-11-12 扩大队列 AMR 评测：冻结 Bacformer-L 相对匹配 ESM++ 均值在两药上有收益；环丙沙星另有 StrainAMR 发布模型对照。来源组件跨分区，外部模型未按共同训练集重训，不能据此声称来源独立泛化或受控的三模型架构优劣。历史小队列的负结果与敏感性分析均保留。

2026-11-13 归整发现并由协作者修正了外部对照的跨药索引错误；[更正结果](data/results/escherichia_coli_ast_amr_external_eval_20261113/README.md)已纳入 HTML。原 11 月 12 日外部表保留为失效记录。更正版统一采用逐种子指标平均，433 例子集的 StrainAMR AP 高于 Bacformer，但仍不是共同训练预算下的比较。独立分药的 `amr-eval` 主分析不受该索引错误影响。

直接打开[模型与历次结果 HTML](data/results/index.html)：按运行、任务、分析集筛选，查看架构、重点候选与原始表格。阅读[参考文献 HTML](ref/index.html)可查看全部收录资料的简要讲解、版本及适用范围。它们都可离线打开；原始来源链接需要联网。

物种适配、数据规模与生成评价仍是后续研究方向。研究计划中的 AMR-01 段落保留为当时阶段安排；实际后续执行依据 [DECISIONS](docs/DECISIONS.md) 的 D0-30/D0-31 与对应结果目录，不能把旧阶段状态当作最新完成状态。

## 目录与报告

`ref/` 保存书目、文献讲解与可获取全文；`data/results/` 保留按运行命名的原始结果，并以统一 [目录说明](data/results/README.md) 和 HTML 索引组织。`data/raw/`、`data/processed/`、`config/`、`src/` 分别承载输入、派生数据、配置和正式实现。历史目录不重命名，避免破坏结果契约、哈希和既有引用。

更新结果后，在仓库根目录运行：

```bash
PYTHONPATH=src python3 -m strainfm report
```

该命令只读取聚合结果并生成 `data/results/index.html`、`catalog.json`、`ref/index.html` 和 `references.bib`，不训练模型。文献说明以 `ref/references.json` 为唯一维护入口。新增任务若使用不同结果格式，需在报告入口加入对应读取规则。

## 文档

- [AGENTS.md](AGENTS.md)：协作、执行与治理规范
- [docs/RESEARCH_PLAN.md](docs/RESEARCH_PLAN.md)：研究领域、假设、路线、评价与阶段门槛
- [docs/DATA_GOVERNANCE.md](docs/DATA_GOVERNANCE.md)：数据来源、权限、质量、去重与版本
- [docs/DECISIONS.md](docs/DECISIONS.md)：已决/待决事项、负责人和证据

用户/导师保留目标物种、研究用途、数据权限、资源和最终验收的决定权。

## 数据工具

在仓库根目录使用统一入口（Python 3.11+）：

```bash
PYTHONPATH=src python3 -m strainfm --help
PYTHONPATH=src python3 -m strainfm download --manifest /path/to/manifest.tsv --root /path/to/genomes
PYTHONPATH=src python3 -m strainfm split-virus --manifest /path/to/manifest.tsv --input /path/to/uhgv.fna.zst --output /path/to/new_genomes
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

`download` 使用清单的 `target_relpath` 和 `gwh_url`（或 `ncbi_url`），检查已有 FASTA；新文件下载、完整解压和序列验证成功后才替换目标，任一下载失败则命令返回非零。它替代旧的 `gwh_retry.py`、`gwh_retry_map.py`；保留的 GWH 发现脚本也使用相同下载实现，运行时须设置 `PYTHONPATH=src`。

`split-virus` 替代 `split_uhgv.py`，需要 `zstd`；使用含 `genome_id`、`species_id` 的清单，目标目录必须尚不存在。仅当整个解压进程成功时发布目录，失败保留原始数据。

`metadata` 与 `metadata-bulk` 统一由 `src/strainfm/metadata.py` 实现，共享 NCBI 字段解析。抽样清单和批量 inventory 在读取时检查必需字段与版本化 accession；批量入口仍校验缓存哈希和续跑契约。模块合并会改变代码哈希，已有批量输出应按原版本续跑，新代码使用新输出目录。

`snapshot`、`metadata-bulk`、`diversity` 拒绝占用没有匹配契约的非空输出目录。恢复必须使用原始输入、配置和匹配代码；修复前的契约不自动迁移，请为新代码选择新输出目录。`diversity` 校验输入 SHA-256，完成后重复调用只核验并读取结果；中断后的距离阶段会重新计算。历史 E. coli 结果仍对应生成时的代码，以上回归测试不构成全量数据重验证。
