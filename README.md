# StrainFM

当前主线：评估 Bacformer 在大肠杆菌真实标签下游任务上的能力，再依据证据决定模型研究切入点。

## 项目定位与当前状态

长期目标是单物种菌株级基础模型。近期比较冻结 Bacformer 与匹配 ESM 表示在可靠标签上的预测收益；表示、功能预测、检索与蛋白序列生成分别取证。

DH5α operon 官方输入和本地输入的评价闭环已完成，保留为单参考菌株结果，不重复开展，也不外推多菌株泛化。数据快照、质量检查、身份与 Mash 分组继续作为支撑；停止扩展无标签检索和全量编码，保留已落盘产物。

当前交付 **AMR-01**：解释旧报告 609→37 的损失，统计逐药 BioSample、S/R、冲突和质量情况，再判断原协议门槛与独立评价是否成立。详见[AMR-01 结果](data/results/escherichia_coli_amr_feasibility_20261108/README.md)。本次不新选训练集、不新划分、不运行模型。

只有在任务成立后，才针对所需菌株提取匹配表示并比较预测；物种适配、数据规模与生成评价均依据后续能力证据决定。研究路线统一见[研究计划](docs/RESEARCH_PLAN.md)，已有 operon 证据见[运行记录](data/results/escherichia_coli_bacformer_frozen_eval_20261106/README.md)。

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

`snapshot`、`metadata-bulk`、`diversity` 拒绝占用没有匹配契约的非空输出目录。恢复必须使用原始输入、配置和匹配代码；修复前的契约不自动迁移，请为新代码选择新输出目录。`diversity` 校验输入 SHA-256，完成后重复调用只核验并读取结果；中断后的距离阶段会重新计算。历史 E. coli 结果仍对应生成时的代码，以上回归测试不构成全量数据重验证。
