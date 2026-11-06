# StrainFM

当前主线：实际数据审计与划分 → 冻结表示基线 → 条件性的物种适配。

## 项目定位

本项目研究**细菌种内精细遗传多样性的基因组基础模型与评价方法**。技术起点是 Bacformer v3 的全基因组有序蛋白质上下文表示；核心目标不是复刻一个更小的 Bacformer，而是检验并改善跨物种模型对菌株级差异（核心基因等位变异、accessory genome、基因排列、移动遗传元件）的表征能力。

暂定英文方向：**Species-adaptive genome representation learning for fine-grained bacterial strain diversity**。当前三个阶段只在**单一目标物种**内验证。

## 核心研究问题

1. **表征缺口**：跨物种 Bacformer 在种内核心变异、可变基因和基因组组织上，保留了什么、遗漏了什么？
2. **物种适配**：目标物种继续预训练或专用监督，能否在未见菌株与独立来源上超过冻结模型和传统基线？
3. **数据效率**：新增基因组数量与新增种内多样性，哪个更能解释性能增益？
4. **评价有效性**：如何用独立标签、谱系感知划分和公平基线，避免循环论证与数据泄漏？

第一轮以固定的遗传结构与检索指标检验这些问题，不能只用 UMAP、聚类图或预训练 loss 证明价值；表型任务按后续研究用途增加。

## 三个连续阶段

1. **实际数据审计**：固定下载版本，确认目标物种、真实文件状态、元数据覆盖、独立样本数量和种内多样性；先划分开发集与独立测试集，再做模型预处理。
2. **可信表示基线**：同批基因组对比冻结 Bacformer、匹配的 ESM 平均蛋白表示和传统基因组方法，检验种内结构恢复与未见菌株检索。
3. **物种适配价值**：基线支持后比较物种专属继续预训练与冻结模型，控制训练数据和计算预算，以数据规模实验判断泛化收益。

当前处于第一阶段：E. coli 候选数据的全量快照、元数据、身份/近克隆合并、多样性分组和来源隔离划分已完成并核验；开发集含训练 26,462、验证 2,113 个代表，独立测试候选 2,112 个。正式物种选择待用户确认，尚无模型预处理或实验结果。详见[候选数据交付](data/results/escherichia_coli_data_preparation_20261105/README.md)。审计仅保留会影响纳入、划分和结果解释的检查；完成条件见[研究计划](docs/RESEARCH_PLAN.md)。

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
