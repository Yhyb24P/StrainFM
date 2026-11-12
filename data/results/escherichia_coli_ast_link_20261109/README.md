# E. coli 实验药敏标签与基因组候选链接

日期：2026-11-09。负责人：gpt；CLI 接入：ds。状态：候选链接已执行、定向测试通过；基因组补下载进行中，尚未形成模型就绪队列。代码基线 `0d93927` 加本目录 `artifacts_sha256.json` 所列工作树源码。

本工作落实用户附件中的“先补菌株—基因组—真实标签对应表，再扩大共同评测”。资源调研与实验方案见[来源调研](../escherichia_coli_label_sources_survey_20261109/README.md)。不重新开展旧 AMR-02 或 DH5α 评价。

## 实际输入与结果

输入为固定 BV-BRC 查询快照：公开、taxon_id 562、`Laboratory Method`、环丙沙星与庆大霉素，共 30,961 条记录；ENA 别名映射；已有本地质量队列；ATB 2025-05 聚合组装元数据；ENA run/study 查询。文件级哈希见 `summary.json`，下载日期、元数据 URL、字节数及哈希见 `atb_metadata_source.json`。ATB 压缩元数据为 68,738,428 字节，已完整解压验证。源文件保留在本地，逐样本表不提交。

3,945 个规范 BioSample 至少具有一药无冲突 S/R 记录。它们不是已经确认相互独立的菌株；近克隆和来源分组尚未完成。

| 动作 | BioSample 数 |
|---|---:|
| 复用本地合格文件且实际 SHA-256 一致 | 594 |
| ATB 元数据通过且无已知来源物种冲突，进入下载清单 | 2,950 |
| ATB 质量/物种等元数据不通过 | 390 |
| ATB 元数据未找到 | 8 |
| accession 与 BioSample 身份冲突隔离 | 2 |
| ENA 与 ATB 物种声明冲突隔离 | 1 |

ATB 纳入要求 `asm_pipe_filter`、`sylph_filter`、`hq_filter` 全部 PASS，物种为 E. coli，有组装、测序 run 和预期单样本 URL。此要求不能代替下载后实际 FASTA 质量检查。ENA 中另外两条 K. pneumoniae 声明已经被 ATB 物种检查排除；单独隔离的那条是 ENA 声明 K. pneumoniae、ATB 声明 E. coli，保留双方声明，不裁定哪个正确。

下表是 3,544 个可复用或待下载候选的标签覆盖，不是最终建模计数；S/R 外状态不参与该药的二元训练。

| 药物 | S | R | I | unknown | conflict |
|---|---:|---:|---:|---:|---:|
| ciprofloxacin | 1,950 | 1,356 | 47 | 164 | 27 |
| gentamicin | 2,536 | 813 | 56 | 138 | 1 |

缺失记录不投票，不解释为敏感；同一 BioSample 出现多个不同非缺失调用（包括 I）标为冲突。不同 BioSample 之间不转移标签，accession 主体相同不替代准确版本对应。此处按 BioSample 汇总，还不能证明实验材料与 ATB 同样本重组装在所有时间点完全等同。

ENA 查询覆盖过滤前 3,943 个身份无冲突候选：3,754 个返回 run，共 3,850 条、42 个研究；14 个样本有多个研究编号，全部保留。其余 189 个未返回 run；没有将缺失来源填成某个现有研究。完整聚合见 `ena_source_summary.json`；这些数量先于 ATB 筛选。

## 重现

在仓库根目录使用统一入口，输出目录须不存在：

```bash
PYTHONPATH=src python3 -m strainfm amr-sources \
  --source-dir data/processed/ecoli_external_sources_20261109 \
  --cohort-manifest data/processed/escherichia_coli_identity_cohort_20261105/cohort_manifest.tsv \
  --local-genome-root /mnt/e/StrainFM/data/raw/genomes \
  --atb-assemblies /mnt/e/StrainFM/data/raw/allthebacteria_metadata_202505/atb.metadata.202505.sqlite.assembly.tsv.xz \
  --ena-runs /mnt/e/StrainFM/data/raw/ecoli_ast_ena_sources_20261109 \
  --output data/processed/escherichia_coli_ast_source_checked_20261109

PYTHONPATH=src python3 -m strainfm download \
  --manifest data/processed/escherichia_coli_ast_source_checked_20261109/download_manifest.tsv \
  --root /mnt/e/StrainFM/data/raw/allthebacteria_ecoli_ast_20261109
```

下载使用现有实现；完整解压和 FASTA 验证通过才发布目标文件，失败命令非零退出，重试时复用已验证文件。第二条命令已启动，进行中观测见 `download_progress.json`；本报告不把下载清单计数写作实际文件数。

验证：`test_amr_sources.py` 6 项和 `test_download_virus.py` 6 项通过，`git diff --check` 通过。测试涵盖缺失/冲突标签、样本错配隔离、快照篡改与路径拒绝、别名解析、ATB 物种/QC、ENA 多研究与物种冲突，以及原下载入口兼容。实际 CLI 已生成上述计数。不是新模型效果的验证。

## 协作与未完成事项

AIA scope 为现有 `.state/aia`；沿用当前 gpt 会话。ds 消息 `UPD-ds-20261109-link-verified-048` 核验初始链接的 4 个哈希与计数；`UPD-ds-20261109-ena-cli-049` 完成 `--ena-runs` CLI 接线并核对 ATB 筛选计数。最终加入 ENA 隔离的输出与源码哈希列于 `artifacts_sha256.json`，不将先前复核泛化为对所有后续改动的签核。ds 持有 `__main__.py` 写入范围，gpt 负责链接实现和测试。

尚需下载后的文件哈希与质量纳入、近克隆/来源分组、训练验证测试共同划分，以及扩大后的三方法评价。旧测试集已看过的结果不能充作新的确认性证据；未知预训练暴露仍为 unknown。StrainAMR 固定源码已落盘，独立环境补充 SHAP 后 `build-train --help` 通过；外部工具和数据库尚未就绪，环境版本差异及一次依赖求解失败见 `strainamr_runtime_preparation.json`。入口检查不等于依赖完整或训练成功；当前尚未运行 StrainAMR 特征构建/训练。此阶段没有新的模型优劣结论。
