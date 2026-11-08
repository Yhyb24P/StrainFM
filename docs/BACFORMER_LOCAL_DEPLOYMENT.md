# Bacformer 本地部署与评测环境规范

版本：2026-11-06。本条固定第二阶段（冻结表示）在本机的部署事实与运行约定。研究路线见 [RESEARCH_PLAN.md](RESEARCH_PLAN.md)，决策见 [DECISIONS.md](DECISIONS.md)（D0-17/D0-18）。

## 1. 范围与硬约束

- 只做**冻结推理**：不继续预训练、不微调主干；AMR 类任务仅允许论文式轻量下游头，且须在闭环 A 通过后另行决定（D0-17）。
- 主评测框架为 **Bacformer 官方模型 + BacBench 官方评测框架**；StrainFM 只保留薄适配层（固定划分、固定模型快照、结果落盘），不重写模型或另建平行评测系统。
- GPU 为共享资源：每次 GPU 推理前先 `nvidia-smi` 确认余量；不足时回退 CPU 并在结果中如实记录 `device`。
- 冷数据（权重、数据集、大型中间产物）放 `/mnt/e/StrainFM/models`，不入 Git；Plan A 的注释清单与冻结向量放 `/mnt/e/StrainFM/processed/`，仓库 `data/results/` 只保留小型结果与运行记录。NCBI FTP/HF 下载走本机代理 `127.0.0.1:28080`（HF 可选 hf-mirror，镜像缺 revision 时回退代理）。

## 2. 硬件与环境（本机实测）

| 项目 | 实测值 |
|---|---|
| GPU | NVIDIA GeForce RTX 5090 D，32,607 MiB，Capability (12,0)（Blackwell sm_120） |
| 驱动 | 610.88；CUDA 运行时由 PyTorch wheel 自带（cu13x 构建） |
| 系统 | WSL2，Python 工具链沿用现状 |
| Conda 环境 | **复用 `esm3`**（用户指示优先复用本项目环境，不新建） |
| torch | 2.11.0（CUDA 可用，`torch.cuda.is_available() == True`） |
| transformers | 4.48.1（固定，未升级） |
| bacformer | 0.2.0（pip 补装） |
| 其他 | pyarrow 25.0.1、pytorch-lightning 2.6.6、datasets 5.0.1、tokenizers 0.21.4 |
| flash-attn | **未安装，不作为前置条件**（官方列为可选加速） |

安装约束（复现用）：

```bash
/home/yhshy/miniconda3/envs/esm3/bin/python -m pip install \
  'bacformer==0.2.0' 'transformers==4.48.1' 'torch==2.11.0' \
  pyarrow pytorch_lightning fsspec datasets
```

注意：`bacformer` 的元数据允许 `transformers>=4.38.2`，不 pin 会把共享环境升到 transformers 5.x；本环境固定 4.48.1。原始方案中的 `pip install torch --index-url .../cu121` **不适用于 5090 D**（sm_120 需更新的 CUDA 构建），已按机器实际情况替换。

### 2.1 环境验收标准

1. 模型可从固定本地快照加载；
2. 完整基因组前向计算成功，输出无 NaN/Inf；
3. 重复运行一致：**同设备复跑 embedding 的 max abs diff ≤ 1e-6**（2026-11-06 CPU 实测为 0.0，逐位一致）；CPU 与 GPU 结果不做比特级比较，跨设备只要求指标差异可解释；
4. 记录峰值显存（GPU）或峰值内存（CPU）、软件版本与 `device`。

能 `import bacformer` 不等于部署完成。

### 2.2 已完成的验收记录（2026-11-06，deployment verification）

- 闭环 A 官方输入（DH5α）在 CPU 上跑通：Bacformer-L AUROC 0.7125 / AP 0.2194；ESM++ 均值对照 AUROC 0.6724 / AP 0.1777（4,439 对，574 正例对；单菌株，仅证明部署与官方评测代码可复算）。明细见 `data/results/escherichia_coli_bacformer_frozen_eval_20261106/`。
- 资源：Bacformer 运行 37.13s / 峰值 RSS 9.83GB；ESM 运行 73.21s / 峰值 RSS 5.27GB（CPU、OMP 16）。
- 确定性：同设备复跑逐位一致（max abs diff = 0.0）。
- 本地 DH5α 输入（同版本 GBFF，4,444 蛋白，显式映射 4,429/4,440）在 4,408 个对齐邻接对上复算：Bacformer-L 0.71262/0.22060，官方同对集合 0.71186/0.22041；ESM++ 0.67177/0.17868 vs 0.67176/0.17868。注释漂移无可观测影响（详见结果目录 `operon/local_dh5a/`）。

## 3. 版本锁定（正式评测必须使用）

| 组件 | 固定标识 | 本地快照 |
|---|---|---|
| Bacformer-L（complete genomes） | HF revision `84f85c37a05e55559142b842c1569e265183b554`（commit message `fix-rope-buffer-reinit`，2025-05-15 embedding 质量 bug 修复后） | `/mnt/e/StrainFM/models/hf/models--macwiatrak--bacformer-large-masked-complete-genomes/snapshots/84f85c37a05e55559142b842c1569e265183b554` |
| ESM++ small（匹配 pLM 对照） | HF revision `1e40c1b8ef46c33f93a9b817eb0bd81279ab4088`（BacBench `ESMPlusPlusEmbedder._load` 内嵌 pin，为 Bacformer-L 训练时版本） | `/mnt/e/StrainFM/models/hf/models--Synthyra--ESMplusplus_small/snapshots/1e40c1b8ef46c33f93a9b817eb0bd81279ab4088` |
| BacBench | git commit `64b7d45823b2ff2ef88458e3aa235caa461a12aa` | `data/results/bacbench_pinned_20261106/BacBench`（本地，不入 Git） |
| Bacformer 源码（审计参照） | commit `61dca9537429b95365d73d72a3e11c207df1ed6c` | `data/results/bacformer_source_audit_20261104/Bacformer` |
| operon 长读长数据集 | dataset revision `d77df90f19544475a4da28c0cb80cda2e4d00652` | `/mnt/e/StrainFM/models/hf_datasets/datasets--macwiatrak--operon-identification-long-read-rna-sequencing-protein-sequences/snapshots/d77df90f19544475a4da28c0cb80cda2e4d00652` |

- 权重 SHA-256、全部 revision 与路径记录在 `/mnt/e/StrainFM/models/revisions.json`；正式结果目录须复制该清单或引用其哈希。
- `trust_remote_code=True` 加载的模型实现属于实验的一部分（随快照固定），不得隐式追踪远程 `main`。
- 运行固定快照时使用离线解析：

```bash
export HF_HOME=/mnt/e/StrainFM/models/hf
export HF_HUB_CACHE=/mnt/e/StrainFM/models/hf   # 与下载时 cache_dir 一致，缺一不可
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
```

## 4. 已知上游行为与坑（固定 commit 上核查）

| 位置 | 行为 | 处理 |
|---|---|---|
| `run_embed_prot_seqs.py` | `model_path` 含 `bacformer-large` 时**硬编码**把 pLM 换成 `Synthyra/ESMplusplus_small`（revision 由 embedder 固定），不暴露独立 revision 参数 | 首次功能验证按官方运行；正式 runner 必须支持分别传入 Bacformer 与 pLM 的本地快照 |
| `run_embed_prot_seqs.py` | `load_seq_embedder(model_path)` **不传 device**，pLM 默认落在 CUDA；`--device cpu` 只影响 Bacformer | 本仓 runner 用 `--device`（Bacformer）与 `--esm-device auto`（ESM）分别控制，两者都写入 contract；CPU-only 机器上 ESM 会显著变慢 |
| 共享 GPU | 运行中外部进程可抢占显存/算力（实测 ESM 阶段因外部作业 97% 占用而变慢） | 运行记录如实报告墙钟时间与 `device`；不把共享排队时间当作模型性能 |
| `run_embed_prot_seqs.py` | `--max-n-proteins` 默认 9000，模型训练上限为 6000 | 显式传 6000；超长基因组单独统计截断 |
| `run_embed_prot_seqs.py` | 输出时 `df["split"] = split_name` | 划分永远以独立 split manifest 为准，按 `genome_id` 重连并断言 |
| `embedder.py` | 按路径字符串匹配模型类型（`esmplusplus`/`bacformer`） | 本地快照路径保留可识别子串；换名需同步检查 |
| `preprocess.py`（Bacformer 源码） | GenBank 分支跳过无 `translation` 的 CDS、不保存 strand；GFF 分支不完整；`strain_name` 取 `organism` | 不直接作为正式输入；适配层保留 strand/protein_id/locus_tag/坐标与跳过计数 |
| 分组划分 | 失败后静默回退样本随机划分 | 正式分组实验失败即停止 |
| 二进制标签 | 字符串标签按排序编码 | 显式固定 S=0、R=1 并检查阳性方向 |
| `test_after_train` | 默认关闭 | 不把 `val_*` 当最终测试结果 |

## 5. 两条闭环与输入契约

### 5.1 闭环 A：DH5α operon（官方标签）

- 参考：`GCF_002899475.1`（DH5α），contig `NZ_CP026085.1`，4,440 个蛋白、369 个 operon；本项目清单与固定划分中均包含该 accession。
- **第一先跑作者原始蛋白列表与 `operon_prot_indices`**：`data/processed/escherichia_coli_operon_reference_20261106/ecoli_operon_reference.parquet`（SHA-256 `b54cdb7c1dbd003ab51a01175d42af810156bec128a5f24eabbb9d055134c487`）。
- 官方流程（不聚合全基因组，需逐蛋白向量）：

```bash
cd data/results/bacbench_pinned_20261106/BacBench
PYTHONPATH=$PWD HF_HOME=/mnt/e/StrainFM/models/hf HF_HUB_CACHE=/mnt/e/StrainFM/models/hf \
HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 \
python bacbench/modeling/run_embed_prot_seqs.py \
  --input-parquet-path data/processed/escherichia_coli_operon_reference_20261106/ecoli_operon_reference.parquet \
  --model-path <bacformer 本地快照> --batch-size 8 --max-prot-seq-len 1024 --max-n-proteins 6000 \
  --device <cuda:0|cpu> --output-filepath <结果>/operon/official_input/bacformer/ecoli_operon_bacformer.parquet
python bacbench/tasks/operon/run_evaluation_long_read_rna_seq.py \
  --input-filepath <上一步输出> --output-filepath <结果>.csv
```

- ESM 对照用同一命令、`--model-path <ESM++ 本地快照>`（无 Bacformer 阶段），同一标签与链规则。
- 主结果用官方 benchmark 实现：`(cos+1)/2`、异链置零；代码中的 `AUPRC` 实为 `average_precision_score`，报告时写明。论文带间距的分数如需复现，单独命名报告。
- 结果表述限定为“DH5α 实验标注参考上的 operon 边界评价”，不外推为多菌株泛化。

### 5.2 闭环 B：多菌株 AMR（先完成 AMR-01）

2026-11-06 的交集报告在 accession 与 genome_id 的键映射上漏计，不能据“37 个代表”认定 AMR 不可行。2026-11-08 按用户要求逐记录追踪质量、代表选择和键映射，并统计 BioSample 与逐药类别/分区支持；现役结论见 [AMR-01 报告](../data/results/escherichia_coli_amr_feasibility_20261108/README.md)。旧结果保留为历史，尚未进行多菌株 AMR 模型训练或测试。

下游评价是否成立，先看真实标签、质量合格的原始组装、样本与近克隆/来源隔离及逐药支持。不自动扩大检索，也不因标签总数足够就直接启动全量推理。不同任务可有专用集合，须在查看模型测试表现前固定；不将标签复制给近缘组装。

## 6. 数据链约定（E. coli 自有基因组）

- 优先使用同版本配套注释（GBFF 或 FAA+GFF）；仅在无配套注释时用固定版本 pyrodigal（translation table 11）。
- GBFF 查找同时尝试 `genome_id`/`file_accession`/`paired_accession`（GCA/GCF 配对命名）；0 CDS 的 WGS 主记录 GBFF 自动回退 pyrodigal 并在 `summary.json` 记 `n_gbff_fallbacks`。
- 源文件在慢速/共享挂载（如 WSL2 `/mnt/e` 的 9p）上解析前先复制到本地临时目录：实测并发读同一 3MB GBFF 时出现重复读放大（累计读取 >70GB、分钟级停滞），复制后解析约 0.7s。
- 注释产物记录 `n_cds_without_translation`、`n_pseudo`、`n_proteins_over_1024aa`、`n_over_max_proteins`，并在 `genome_input_summary.tsv` 写入 `source_file` 以回溯。
- 有实验标签的参考菌株**不得重注释**后按行号贴回标签；须建立显式映射并列出失败对象。
- 输入排序固定为：contig 在 FASTA 中的顺序 → CDS start → end；`max_n_proteins=6000` 截断规则固定并记录 `truncated_protein_count`。
- 蛋白质级字段必须可回溯到“哪个组装、哪个 contig、哪个 CDS”：保留 genome_id、assembly accession（含版本）、BioSample、contig、contig 内索引、坐标、strand、locus tag、protein_id、序列哈希。

## 7. 测试集纪律

1. 官方 DH5α 输入先跑通（本文件 §2.1 验收 + operon 指标可复算）；
2. 再经显式映射跑本地 DH5α 输入并解释差异；
3. 验证流程稳定后冻结配置（写入 DECISIONS/结果 README）；
4. 测试集一次性开封；开封后不改主指标、k、正例或链规则。

## 8. 运行记录清单

正式结果必须记录：`model_revision`（两级）、`bacbench_commit`、`bacformer_version/transformers_version/torch_version/cuda_version`、`gpu_name` 或 `cpu`、完整命令、`input_manifest_sha256`/`split_manifest_sha256`、`code_commit`、`device/dtype`、峰值内存与运行时长。
