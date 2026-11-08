# AMR-01：标签损失与多菌株评价可行性

日期：2026-11-08。负责人：gpt。状态：实际数据核验完成；独立只读复核通过。此为数据/协议可行性结果，不是模型实验。研究优先级来自用户本轮 AMR-01 指示（D0-25）。

## 结论与下一步

旧“609→37”不是单纯的标签稀缺或近克隆去重：**305 条已入选记录被实现中的 accession/genome_id 混键漏计**。原有标签与质量口径下有两药达到 500 总样本、每类 50 的计数门槛，不能据旧报告否决 AMR。

但是，**尚无药物在本交付中被批准进入确认性模型评价**。两药的门槛通过依赖部分 GCA/GCF 配对映射；只用直接 file_accession 匹配时两药均为 496。原来源隔离划分还造成测试样本极少/缺类。下一步先确认配对标签对应的可接受证据与任务隔离要求，再决定是否需要 AMR 专用集合/划分；不直接编码，不放宽阈值凑数。

## 1. 609→37 的逐记录去向

| 互斥类别 | 记录数 | 含义 |
|---|---:|---|
| 旧报告已识别的入选代表 | 37 | genome_id 恰与匹配 file_accession 相同 |
| 已入选但旧键映射漏计 | 305 | genome_id 为 GCF，实际 file_accession 为 GCA；同一清单行关联，无需将标签迁给别的组装 |
| 近克隆非代表 | 247 | 质量合格、有标签原始组装，未入选原检索代表集合；不等于与代表为同一 BioSample |
| 原质量约束排除 | 20 | ambiguous_fraction 超限 10、N50 低于门槛 7、assembly suppressed 3 |
| 合计 | 609 | 37 + 305 + 247 + 20 |

因此，质量合格为 589，实际入选原划分代表为 342，身份去重阶段损失为 0。609 条分别对应 **609 个已知 BioSample**，没有同一 BioSample 的多组装；此处 BioSample 去重后的样本单位不等于流行病学或统计独立。合格 589 条覆盖 **424 个近克隆操作组**，不能按 589 个独立遗传重复解释。没有匹配 accession 版本的降级合并。

旧实现以 file_accession 建立 matched_cohort_rows，再与 genome_id 键的 representatives 求交，导致漏计；修复统一以清单 genome_id 追踪，而标签来源始终保留原 versioned accession。旧 20261106 结果保持原样，D0-18 的否决解释由本报告修正，不将新检查追记为旧实验已完成。

609 个当前压缩 FASTA 均与冻结快照 SHA-256 一致，复用其完整解压/格式检查和已固定质量判据；未重新跑污染度/完整度工具。逐记录保留质量原因、原始身份代表、近克隆代表、继承分区、标签来源与当前文件哈希。

## 2. 逐药结果

下表 N 为质量合格、BioSample 已知、该药非冲突二元标签的样本数；S=0、R=1，沿用发布标签编码。逐药缺失不补负例，不跨药相加。完整 37 药见 [drug_feasibility.tsv](drug_feasibility.tsv)。

| 药物 | 合格 BioSample | S | R | 原划分 train S/R | val S/R | test S/R | 计数门槛 |
|---|---:|---:|---:|---|---|---|---|
| gentamicin | 546 | 440 | 106 | 431/99 | 8/5 | 1/2 | 通过（含配对映射） |
| ciprofloxacin | 545 | 268 | 277 | 263/268 | 5/6 | 0/3 | 通过（含配对映射） |
| ampicillin | 534 | 22 | 512 | 22/504 | 0/6 | 0/2 | 未通过 |
| meropenem | 473 | 437 | 36 | 436/26 | 1/7 | 0/3 | 未通过 |
| amikacin | 470 | 448 | 22 | 442/18 | 4/4 | 2/0 | 未通过 |
| cefepime | 462 | 123 | 339 | 122/330 | 1/7 | 0/2 | 未通过 |

37 药中，2 药通过计数门槛且两类均满足近克隆/来源三分区的必要组数，协议仍待确定；仅 cefoxitin 是未达门槛但通过这些必要条件的探索性候选。其余 34 药在保持当前来源隔离的三分区要求下不可评价：29 药虽有近克隆支持，但某类仅涉及不足 3 个来源分量；另 5 药还缺少近克隆/双类支持（clindamycin、erythromycin 仅 1 个 R；oxacillin、penicillin、telithromycin 无标签）。例如 sulfamethoxazole 的 317 个有效样本全部属于同一来源分量，不能暗中拆开。探索性候选不表示已验证统计功效或已授权降低原协议要求。

标签存在两种命中：554 条候选直接命中 file_accession（也等于 current_accession），55 条只命中同一记录的 paired_accession。两种证据在逐记录表中区分，不将数据库配对声明写成独立验证的序列身份。

| 药物 | 直接 file_accession 命中 S/R（合格 BioSample） | 仅配对映射 S/R | 全口径 S/R |
|---|---|---|---|
| gentamicin | 416 / 80（496） | 24 / 26 | 440 / 106 |
| ciprofloxacin | 259 / 237（496） | 9 / 40 | 268 / 277 |

在可观察的发布二元 CSV 中，未见同一 BioSample/药物的 S/R 冲突、非二元标签或未知 BioSample。该表已折叠原 AST，**原实验内部冲突、I 标签、方法、断点及条件无法由此复核**；“可见冲突 0”不等于原 AST 无冲突。标签来源为作者标注的 NCBI AST，未在本轮重新核查每一条原始 AST。

## 3. 评价协议判定

500/50 的纳入门槛及随机 70/10/20 参数由固定 [BacBench 官方实现](https://github.com/macwiatrak/BacBench/blob/64b7d45823b2ff2ef88458e3aa235caa461a12aa/bacbench/tasks/antibiotic_resistance/train_and_predict_linear.py)核对（filter_drugs 与 train_eval_loop；不是所有 CLI 默认值都等于论文设置）。数据卡另建议五折，不能与随机协议混称。本轮未成功重新获取 v3 PDF（网站 403），不据此声称逐页重核了论文方法。

- **gentamicin**：全匹配口径计数通过；当前 train/val/test 为 530/13/3，测试只有 1 S、2 R。三分区虽都有两类，但远不足以仅凭类别存在就宣称有可靠独立测试；确认性协议未成立。
- **ciprofloxacin**：全口径计数通过；当前 test 为 0 S、3 R，无法在原测试分区评价二类 AUROC，原检索划分不可机械沿用。
- 合格 589 条继承原成员分区后为 train/val/test = 565/16/8。565 条属于同一个旧来源连接分量 PRJDA32509（95.93%）；这是全体身份/近克隆成员与来源相连形成的分量，并不表示所有 565 条自己的 BioProject 都是它。维持该约束不能同时凑成 70/10/20，不能偷偷拆大分量。
- gentamicin 两类分别涉及 336/78 个近克隆组、5/6 个来源分量；ciprofloxacin 为 237/171 个近克隆组、4/6 个来源分量。每类至少 3 组仅是三分区含类的必要条件，不证明存在满足目标比例、来源隔离和统计功效的划分。本轮没有搜索或生成新划分。

允许进入下一步的是**两药标签映射与任务评价设计**，不是直接开展 AMR 模型训练。若无法形成适宜评价，记录限制后检查其他 BacBench 下游任务，不自动扩展 Mash 检索。任务成立后只处理所需有标签菌株，以同样本、同蛋白集合、同标签与同划分比较冻结 ESM/Bacformer。

## 工件与复现

详细证据仅本地保存，不自动提交临床相关逐样本元数据：

- `data/processed/escherichia_coli_amr_feasibility_source_isolated_20261108/record_disposition.tsv`：609 行；逐组装去向、旧键漏计、来源/分组、排除原因和文件校验。
- 同目录 `sample_drug_evidence.tsv`：609 × 37 个 BioSample/药物单元，含缺失/冲突状态及真实标签所在组装 ID；没有训练集选择结果。
- 本目录聚合统计、contract、summary 是该次运行产物的同字节副本；原始标签 CSV、FNA、逐记录证据不发布。`verification.json` 记录独立计数复核、产物哈希和测试。

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m strainfm amr-intersection \
  --labels data/results/bacformer_source_audit_20261104/hf_amr_binary_labels.csv \
  --cohort-manifest data/processed/escherichia_coli_identity_cohort_20261105/cohort_manifest.tsv \
  --split-manifest data/processed/escherichia_coli_source_split_20261105/split_manifest.tsv \
  --groups data/processed/escherichia_coli_mash_groups_20261105/sequence_groups.tsv \
  --genome-root /mnt/e/StrainFM/data/raw/genomes \
  --output <新的非raw输出目录>
```

源码基点提交：`e541e746407de2909ca1971d0b3eae0b95ce10b1`；运行包含本轮未提交 AMR 修复，精确源码与输入 SHA-256 见 [contract.json](contract.json)，不能仅用基点提交代表新代码。

标签版本：HF revision `fa5f826a5a6ae28126483a751c6e8acffaca0eec`，CSV SHA-256 `644b91687d53c36e8645f0a5eb560b33cd807276ae41d81ae9a2f411050b0f63`；来源 URL/下载时间见 `data/results/bacformer_source_audit_20261104/sources.json`。快照、cohort、split 与 Mash 分组哈希固定于 contract。

全量 Plan A 编码 PID 1196785 在本轮用户指示后以 SIGINT 中断，随后确认进程已不存在；85 个 rows 分片和 85 个 done 标记保留，原进度 1360/24973，没有 completion。未删除原始数据、历史结果、权重或已落盘模型产物。操作只停止继续编码，不能撤销之前已耗费的计算。

本轮未重新运行 operon，未读取 AMR 模型测试结果，未训练下游头或主干，未开始继续预训练、额外下载或蛋白生成任务。编码、标签可行性与模型能力结果严格分开。
