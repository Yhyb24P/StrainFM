# Plan A1 运行状态：正式检索基线（进行中）

更新：2026-11-08 18:30 CST。本文件记录方案 A1 的配置、当前进度与待办，供核对是否偏离研究主线（D0-14/D0-20/D0-22/D0-23/D0-24）。结果产物完成后本目录补充指标与验收文件。

## 范围（用户选定，D0-22）

- 查询：validation 2,113（按 support 文件封闭集 1,953 / 开放集 160）
- Gallery：验证查询 452 个 Mash 邻域内的全部训练基因组 22,860
- 合计 24,973；负例池限于该 gallery；测试集封存（审计：A1 清单 test=0）

## 注释（D0-23）

- 来源：**全量统一 pyrodigal**（固定版本、table 11、≥30aa）。GBFF 下载受网络阻塞（仅 8,373/24,973 落盘），已按用户决定停止；GBFF 子集保留，后续作敏感性分析，不混入本轮
- 完成：24,973/24,973，0 失败；划分 22,860 train + 2,113 validation + 0 test
- 蛋白：120,732,992；清单 parquet 38.96GB；SHA-256 `93a8accbdd7b7ac5a167fd024401a6194ba499d31a47bfa811dbd528ffd928bb`
- 覆盖统计：60 个基因组 >6,000 蛋白（推理时截断并计数）；1,549,469 个蛋白 >1,024aa（ESM token 截断；写入 limitations）
- 位置（冷数据）：`/mnt/e/StrainFM/processed/escherichia_coli_bacformer_inputs_planA_20261106/`

## 冻结推理

- 模型：Bacformer-L `84f85c37…`（bf16）+ ESM++ `1e40c1b8…`（bf16）；ESM 单次计算复用为 Bacformer 输入，两向量共享前 6,000 蛋白
- 设备：`device=cuda:0`，`esm_device=cuda:0`（启动时 GPU 空闲 32GB）；契约记录
- 参数：batch 8、max_seq_len 1024、max_n_proteins 6000、max_n_contigs 1000、chunk 16
- 输出：`/mnt/e/StrainFM/processed/escherichia_coli_bacformer_embeddings_planA_20261106/`（分片断点续跑，done.json 含 rows 哈希与设备）
- 进度（18:30）：944/24,973；实测 ~291 基因组/小时；预计 2026-11-12 凌晨完成

## 指标口径修订（D0-24，gpt 对抗性评审后）

- AP@k 分母 = min(k, 正例数)；protocol A 严格按 `--split` 过滤；support 三类偏差逐查询计数；向量读入强校验；ESM auto 回退 CPU；resume 内容哈希校验
- 测试：`tests/test_eval_retrieval_distance.py` 增加 split 隔离测试（合成数据验证 validation/test 互不混入）

## 待办

1. 推理完成后跑协议 A（Bacformer 与 ESM 各一次）
2. 按 gpt 最小验收清单核对：覆盖 24,973、同蛋白集合、NaN/维度/零范数、截断统计、分母 1,953/160、单查询独立复算、tie 样例
3. 结果 README 写入 limitations：Mash 第一层诊断标签、gallery 条件化、伪开放集、统一 pyrodigal 偏差、8,373 GBFF 子集选择偏差、预训练暴露 unknown
4. 提交结果并请 A/gpt 只读评审
