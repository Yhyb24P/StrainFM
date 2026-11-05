# Bacformer 与 BacBench 固定版本证据

负责人：Codex 当前会话。日期：2026-11-04。状态：公开源码和元数据已取回并做静态审查；未安装或执行第三方模型、未下载权重。

`sources.json` 记录每个文件的确切 URL、抓取时间、revision（若适用）、字节数和 SHA-256。`Bacformer/`、`BacBench/` 是用于审计的上游文件快照（未随仓库发布；URL、revision 和 SHA-256 见 `sources.json`），附各自 LICENSE，不是项目的第二套正式实现。`local_context.json` 记录读取的项目文档哈希、当前 PATH 工具可用性和资源快照。

版本：

- Bacformer `61dca9537429b95365d73d72a3e11c207df1ed6c`；BacBench `64b7d45823b2ff2ef88458e3aa235caa461a12aa`。
- Large complete-genomes `84f85c37a05e55559142b842c1569e265183b554`；Base `d1c461f1db75902e6502e5ba7588ea1a4ea96b29`。参数量取 HF 元数据，不是本机加载结果。
- 药敏标签 revision `fa5f826a5a6ae28126483a751c6e8acffaca0eec`，只下载标签 CSV 和卡片，未下载大型序列 parquet。

复核静态结论时查看：

- `BacBench/bacbench/modeling/utils/utils.py` L128–211：9000 默认、前缀截断、Base 的 SEP 默认引用 CLS。
- `BacBench/bacbench/modeling/embed_prot_seqs.py` L128–207：6000 默认、调用参数和 pooling。
- `BacBench/bacbench/tasks/antibiotic_resistance/train_and_predict_linear.py` L131–146：group split 异常后的随机回退。
- `Bacformer/bacformer/modeling/modeling_large.py`：contig embedding、旋转位置编码、masked 蛋白簇分类头。
- `hf_large_config.json`：6000、960 维、50001 类等 checkpoint 配置。
- `biorxiv_versions.json`：三版日期及 published=NA；不能据此核验论文正文的结果数字。

复核某个公开文件：按 `sources.json` 中的 URL 下载到新临时文件，然后计算 SHA-256，与记录比较；不要直接覆盖快照。原始网页正文有更新时，以固定 commit URL 或记录的抓取文件为准。

访问限制与未完成项：bioRxiv v3 正文返回 403，未确认 85 页/九群/ANI-Panaroo 细节；部分出版社/PMC 页面重定向或触发访问检查，已使用能访问的官方来源支撑对应论点；kimi-webbridge 健康检查未成功，网页核验使用 web 工具与 urllib 官方 API。没有通过这些失败推断论文内容，也没有联系作者或协作者。

阅读入口：[研究报告](../../../docs/STRAIN_REPRESENTATION_RESEARCH_20261104.md)。
