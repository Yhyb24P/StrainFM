# 参考文献与导读

打开 [index.html](index.html) 阅读全部收录参考的简要讲解、项目用途与解释边界。书目维护入口是 [references.json](references.json)，[references.bib](references.bib) 由统一报告命令生成。

`papers/` 收录实际下载成功的开放获取全文；来源 URL、文件大小和 SHA-256 见 [downloads.json](downloads.json)。Bacformer v3 原站当前返回 403，未伪装成已获得全文；其 DOI、固定模型配置与官方实现入口均保留。StrainAMR 依作者 README 按软件引用记录，未编造 DOI，未复制其许可未定的源码和权重。

论文、官方模型说明、软件实现和数据资源在 HTML 中明确区分。背景候选（如 gLM2、Panaroo）没有被列为本项目已运行模型。新增参考必须同时填写简要说明、与项目关系和适用边界；重建命令：`PYTHONPATH=src python3 -m strainfm report`。
