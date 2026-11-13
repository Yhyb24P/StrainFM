"""Build offline literature and experiment views from the saved result tables."""

import csv
import hashlib
import html
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


STYLE = """
.controls label{min-width:0;max-width:100%}.controls select{min-width:0;flex:1}.barrow span{overflow-wrap:anywhere}
:root{--ink:#182c39;--muted:#546979;--paper:#f4f5ef;--teal:#087e82;--blue:#4677b0;--orange:#b45e24}
*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.65 system-ui,'Noto Sans CJK SC',sans-serif}header{background:#173642;color:white;padding:44px max(5vw,20px)}main{max-width:1200px;margin:auto;padding:26px 20px 70px}h1{font-size:clamp(28px,4vw,46px);line-height:1.2}h2{margin-top:0;font-size:26px}h3{margin:10px 0}a{color:#096c87}header a{color:#b8eef0}nav{display:flex;gap:20px;flex-wrap:wrap}section,.card{background:white;border:1px solid #dce3df;border-radius:14px;padding:24px;margin:20px 0}p{max-width:95ch}.muted,small{color:var(--muted)}.eyebrow{letter-spacing:.15em;font-size:13px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:18px}.grid .card{margin:0}.badge{display:inline-block;font-size:12px;padding:3px 9px;background:#e3f2ef;border-radius:20px}.warning{border-left:4px solid var(--orange);padding:12px 18px;background:#fff6eb}.flow{display:flex;gap:10px;align-items:stretch;flex-wrap:wrap;margin:20px 0}.flow span{flex:1;min-width:150px;padding:16px;background:#eaf4f3;border-radius:10px}.flow b{display:block}table{width:100%;border-collapse:collapse;font-size:14px}th,td{text-align:left;padding:10px;border-bottom:1px solid #e2e8e5;vertical-align:top}th{background:#edf3f0}.scroll{overflow:auto}.controls{display:flex;gap:16px;flex-wrap:wrap;margin:20px 0}select,input,button{font:inherit;padding:8px;border:1px solid #a7b9b9;border-radius:6px;max-width:100%}label{display:flex;gap:8px;align-items:center}button{cursor:pointer;background:#e5f2ee}svg{max-width:100%;height:auto}.barrow{display:grid;grid-template-columns:180px 1fr 75px;gap:12px;align-items:center;margin:12px 0}.track{height:22px;background:#edf1ee;border-radius:5px}.bar{height:100%;border-radius:5px;background:var(--teal)}details{margin:14px 0}summary{cursor:pointer;font-weight:600}code{overflow-wrap:anywhere}footer{padding:20px;color:var(--muted)}:focus-visible{outline:3px solid #ca7c32;outline-offset:3px}@media(max-width:650px){section,.card{padding:16px}.barrow{grid-template-columns:115px 1fr 65px;font-size:13px}}@media print{header{background:white;color:black}.controls,button{display:none}section{break-inside:avoid}a{color:inherit}}
"""

ARCHITECTURE = """
<section id="candidate"><div class="eyebrow">重点候选 / CONTEXTUAL REPRESENTATION</div>
<h2>Bacformer-L：给每个蛋白补上基因组上下文</h2>
<p><b>它是什么：</b>把一个细菌基因组视为有序蛋白序列的 Transformer 编码器。它在已有蛋白向量之上学习跨蛋白关系。本项目把它作为冻结的特征提取器，仅训练耐药分类头。</p>
<p>下图是 AMR 路线；operon 任务直接比较相邻蛋白向量，不使用基因组均值和耐药分类头。</p>
<div class="flow" aria-label="Bacformer 数据流"><span><b>① DNA → 蛋白</b>pyrodigal 统一注释<br>按 contig、起点排序</span><span><b>② 单蛋白编码</b>固定 ESM++ small<br>每条蛋白 → 960 维</span><span><b>③ 上下文编码</b>Bacformer-L<br>蛋白序列 → 上下文向量</span><span><b>④ 任务输出</b>基因组均值池化<br>轻量头 → 耐药分数</span></div>
<div class="grid"><div class="card"><h3>30 层 · 15 头 · 960 维</h3><p>自注意力让一个蛋白的表示结合其他蛋白的信息；旋转位置编码表达顺序，contig 信息区分组装片段。前馈层和残差连接继续变换每个位置的表示。</p></div><div class="card"><h3>预训练学什么？</h3><p>遮住部分蛋白信息，利用上下文预测蛋白家族。固定配置的遮蔽比例为 0.15、家族词表为 50,001 类。下游使用编码器隐藏向量，而不是把家族类别直接当作耐药标签。</p></div><div class="card"><h3>本地实际训练什么？</h3><p>ESM++ 和 Bacformer 主干均冻结。每药训练 LayerNorm → Dropout(0.1) → Linear(960,1)，采用 BCEWithLogits 与 Adam。验证集选学习率与 epoch；三个种子使用相同分区。</p></div></div>
<p><b>直观原理：</b>同一条蛋白出现在不同的邻近基因、共现蛋白或基因组背景中，其上下文向量可以不同。这个机制可能帮助预测功能，但当前实验没有分离“顺序”“共现”和其他结构各自的贡献。</p>
<details><summary>展开：一层 Transformer 怎样改变蛋白向量？</summary><p>每个蛋白向量被映射为查询 Q、键 K、值 V。查询与其他位置的键相乘，形成相关性分数；softmax 把分数归一化成权重，再对值向量加权求和：<code>Attention(Q,K,V) = softmax(QKᵀ / √d)V</code>。15 个头从不同投影观察关系，每头维度为 960/15=64；合并后继续经过前馈网络与残差路径。</p><p>旋转位置编码作用在查询与键上，让模型区分顺序关系；contig 编码提供片段身份。30 层重复组合信息，最终仍输出 N×960 的逐蛋白矩阵。对蛋白轴取均值得到一个 960 维基因组向量，再由分类头映射为一个分数。注意力权重本身不能证明生物学因果关系。</p></details>
<p><b>为什么作为主线候选：</b>在扩大队列的同划分、同蛋白输入、同轻量头设置下，Bacformer 对两药的预测均优于匹配 ESM 均值。它仍有可检验的局限：均值可能稀释局部信号，组装顺序可能带来技术影响；这些是解释方向，不是已经证实的缺陷。</p>
<p><b>固定实现：</b>最多取前 6,000 个蛋白；扩大队列中 1 个基因组触发截断。Bacformer revision <code>84f85c37a05e55559142b842c1569e265183b554</code>；ESM++ revision <code>1e40c1b8ef46c33f93a9b817eb0bd81279ab4088</code>。<a href="bacformer_source_audit_20261104/hf_large_config.json">配置证据</a> · <a href="../../ref/index.html#bacformer-code">文献与模型来源</a></p>
<div class="warning">输出是表示或预测分数，不是新氨基酸序列。当前结果不能证明跨来源泛化、预训练未见菌株泛化或物种专属预训练的必要性。</div></section>
<section id="models"><h2>三种模型，各自看到了什么？</h2><div class="grid">
<article class="card"><span class="badge">匹配基线</span><h3>ESM++ 均值</h3><p>单蛋白氨基酸序列 → 蛋白语言模型 → 蛋白均值向量 → 基因组均值 → 轻量分类头。</p><p>它能建模蛋白内部序列，但基因组均值本身不保留蛋白排列；与 Bacformer 比较的是增加上下文模块后的整体收益。</p></article>
<article class="card"><span class="badge">主要研究候选</span><h3>Bacformer-L</h3><p>相同 ESM++ 蛋白向量 → 带顺序和 contig 信息的 Transformer → 上下文向量均值 → 相同分类头。</p><p>原生预训练为跨物种；本项目未继续预训练。它在当前两药受控表示比较中领先。</p></article>
<article class="card"><span class="badge">专项参照</span><h3>StrainAMR</h3><p>DNA → 蛋白簇 / k-mer / 变异图特征 → 多模态模型 → 耐药分数。发布包中的训练特征库用于构建测试特征。</p><p>本地使用环丙沙星 Fold1 发布模型、model-only 模式；关闭最终机制规则融合不等于特征工程没有领域知识。没有在共同训练集上重新训练，不能只按分数归因于架构。</p></article></div>
<p>2026-11-13 修正后的 433 例子集：Bacformer / ESM / StrainAMR 的 AUROC 为 0.885 / 0.614 / 0.930，AP 为 0.836 / 0.591 / 0.948。Bacformer − StrainAMR 的 ΔAUROC 为 −0.0449，95% CI [−0.1074, 0.0142]，不能判定差异，也不能证明等效；ΔAP 为 −0.1125，区间 [−0.1826, −0.0452]，StrainAMR 更高。以上是当前专项参照的表现，不是受控重训的架构排名。<a href="escherichia_coli_ast_amr_external_eval_20261113/external_paired.tsv">精确配对结果</a></p>
<details><summary>专项候选 StrainAMR：为什么输入方式不同？</summary><p>蛋白簇通道把相似蛋白归类，描述基因内容；k-mer 通道利用短 DNA 片段特征；变异图通道描述相对训练参考图的序列变化。多模态模型学习这些特征与药敏标签的联系。它的测试预处理依赖发布包中的训练特征库，不能把全部队列先构建特征再称作独立测试。</p><p>本地运行 build-test 后调用发布的 Fold1 检查点，并使用 model-only 输出。该结果是已训练专项系统的参照，不是从头在本次训练集上训练的同预算对照。原 11 月 12 日冻结分数跨药覆盖的错误版本仍在目录中；当前解释使用 11 月 13 日更正版。</p></details></section>
"""


def run_report(root):
    root = Path(root).resolve()
    results = root / "data/results"
    references = json.loads((root / "ref/references.json").read_text())
    if len({entry["id"] for entry in references}) != len(references):
        raise ValueError("Duplicate reference ID")
    records, pairs, inputs = [], [], {}
    run_dirs = sorted(path for path in results.iterdir() if path.is_dir())
    for directory in run_dirs:
        for path in sorted(directory.rglob("*")):
            if not path.is_file():
                continue
            if path.name not in {"amr_metrics.tsv", "external_metrics.tsv", "metrics_compact.csv", "retrieval_metrics.tsv", "amr_model_comparison.tsv", "external_paired.tsv"}:
                continue
            relative = path.relative_to(results).as_posix()
            inputs[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
            with path.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="," if path.suffix == ".csv" else "\t"))
            if path.name in {"amr_model_comparison.tsv", "external_paired.tsv"}:
                pairs.extend({"run": directory.name, "source": relative, **row} for row in rows)
                continue
            grouped = defaultdict(list)
            for row in rows:
                if path.name == "amr_metrics.tsv":
                    key = (row["drug"], row["analysis_set"], row["model"])
                elif path.name == "external_metrics.tsv":
                    key = (row["drug"], row["analysis_set"], row["model"])
                elif path.name == "metrics_compact.csv":
                    key = ("DH5α operon", "official_input", row["method"])
                else:
                    key = ("Mash 邻域检索", row["split"] + " / " + row["stratum"], row["model_name"])
                grouped[key].append(row)
            for (task, analysis, model), group in grouped.items():
                first = group[0]
                measures = {}
                for metric in ("AUROC", "average_precision", "hit@1", "hit@5", "hit@10", "recall@1", "recall@5", "recall@10", "mrr", "map@10"):
                    values = [float(row[metric]) for row in group if row.get(metric) not in (None, "")]
                    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in values):
                        raise ValueError(f"Invalid {metric} in {path}")
                    if values:
                        measures[metric] = statistics.mean(values)
                if path.name == "amr_metrics.tsv":
                    estimator = "逐种子指标算术平均"
                elif path.name == "external_metrics.tsv":
                    estimator = {"escherichia_coli_ast_amr_external_eval_20261112": "平均预测后计算指标（跨药混用，已失效）",
                                 "escherichia_coli_ast_amr_external_eval_20261113": "逐种子指标算术平均"}.get(directory.name, "外部表点估计（见运行协议）") if model in {"bacformer", "esm_mean"} else "发布单检查点"
                else:
                    estimator = "单次任务评分"
                records.append({"run": directory.name, "task": task, "analysis": analysis, "model": model,
                                "estimator": estimator, "source": relative, "metrics": measures,
                                "n": first.get("n_test", first.get("n_pairs", first.get("n_queries_scored", ""))),
                                "positive": first.get("test_R", first.get("n_positive_pairs", "")),
                                "seeds": [row["seed"] for row in group if "seed" in row], "raw": group})
    for path in sorted(results.glob("*/operon/local_dh5a/eval_*/summary.json")):
        relative = path.relative_to(results).as_posix()
        inputs[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
        summary = json.loads(path.read_text())
        records.append({"run": path.relative_to(results).parts[0], "task": "DH5α operon", "analysis": "local_aligned_pairs",
                        "model": summary["method"], "estimator": "单次任务评分", "source": relative,
                        "metrics": {key: summary[key] for key in ("AUROC", "average_precision")},
                        "n": summary["n_aligned_pairs"], "positive": summary["n_positive_pairs"], "seeds": [], "raw": []})
    runs = []
    for directory in run_dirs:
        count = sum(record["run"] == directory.name for record in records)
        category = "模型结果" if count else "数据 / 协议 / 部署记录"
        if "smoke" in directory.name:
            category = "部署与检索诊断"
        if directory.name == "bacbench_pinned_20261106":
            category = "本地第三方代码快照"
        runs.append({"id": directory.name, "category": category, "record_count": count,
                     "readme": (directory.name + "/README.md") if (directory / "README.md").exists() else None})
    issues = [{"run": "escherichia_coli_ast_amr_external_eval_20261112", "status": "invalid_historical_comparison",
               "reason": "Historical frozen score lookup omitted drug; scores from different drugs overwrote the same model/seed/genome key. Historical frozen metrics and paired comparisons are not valid ranking evidence.",
               "evidence": "Original 20261112 run used _frozen_scores without drug key (subsequently fixed); correctly drug-filtered ciprofloxacin ensemble on 620 samples: Bacformer AUROC 0.916420515574651, AP 0.8257621812680037; ESM AUROC 0.6265887934121017, AP 0.520028771733032."}]
    catalog = {"records": records, "pairs": pairs, "runs": runs, "input_sha256": inputs, "issues": issues,
               "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "note": "从保存的聚合表读取，不重跑模型；不同数据集、统计口径与任务不得合并排名。"}
    results.joinpath("catalog.json").write_text(json.dumps(catalog, ensure_ascii=False, indent=2) + "\n")
    data = json.dumps(catalog, ensure_ascii=False).replace("<", "\\u003c")
    inventory = "".join(f'<tr><td>{html.escape(run["id"])}</td><td>{run["category"]}</td><td>{run["record_count"]}</td><td>' +
                        (f'<a href="{run["readme"]}">运行说明</a>' if run["readme"] else '本地记录 / 无独立说明') + '</td></tr>' for run in runs)
    body = f'''<header><div class="eyebrow">STRAINFM / RESEARCH ATLAS</div><h1>从蛋白上下文，读懂菌株预测</h1><p>历次运行 · 模型结构 · 证据边界</p><nav><a href="#results">结果浏览</a><a href="#candidate">重点候选</a><a href="#models">架构对照</a><a href="#inventory">目录索引</a><a href="../../ref/index.html">参考文献</a></nav></header>
<main><section><h2>当前判断：Bacformer-L 是主要表示候选</h2><p>扩大 E. coli 队列的两药主分析支持其相对匹配 ESM 均值的收益。旧小队列的负结果保留。11 月 13 日已修正外部比较的跨药索引错误：StrainAMR 的环丙沙星 AP 更高，但并非共同训练预算的架构比较。</p><div class="warning">当前测试按近克隆整组隔离，但来源组件跨分区；旧十余例或数百例的结果与扩大队列不能直接作学习曲线。这里没有新增独立验证，也没有评估蛋白生成。</div></section>
<section id="results"><h2>历次运行结果</h2><p class="muted">先选择运行，再选择任务与分析集。条形图从 0 起、统一上限 1；表格保留原始指标与统计单位。历史复算版本可能复用同一批预测，不计作新增实验。</p>
<div class="controls"><label>运行<select id="run"></select></label><label>任务<select id="task"></select></label><label>分析集<select id="analysis"></select></label><label>指标<select id="metric"><option value="AUROC">AUROC</option><option value="average_precision">AP</option><option value="mrr">MRR</option><option value="hit@1">Hit@1</option><option value="recall@10">Recall@10</option></select></label></div>
<p id="context" class="warning"></p><div id="chart" aria-live="polite"></div><div class="scroll"><table><thead><tr><th>模型</th><th>统计口径</th><th>评价 N</th><th>正例</th><th>AUROC</th><th>AP</th><th>来源</th></tr></thead><tbody id="scores"></tbody></table></div>
<details><summary>逐种子与原始聚合行</summary><div id="raw" class="scroll"></div></details><h3>配对差异与置信区间</h3><p class="muted">差值方向为模型 A − B；区间按近克隆组重采样。11 月 12 日旧外部表存在跨药错误且点估计/区间口径不一致；11 月 13 日更正版统一为逐种子指标平均。以对应运行的精确配对表为准。</p><div id="pairs" class="scroll"></div>
<p><a href="catalog.json">下载全部规范化结果与输入哈希 JSON</a> · <button onclick="window.print()">打印当前页面</button></p></section>
{ARCHITECTURE}
<section id="reading"><h2>怎样读这些指标？</h2><p><b>AUROC：</b>检验耐药样本是否排在敏感样本前；<b>AP：</b>汇总精确率—召回率关系，受阳性比例影响，不能当作概率校准指标。单个分数高不等于差异可靠，配对区间用于描述本次固定训练条件下的测试抽样不确定性。</p><p><b>平均分数与集成：</b>“先算每个种子的 AUROC 再平均”与“先平均预测再算 AUROC”通常不同。按正确的环丙沙星行复算，同一 620 例 Bacformer 的两口径约为 0.902 和 0.916。旧 external 表中的 0.927 另含跨药索引错误，不能只用集成口径解释。</p><p><b>历史结论：</b>11 月 8–9 日小队列主分析尚无稳定优势；direct_refit 同时改变训练/验证/测试对象。11 月 12 日扩大队列有更清楚的收益，但数据组成、标签来源和划分也变化，不能单独归因为增加训练样本。</p><p><b>任务单位：</b>operon 评价单位是同一参考菌株的相邻基因对；AMR 是有实验药敏标签的样本；早期检索是 Mash 邻域诊断。三类任务不能用同一排行榜比较。</p></section>
<section id="inventory"><h2>规范化目录索引</h2><p>历史目录保持原路径和哈希；本索引把模型结果、数据/协议和部署诊断按用途归类。未产生指标的运行明确保留为记录，不补造分数。</p><div class="scroll"><table><thead><tr><th>目录</th><th>用途</th><th>聚合行数</th><th>入口</th></tr></thead><tbody>{inventory}</tbody></table></div><p><code>data/raw/</code> 原始输入；<code>data/processed/</code> 派生输入；<code>data/results/</code> 按运行保存结果；<code>ref/</code> 文献与讲解；<code>src/</code> 唯一正式实现。大型权重和序列继续位于本机冷数据目录。</p></section></main><footer>由 python -m strainfm report 生成 · 离线可用 · 所有数值可追溯至保存的结果表</footer>
<script id="data" type="application/json">{data}</script>'''
    script = """
const data=JSON.parse(document.getElementById('data').textContent), $=id=>document.getElementById(id);
const esc=value=>String(value??'').replace(/[&<>\"']/g,char=>({'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;',"'":'&#39;'}[char]));
const labels={bacformer:'Bacformer-L',esm_mean:'ESM++ 均值','strainamr-cipro-frozen':'StrainAMR'};
const displayNames={ciprofloxacin:'环丙沙星',gentamicin:'庆大霉素',primary:'主分析',with_conflicts:'包含冲突组',direct_refit:'直接匹配标签 · 重拟合',primary_all:'外部对照 · 含已知重叠',primary_disjoint:'外部对照 · 排除已知 ID 重叠',official_input:'官方输入',local_aligned_pairs:'本地输入 · 对齐基因对',escherichia_coli_amr_frozen_eval_20261108:'2026-11-08 · 小队列 AMR',escherichia_coli_amr_frozen_eval_20261109:'2026-11-09 · 小队列统计补全与敏感性',escherichia_coli_ast_amr_frozen_eval_20261112:'2026-11-12 · 扩大队列 AMR',escherichia_coli_ast_amr_external_eval_20261112:'2026-11-12 · 外部模型对照（旧版已失效）',escherichia_coli_ast_amr_external_eval_20261113:'2026-11-13 · 外部模型对照（修正）',escherichia_coli_bacformer_frozen_eval_20261106:'2026-11-06 · 单菌株 operon',escherichia_coli_bacformer_smoke_retrieval_20261106:'2026-11-06 · 检索部署诊断',escherichia_coli_bacformer_phase2a_smoke22_20261106:'2026-11-06 · 22 株检索诊断'};
const fmt=value=>value===undefined||value===null||value===''?'—':Number(value).toFixed(4);
function options(id,values){$(id).replaceChildren(...[...new Set(values)].map(value=>{let option=document.createElement('option');option.value=value;option.textContent=displayNames[value]||value;return option}));}
function table(rows){if(!rows.length)return '<p class="muted">此集合没有保存的配对区间；不补算或推定显著性。</p>';const keys=Object.keys(rows[0]);return '<table><thead><tr>'+keys.map(k=>'<th>'+esc(k)+'</th>').join('')+'</tr></thead><tbody>'+rows.map(r=>'<tr>'+keys.map(k=>'<td>'+esc(r[k])+'</td>').join('')+'</tr>').join('')+'</tbody></table>';}
function setTasks(){options('task',data.records.filter(r=>r.run===$('run').value).map(r=>r.task));setAnalyses();}
function setAnalyses(){options('analysis',data.records.filter(r=>r.run===$('run').value&&r.task===$('task').value).map(r=>r.analysis));render();}
function render(){const run=$('run').value,task=$('task').value,analysis=$('analysis').value,metric=$('metric').value;const rows=data.records.filter(r=>r.run===run&&r.task===task&&r.analysis===analysis);
$('context').textContent=task==='Mash 邻域检索'?'部署诊断：N 为实际计分查询数；空指标表示无可计分查询。':task==='DH5α operon'?'单参考菌株：N 为相邻基因对，不代表跨菌株泛化。':data.issues.some(issue=>issue.run===run)?'已失效的历史比较：冻结模型索引遗漏 drug，部分分数被另一药覆盖；这些柱值与配对区间不能支持模型排名。StrainAMR 为发布检查点，且来源不独立、训练预算不匹配。':run.includes('external')?'修正后的发布模型对照：完整药物索引；仍非共同训练预算或来源独立测试。点估计与区间口径请结合原表说明。':run.includes('20261112')?'扩大队列：相同分区和轻量头，近克隆隔离；来源组件仍跨分区。':'历史小队列：primary 在 11 月 8–9 日复用预测；direct_refit 是改变队列后的重拟合。';
$('chart').innerHTML=rows.map(r=>{const value=r.metrics[metric];return '<div class="barrow"><span>'+esc(labels[r.model]||r.model)+'</span><div class="track"><div class="bar" style="width:'+(value===undefined?0:Math.max(0,Math.min(1,value))*100)+'%;background:'+(r.model.includes('strainamr')?'#b45e24':r.model.includes('esm')?'#4677b0':'#087e82')+'"></div></div><b>'+fmt(value)+'</b></div>'}).join('');
$('scores').innerHTML=rows.map(r=>'<tr><td>'+esc(labels[r.model]||r.model)+'</td><td>'+esc(r.estimator)+'</td><td>'+esc(r.n)+'</td><td>'+esc(r.positive)+'</td><td>'+fmt(r.metrics.AUROC)+'</td><td>'+fmt(r.metrics.average_precision)+'</td><td><a href="'+esc(r.source)+'">原表</a></td></tr>').join('');
$('raw').innerHTML=table(rows.flatMap(r=>r.raw));$('pairs').innerHTML=table(data.pairs.filter(r=>r.run===run&&r.drug===task&&r.analysis_set===analysis).map(({run,source,...row})=>row));}
options('run',data.records.map(r=>r.run));$('run').value='escherichia_coli_ast_amr_frozen_eval_20261112';$('run').addEventListener('change',setTasks);$('task').addEventListener('change',setAnalyses);$('analysis').addEventListener('change',render);$('metric').addEventListener('change',render);setTasks();
"""
    results.joinpath("index.html").write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>StrainFM · 模型与结果</title><style>' + STYLE + '</style><body>' + body + '<script>' + script + '</script></body></html>')
    reference_cards, bib = [], []
    downloads_path = root / "ref/downloads.json"
    downloads = {entry["id"]: entry for entry in json.loads(downloads_path.read_text())} if downloads_path.exists() else {}
    for entry in references:
        local = downloads.get(entry["id"], {})
        pdf = f'<a href="papers/{entry["id"]}.pdf">本地全文 PDF</a>' if local.get("status") == "downloaded" else '全文未本地收录；可访问原始来源'
        reference_cards.append(f'<article class="card reference" id="{html.escape(entry["id"])}"><span class="badge">{html.escape(entry["kind"])}</span><h2>{html.escape(entry["title"])}</h2><p class="muted">{html.escape(entry["authors"])} · {html.escape(entry["year"])}</p><p>{html.escape(entry["explanation"])}</p><p><b>与本项目的关系：</b>{html.escape(entry["project"])}</p><p><b>解释边界：</b>{html.escape(entry["limit"])}</p><p><a href="{html.escape(entry["url"])}">原始来源 / 固定版本</a> · {pdf}</p><small>{html.escape(entry.get("version", ""))}</small></article>')
        bib.append('@misc{' + entry['id'] + ',\n  title = {' + entry['title'] + '},\n  author = {' + entry['authors'] + '},\n  howpublished = {' + entry['kind'] + '},\n  url = {' + entry['url'] + '}\n}')
    root.joinpath("ref/references.bib").write_text('\n\n'.join(bib) + '\n')
    root.joinpath("ref/index.html").write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>StrainFM · 参考文献导读</title><style>' + STYLE + '</style><body><header><div class="eyebrow">STRAINFM / READING ROOM</div><h1>参考文献与模型导读</h1><p>从原始来源理解方法，再回到本地证据。</p><nav><a href="../data/results/index.html">模型与结果</a><a href="references.bib">BibTeX</a><a href="references.json">完整书目与讲解</a></nav></header><main><section><h2>阅读路线</h2><p>先读 Bacformer 与匹配 ESM，再读 BacBench 的任务定义；随后看 StrainAMR 的专项输入，最后理解 Mash 分组和数据来源。每条参考均有简要说明，软件和数据库明确标识，不伪装成同行评审论文。</p><label>筛选文献 <input id="search" type="search" placeholder="输入模型、方法或数据源名称"></label><p id="count">共 ' + str(len(references)) + ' 条参考</p></section>' + ''.join(reference_cards) + '</main><script>document.getElementById("search").addEventListener("input",event=>{let count=0;document.querySelectorAll(".reference").forEach(card=>{card.hidden=!card.textContent.toLowerCase().includes(event.target.value.toLowerCase());if(!card.hidden)count++});document.getElementById("count").textContent="显示 "+count+" 条参考"});</script></body></html>')
    return {"runs": len(runs), "metric_records": len(records), "references": len(references), "results_html": str(results / "index.html"), "references_html": str(root / "ref/index.html")}
