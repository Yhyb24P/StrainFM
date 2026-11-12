"""Score an external pretrained AMR model against frozen ESM/Bacformer heads (D0-31).

The external model's own training data may overlap the frozen test partition, so
two evaluation sets are pre-registered: `primary_all` (label-conflict-free test
units) and `primary_disjoint` (primary minus genomes listed in --overlap-ids, i.e.
those present in the external model's training set). Scores from the frozen heads
come from the amr-eval `primary` analysis; subsetting the test rows does not
refit any head.
"""

import json
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

from .amr_eval import _group_bootstrap_multi, _metrics, _paired_group_bootstrap_multi
from .audit import sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output

METRIC_COLUMNS = ["drug", "analysis_set", "model", "n_test", "test_S", "test_R",
                  "AUROC", "average_precision", "AUROC_ci95", "AP_ci95", "n_groups", "n_valid_repeats"]
PAIR_COLUMNS = ["drug", "analysis_set", "model_a", "model_b", "n_test_units",
                "delta_AUROC_mean", "delta_AUROC_sd", "delta_AUROC_ci95",
                "delta_AP_mean", "delta_AP_sd", "delta_AP_ci95"]
PREDICTION_COLUMNS = ["drug", "analysis_set", "genome_id", "near_clone_group", "label", "score"]


def _frozen_scores(eval_dir):
    """Per-model, per-seed test scores from the amr-eval primary analysis."""
    rows = [r for r in read_table(Path(eval_dir) / "amr_predictions.tsv")
            if r["analysis_set"] == "primary" and r["split"] == "test"]
    by_model = {}
    for row in rows:
        by_model.setdefault(row["model"], {}).setdefault(int(row["seed"]), {})[row["genome_id"]] = float(row["score"])
    return by_model


def _score_sets(units, overlap_ids):
    sets = {"primary_all": lambda unit: unit["label_conflict"] != "True"}
    if overlap_ids is not None:
        sets["primary_disjoint"] = lambda unit: (unit["label_conflict"] != "True"
                                                 and unit["genome_id"] not in overlap_ids)
    return sets


def _read_external_scores(path):
    """Accept `genome_id,score` or StrainAMR's native `Sample_ID,Label,Pred,Prob` table."""
    rows = read_table(path)
    if not rows:
        raise ValueError(f"Empty external scores file: {path}")
    if "Sample_ID" in rows[0]:
        id_field, score_field = "Sample_ID", "Prob"
    elif "genome_id" in rows[0] and "score" in rows[0]:
        id_field, score_field = "genome_id", "score"
    else:
        raise ValueError("External scores need genome_id+score or Sample_ID+Prob columns")
    scores = {}
    for row in rows:
        genome_id = row[id_field]
        if genome_id in scores:
            raise ValueError(f"Duplicate external score for {genome_id}")
        scores[genome_id] = float(row[score_field])
    return scores


def run_amr_external_eval(eval_dir, units_path, external_scores_path, output, external_name,
                          overlap_ids_path=None, bootstrap_repeats=1000, seed=20261110):
    eval_dir, units_path, external_scores_path, output = map(
        Path, (eval_dir, units_path, external_scores_path, output))
    if "raw" in output.resolve().parts:
        raise ValueError("Non-raw output directory required")
    if output.exists():
        raise ValueError("Use a new output directory")
    units = read_table(units_path)
    external = _read_external_scores(external_scores_path)
    overlap_ids = None
    if overlap_ids_path is not None:
        overlap_ids_path = Path(overlap_ids_path)
        overlap_ids = {line.strip() for line in overlap_ids_path.read_text().splitlines() if line.strip()}
    frozen = _frozen_scores(eval_dir)
    model_names = sorted(frozen) + [external_name]
    seeds = sorted(next(iter(frozen.values())))
    for model, by_seed in frozen.items():
        if sorted(by_seed) != seeds:
            raise ValueError(f"{model} seed sets differ: {sorted(by_seed)} vs {seeds}")
    contract = {
        "schema": 1,
        "purpose": "External pretrained comparator on the frozen split (D0-31)",
        "external_name": external_name,
        "eval_dir": str(eval_dir), "units": str(units_path),
        "external_scores_sha256": sha256(external_scores_path),
        "units_sha256": sha256(units_path),
        "eval_predictions_sha256": sha256(eval_dir / "amr_predictions.tsv"),
        "overlap_ids": str(overlap_ids_path) if overlap_ids_path else None,
        "overlap_ids_sha256": sha256(overlap_ids_path) if overlap_ids_path else None,
        "frozen_analysis_set": "primary",
        "frozen_seeds": seeds,
        "bootstrap": f"near-clone-group resampling, {bootstrap_repeats} repeats; external scores are constant across seeds",
        "code_sha256": sha256(__file__),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    prepare_output(output, contract)
    drugs = sorted({unit["drug"] for unit in units if unit["genome_id"] in external})
    if not drugs:
        raise ValueError("No units carry an external score")
    metric_rows, pair_rows, prediction_rows, summaries = [], [], [], {}
    for drug in drugs:
        drug_units = [u for u in units if u["drug"] == drug]
        summaries[drug] = {}
        for set_name, keep in _score_sets(drug_units, overlap_ids).items():
            selected = [u for u in drug_units if keep(u)]
            test = sorted((u for u in selected if u["split"] == "test"), key=lambda u: u["genome_id"])
            missing = [u["genome_id"] for u in test if u["genome_id"] not in external]
            if missing:
                raise ValueError(f"{drug}/{set_name}: {len(missing)} test units lack external scores, e.g. {missing[:5]}")
            ids = [u["genome_id"] for u in test]
            labels = [1 if u["label"] == "R" else 0 for u in test]
            groups = [u["near_clone_group"] for u in test]
            if len(set(labels)) < 2:
                summaries[drug][set_name] = {"n_test": len(test), "status": "single_class", "labels": labels[:1]}
                continue
            per_seed_scores = {}
            for model in sorted(frozen):
                per_seed_scores[model] = [[frozen[model][s][g] for g in ids] for s in seeds]
                point_auroc, point_ap = _metrics(labels, _mean_over_seeds(per_seed_scores[model]))
                bootstrap = _group_bootstrap_multi(groups, labels, per_seed_scores[model],
                                                   bootstrap_repeats, seed)
                metric_rows.append({
                    "drug": drug, "analysis_set": set_name, "model": model, "n_test": len(test),
                    "test_S": labels.count(0), "test_R": labels.count(1),
                    "AUROC": point_auroc, "average_precision": point_ap,
                    "AUROC_ci95": json.dumps(bootstrap["AUROC_ci95"]),
                    "AP_ci95": json.dumps(bootstrap["AP_ci95"]),
                    "n_groups": bootstrap["n_groups"], "n_valid_repeats": bootstrap["n_valid_repeats"],
                })
            external_scores = [external[g] for g in ids]
            point_auroc, point_ap = _metrics(labels, external_scores)
            bootstrap = _group_bootstrap_multi(groups, labels, [external_scores], bootstrap_repeats, seed)
            metric_rows.append({
                "drug": drug, "analysis_set": set_name, "model": external_name, "n_test": len(test),
                "test_S": labels.count(0), "test_R": labels.count(1),
                "AUROC": point_auroc, "average_precision": point_ap,
                "AUROC_ci95": json.dumps(bootstrap["AUROC_ci95"]),
                "AP_ci95": json.dumps(bootstrap["AP_ci95"]),
                "n_groups": bootstrap["n_groups"], "n_valid_repeats": bootstrap["n_valid_repeats"],
            })
            for unit, score in zip(test, external_scores):
                prediction_rows.append({"drug": drug, "analysis_set": set_name,
                                        "genome_id": unit["genome_id"],
                                        "near_clone_group": unit["near_clone_group"],
                                        "label": 1 if unit["label"] == "R" else 0, "score": score})
            pairs = list(combinations(model_names, 2))
            for model_a, model_b in pairs:
                per_seed = []
                for index, _ in enumerate(seeds):
                    scores_a = (per_seed_scores[model_a][index] if model_a in per_seed_scores
                                else [external[g] for g in ids])
                    scores_b = (per_seed_scores[model_b][index] if model_b in per_seed_scores
                                else [external[g] for g in ids])
                    per_seed.append({model_a: scores_a, model_b: scores_b})
                auroc_deltas, ap_deltas = [], []
                for scores_by_model in per_seed:
                    a = _metrics(labels, scores_by_model[model_a])
                    b = _metrics(labels, scores_by_model[model_b])
                    if a[0] is None or b[0] is None:
                        continue
                    auroc_deltas.append(a[0] - b[0])
                    ap_deltas.append(a[1] - b[1])
                if not auroc_deltas:
                    continue
                paired = _paired_group_bootstrap_multi(groups, labels, per_seed, bootstrap_repeats,
                                                       seed, names=(model_a, model_b))
                pair_rows.append({
                    "drug": drug, "analysis_set": set_name, "model_a": model_a, "model_b": model_b,
                    "n_test_units": len(test),
                    "delta_AUROC_mean": sum(auroc_deltas) / len(auroc_deltas),
                    "delta_AUROC_sd": _sd(auroc_deltas), "delta_AUROC_ci95": json.dumps(paired["delta_AUROC_ci95"]),
                    "delta_AP_mean": sum(ap_deltas) / len(ap_deltas),
                    "delta_AP_sd": _sd(ap_deltas), "delta_AP_ci95": json.dumps(paired["delta_AP_ci95"]),
                })
            summaries[drug][set_name] = {
                "n_test": len(test), "n_S": labels.count(0), "n_R": labels.count(1),
                "n_groups": len({u["near_clone_group"] for u in test}),
                "models": [name for name in model_names],
                "overlap_examples": sorted(
                    {u["genome_id"] for u in drug_units
                     if overlap_ids is not None and u["genome_id"] in overlap_ids}
                )[:5],
            }
    write_tsv(output / "external_metrics.tsv", metric_rows, METRIC_COLUMNS)
    write_tsv(output / "external_paired.tsv", pair_rows, PAIR_COLUMNS)
    write_tsv(output / "external_predictions.tsv", prediction_rows, PREDICTION_COLUMNS)
    summary = {
        "status": "amr_external_eval_complete",
        "external_name": external_name, "drugs": summaries, "frozen_seeds": seeds,
        "limits": [
            "The external model is pretrained on its own data and budget; this is not a controlled retrain comparison.",
            "primary_disjoint removes test genomes listed as present in the external model's training data; the primary_all set keeps them and is confounded.",
            "StrainAMR runs in model-only mode (no known-mechanism rule layer) for comparability with rule-free representation models.",
            "External predictions cover one released checkpoint; frozen heads average three optimization seeds.",
        ],
        "metrics_sha256": sha256(output / "external_metrics.tsv"),
        "paired_sha256": sha256(output / "external_paired.tsv"),
        "predictions_sha256": sha256(output / "external_predictions.tsv"),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in (
        "summary.json", "external_metrics.tsv", "external_paired.tsv", "external_predictions.tsv", "contract.json")})
    (output / "README.md").write_text(
        "# External comparator evaluation (D0-31)\n\n"
        f"- External model: {external_name}\n"
        f"- Frozen scores: amr-eval `primary` analysis under `{eval_dir}`\n"
        "- `primary_all`: label-conflict-free test units; `primary_disjoint`: additionally excludes test "
        "genomes present in the external model's training set.\n"
        "- Near-clone-group bootstrap; frozen heads average seeds, external scores are constant across seeds.\n"
    )
    return summary


def _mean_over_seeds(scores_by_seed):
    return [sum(values) / len(values) for values in zip(*scores_by_seed)]


def _sd(values):
    import numpy as np
    return float(np.std(values))
