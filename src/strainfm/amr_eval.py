"""Frozen-representation AMR head training and evaluation (D0-27/D0-28).

Backbones stay frozen; only the paper-style light head is trained per drug:
LayerNorm -> Dropout -> Linear. Model selection (learning rate, early stopping)
uses the validation partition only; the test partition is scored once per
frozen configuration. Three analysis sets are reported:
- primary: mixed-label clone groups excluded, paired-accession labels;
- with_conflicts: mixed-label clone groups included (sensitivity);
- direct_refit: evidence-restricted refit on direct file-accession labels only
  (D0-27 sensitivity; frozen partitions kept, no new split).
Group-level bootstrap resamples near-clone groups; seed-averaged statistics
average the three optimization seeds on the same resampled groups before the
percentile interval is taken.
"""

import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output

METRIC_COLUMNS = ["drug", "model", "analysis_set", "seed", "n_train", "n_validation", "n_test",
                  "test_S", "test_R", "test_prevalence", "AUROC", "average_precision",
                  "val_AUROC", "selected_lr", "selected_epoch", "selected_threshold"]
ANALYSIS_SETS = (
    ("primary", lambda unit: unit["label_conflict"] != "True"),
    ("with_conflicts", lambda unit: True),
    ("direct_refit", lambda unit: unit["label_conflict"] != "True" and unit["label_source"] == "file_accession"),
)


def _load_vectors(embeddings_dir, model):
    import pyarrow.parquet as pq
    name = {"bacformer": "genome_embeddings.parquet", "esm_mean": "esm_mean_embeddings.parquet"}[model]
    path = Path(embeddings_dir) / name
    if not path.exists():
        raise FileNotFoundError(f"Missing embeddings: {path}")
    table = pq.read_table(path, columns=["genome_id", "embedding"])
    return {row["genome_id"]: row["embedding"] for row in table.to_pylist()}


def _select_units(drug_units, vectors, keep):
    selected = [unit for unit in drug_units if keep(unit)]
    present = [unit for unit in selected if unit["genome_id"] in vectors]
    return present, len(selected) - len(present)


def _train_head(x_train, y_train, x_val, y_val, lr, max_epochs, patience, seed):
    import torch
    torch.manual_seed(seed)
    head = torch.nn.Sequential(torch.nn.LayerNorm(x_train.shape[1]), torch.nn.Dropout(0.1),
                               torch.nn.Linear(x_train.shape[1], 1))
    optimizer = torch.optim.Adam(head.parameters(), lr=lr)
    loss_fn = torch.nn.BCEWithLogitsLoss()
    xt = torch.tensor(x_train, dtype=torch.float32)
    yt = torch.tensor(y_train, dtype=torch.float32).unsqueeze(1)
    xv = torch.tensor(x_val, dtype=torch.float32)
    yv = torch.tensor(y_val, dtype=torch.float32).unsqueeze(1)
    best_state, best_score, best_epoch, stale = None, -1.0, 0, 0
    val_single_class = len(set(y_val)) < 2
    for epoch in range(1, max_epochs + 1):
        head.train()
        optimizer.zero_grad()
        loss_fn(head(xt), yt).backward()
        optimizer.step()
        head.eval()
        if val_single_class:
            with torch.no_grad():
                score = -float(loss_fn(head(xv), yv))
        else:
            from sklearn.metrics import roc_auc_score
            with torch.no_grad():
                score = float(roc_auc_score(y_val, head(xv).squeeze(1).numpy()))
        if score > best_score + 1e-9:
            best_score, best_epoch, stale = score, epoch, 0
            best_state = {k: v.clone() for k, v in head.state_dict().items()}
        else:
            stale += 1
            if stale >= patience:
                break
    head.load_state_dict(best_state)
    head.eval()
    return head, best_score, best_epoch


def _predict(head, x):
    import torch
    with torch.no_grad():
        return head(torch.tensor(x, dtype=torch.float32)).squeeze(1).numpy()


def _metrics(labels, scores):
    from sklearn.metrics import average_precision_score, roc_auc_score
    if len(set(labels)) < 2:
        return None, float(average_precision_score(labels, scores)) if labels else None
    return float(roc_auc_score(labels, scores)), float(average_precision_score(labels, scores))


def _group_bootstrap(groups, labels, scores, repeats, seed):
    import numpy as np
    unique = sorted(set(groups))
    rng = np.random.default_rng(seed)
    aurocs, aps = [], []
    by_group = defaultdict(list)
    for group, label, score in zip(groups, labels, scores):
        by_group[group].append((label, score))
    for _ in range(repeats):
        sampled = rng.choice(len(unique), size=len(unique), replace=True)
        labels_b, scores_b = [], []
        for index in sampled:
            for label, score in by_group[unique[index]]:
                labels_b.append(label)
                scores_b.append(score)
        auroc, ap = _metrics(labels_b, scores_b)
        if auroc is not None:
            aurocs.append(auroc)
        if ap is not None:
            aps.append(ap)
    def interval(values):
        if not values:
            return None
        return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]
    return {"n_repeats": repeats, "n_valid_repeats": len(aurocs), "n_groups": len(unique),
            "AUROC_ci95": interval(aurocs), "AP_ci95": interval(aps)}


def _group_indices(groups):
    by_group = defaultdict(list)
    for index, group in enumerate(groups):
        by_group[group].append(index)
    return sorted(by_group), by_group


def _group_bootstrap_multi(groups, labels, scores_by_seed, repeats, seed):
    """CI of the seed-averaged metric; groups resampled once, seeds averaged per replicate."""
    import numpy as np
    unique, by_group = _group_indices(groups)
    rng = np.random.default_rng(seed)
    aurocs, aps = [], []
    for _ in range(repeats):
        sampled = rng.choice(len(unique), size=len(unique), replace=True)
        indices = [i for g in sampled for i in by_group[unique[g]]]
        y = [labels[i] for i in indices]
        if len(set(y)) < 2:
            continue
        seed_aurocs, seed_aps = [], []
        for scores in scores_by_seed:
            auroc, ap = _metrics(y, [scores[i] for i in indices])
            if auroc is not None:
                seed_aurocs.append(auroc)
            if ap is not None:
                seed_aps.append(ap)
        if seed_aurocs:
            aurocs.append(float(np.mean(seed_aurocs)))
        if seed_aps:
            aps.append(float(np.mean(seed_aps)))
    def interval(values):
        if not values:
            return None
        return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]
    return {"n_repeats": repeats, "n_valid_repeats": len(aurocs), "n_groups": len(unique),
            "AUROC_ci95": interval(aurocs), "AP_ci95": interval(aps)}


def _paired_group_bootstrap_multi(groups, labels, scores_by_model_seed, repeats, seed, names):
    """CI of the seed-averaged paired delta (names[0] - names[1]) on identical groups."""
    import numpy as np
    if len(names) != 2:
        raise ValueError("Paired bootstrap requires exactly two model names")
    for scores in scores_by_model_seed:
        if sorted(scores) != sorted(names):
            raise ValueError("Model score keys differ from the requested comparison order")
    unique, by_group = _group_indices(groups)
    rng = np.random.default_rng(seed)
    aurocs, aps = [], []
    for _ in range(repeats):
        sampled = rng.choice(len(unique), size=len(unique), replace=True)
        indices = [i for g in sampled for i in by_group[unique[g]]]
        y = [labels[i] for i in indices]
        if len(set(y)) < 2:
            continue
        seed_deltas_auroc, seed_deltas_ap = [], []
        for scores_by_model in scores_by_model_seed:
            a = _metrics(y, [scores_by_model[names[0]][i] for i in indices])
            b = _metrics(y, [scores_by_model[names[1]][i] for i in indices])
            if a[0] is not None and b[0] is not None:
                seed_deltas_auroc.append(a[0] - b[0])
            if a[1] is not None and b[1] is not None:
                seed_deltas_ap.append(a[1] - b[1])
        if seed_deltas_auroc:
            aurocs.append(float(np.mean(seed_deltas_auroc)))
        if seed_deltas_ap:
            aps.append(float(np.mean(seed_deltas_ap)))
    def interval(values):
        if not values:
            return None
        return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]
    return {"n_repeats": repeats, "n_valid_repeats": len(aurocs), "n_groups": len(unique),
            "delta_AUROC_ci95": interval(aurocs), "delta_AP_ci95": interval(aps)}


def run_amr_eval(embeddings_dir, units_path, output, models=("bacformer", "esm_mean"),
                 lr_grid=(0.0001, 0.001, 0.01), max_epochs=100, patience=10,
                 seeds=(20261108, 20261109, 20261110), bootstrap_repeats=1000):
    embeddings_dir, units_path, output = map(Path, (embeddings_dir, units_path, output))
    if "raw" in output.resolve().parts:
        raise ValueError("Non-raw output directory required")
    if output.exists():
        raise ValueError("Use a new output directory")
    units = read_table(units_path)
    contract = {
        "schema": 1,
        "purpose": "Frozen-backbone AMR light-head evaluation per D0-27/D0-28",
        "models": list(models),
        "analysis_sets": {
            "primary": "mixed-label clone groups excluded; paired-accession labels",
            "with_conflicts": "mixed-label clone groups included (sensitivity)",
            "direct_refit": "evidence-restricted refit on direct file-accession labels; frozen partitions kept",
        },
        "lr_grid": list(lr_grid), "max_epochs": max_epochs, "patience": patience, "seeds": list(seeds),
        "head": "LayerNorm -> Dropout(0.1) -> Linear(960,1); BCEWithLogits; Adam",
        "optimizer_note": ("torch.optim.Adam without explicit weight_decay; the pinned BacBench AMR script uses "
                           "AdamW without explicit weight_decay, so this is a documented deviation, not an exact "
                           "paper reproduction claim"),
        "selection": "learning rate and epoch selected on validation AUROC; test scored once per seed",
        "bootstrap": (f"group-level resampling of near-clone groups, {bootstrap_repeats} repeats; "
                      "seed-averaged statistics use the same resampled groups across the three seeds"),
        "units_sha256": sha256(units_path),
        "code_sha256": sha256(__file__),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    prepare_output(output, contract)
    import numpy as np
    vectors = {model: _load_vectors(embeddings_dir, model) for model in models}
    drugs = sorted({u["drug"] for u in units})
    metric_rows, prediction_rows, summaries = [], [], {}
    for drug in drugs:
        drug_units = [u for u in units if u["drug"] == drug]
        summaries[drug] = {}
        for analysis_set, keep in ANALYSIS_SETS:
            selected, missing = _select_units(drug_units, vectors[models[0]], keep)
            by_split = defaultdict(list)
            for row in selected:
                by_split[row["split"]].append(row)
            if not by_split.get("train") or not by_split.get("validation") or not by_split.get("test"):
                summaries[drug][analysis_set] = {"n_units": len(selected), "n_missing_embeddings": missing,
                                                 "partitions": {s: len(v) for s, v in by_split.items()},
                                                 "status": "insufficient_partition"}
                continue
            analysis = {"n_units": len(selected), "n_missing_embeddings": missing,
                        "partitions": {s: len(v) for s, v in by_split.items()}}
            for model in models:
                seed_scores, run_rows = {}, []
                for seed in seeds:
                    train, validation, test = by_split["train"], by_split["validation"], by_split["test"]
                    x_train = np.asarray([vectors[model][u["genome_id"]] for u in train], dtype=np.float32)
                    y_train = np.asarray([1 if u["label"] == "R" else 0 for u in train], dtype=np.float32)
                    x_val = np.asarray([vectors[model][u["genome_id"]] for u in validation], dtype=np.float32)
                    y_val = np.asarray([1 if u["label"] == "R" else 0 for u in validation], dtype=np.float32)
                    x_test = np.asarray([vectors[model][u["genome_id"]] for u in test], dtype=np.float32)
                    y_test = np.asarray([1 if u["label"] == "R" else 0 for u in test], dtype=np.float32)
                    best = None
                    for lr in lr_grid:
                        head, val_score, epoch = _train_head(x_train, y_train, x_val, y_val, lr,
                                                             max_epochs, patience, seed)
                        if best is None or val_score > best[0]:
                            best = (val_score, lr, epoch, head)
                    val_score, lr, epoch, head = best
                    scores = _predict(head, x_test)
                    auroc, ap = _metrics(y_test, scores)
                    seed_scores[seed] = scores
                    metric_rows.append({
                        "drug": drug, "model": model, "analysis_set": analysis_set, "seed": seed,
                        "n_train": len(train), "n_validation": len(validation), "n_test": len(test),
                        "test_S": int((y_test == 0).sum()), "test_R": int((y_test == 1).sum()),
                        "test_prevalence": float(y_test.mean()) if len(y_test) else "",
                        "AUROC": "" if auroc is None else auroc, "average_precision": ap,
                        "val_AUROC": val_score, "selected_lr": lr, "selected_epoch": epoch,
                        "selected_threshold": "",
                    })
                    for unit, score in zip(test, scores.tolist()):
                        run_rows.append({"drug": drug, "model": model, "analysis_set": analysis_set,
                                         "seed": seed, "genome_id": unit["genome_id"],
                                         "near_clone_group": unit["near_clone_group"],
                                         "label": 1 if unit["label"] == "R" else 0,
                                         "score": score, "split": "test"})
                prediction_rows.extend(run_rows)
                groups = [unit["near_clone_group"] for unit in by_split["test"]]
                labels = [1 if unit["label"] == "R" else 0 for unit in by_split["test"]]
                analysis.setdefault("models", {})[model] = {
                    "bootstrap": _group_bootstrap_multi(groups, labels, [seed_scores[s] for s in seeds],
                                                        bootstrap_repeats, seeds[0]),
                    "runs": [r for r in metric_rows if r["drug"] == drug and r["model"] == model
                             and r["analysis_set"] == analysis_set],
                }
            summaries[drug][analysis_set] = analysis
    write_tsv(output / "amr_metrics.tsv", metric_rows, METRIC_COLUMNS)
    write_tsv(output / "amr_predictions.tsv", prediction_rows,
              ["drug", "model", "analysis_set", "seed", "genome_id", "near_clone_group", "label", "score", "split"])
    comparison_rows = []
    for drug in drugs:
        for analysis_set, _ in ANALYSIS_SETS:
            per_seed = []
            common = None
            labels = None
            groups = None
            for seed in seeds:
                rows = [r for r in prediction_rows if r["drug"] == drug and r["analysis_set"] == analysis_set and r["seed"] == seed]
                if not rows:
                    continue
                by_model = defaultdict(dict)
                meta = {}
                for row in rows:
                    by_model[row["model"]][row["genome_id"]] = row["score"]
                    meta[row["genome_id"]] = (row["label"], row["near_clone_group"])
                ids = sorted(set(by_model[models[0]]) & set(by_model[models[1]]))
                if not ids:
                    continue
                per_seed.append({model: [by_model[model][g] for g in ids] for model in models})
                if common is None:
                    common = ids
                    labels = [meta[g][0] for g in ids]
                    groups = [meta[g][1] for g in ids]
                elif ids != common:
                    raise ValueError("Seed prediction unit sets differ; refusing inconsistent comparison")
            if not per_seed:
                continue
            auroc_deltas, ap_deltas = [], []
            for scores_by_model in per_seed:
                a = _metrics(labels, scores_by_model[models[0]])
                b = _metrics(labels, scores_by_model[models[1]])
                if a[0] is None or b[0] is None:
                    continue
                auroc_deltas.append(a[0] - b[0])
                ap_deltas.append(a[1] - b[1])
            if not auroc_deltas:
                continue
            paired = _paired_group_bootstrap_multi(groups, labels, per_seed, bootstrap_repeats, seeds[0],
                                                   names=(models[0], models[1]))
            comparison_rows.append({
                "drug": drug, "analysis_set": analysis_set,
                "model_a": models[0], "model_b": models[1], "n_test_units": len(common),
                "delta_AUROC_mean": float(np.mean(auroc_deltas)),
                "delta_AUROC_sd": float(np.std(auroc_deltas)),
                "delta_AUROC_ci95": json.dumps(paired["delta_AUROC_ci95"]),
                "delta_AP_mean": float(np.mean(ap_deltas)),
                "delta_AP_sd": float(np.std(ap_deltas)),
                "delta_AP_ci95": json.dumps(paired["delta_AP_ci95"]),
            })
    write_tsv(output / "amr_model_comparison.tsv", comparison_rows, [
        "drug", "analysis_set", "model_a", "model_b", "n_test_units",
        "delta_AUROC_mean", "delta_AUROC_sd", "delta_AUROC_ci95",
        "delta_AP_mean", "delta_AP_sd", "delta_AP_ci95"])
    summary = {
        "status": "amr_eval_complete",
        "models": list(models), "seeds": list(seeds),
        "analysis_sets": [name for name, _ in ANALYSIS_SETS],
        "drugs": summaries,
        "metrics_sha256": sha256(output / "amr_metrics.tsv"),
        "predictions_sha256": sha256(output / "amr_predictions.tsv"),
        "comparison_sha256": sha256(output / "amr_model_comparison.tsv"),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "limits": [
            "Frozen backbones with a paper-style light head; no backbone training.",
            "Label evidence is paired-accession primary; direct_refit is an evidence-restricted refit (D0-27), not a gate pass.",
            "Source components may cross partitions (D0-28); group bootstrap respects near-clone groups only.",
            "Single-class test partitions cannot yield AUROC and are reported as such.",
            "Cross-drug partitions differ for shared BioSamples; a future shared trained encoder needs a unified split.",
            "Optimizer differs from the pinned BacBench AMR script (Adam here vs AdamW there, neither with explicit weight_decay); no exact paper reproduction is claimed.",
        ],
    }
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in (
        "summary.json", "amr_metrics.tsv", "amr_predictions.tsv", "amr_model_comparison.tsv", "contract.json")})
    (output / "README.md").write_text(
        "# AMR-02 frozen-representation evaluation (2026-11-09 completion run)\n\n"
        "Frozen Bacformer-L and matched ESM++ genome vectors, pre-registered near-clone-group split "
        "(D0-27/D0-28), paper-style light head per drug. `primary` excludes mixed-label clone groups; "
        "`with_conflicts` includes them; `direct_refit` refits on direct file-accession labels only "
        "(evidence-restricted sensitivity). Seed-averaged statistics use the same resampled groups "
        "across the three seeds before the percentile interval is taken.\n"
    )
    return summary
