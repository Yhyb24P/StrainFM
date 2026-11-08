"""Frozen-representation AMR head training and evaluation (D0-27/D0-28).

Backbones stay frozen; only the paper-style light head is trained per drug:
LayerNorm -> Dropout -> Linear. Model selection (learning rate, early stopping)
uses the validation partition only; the test partition is scored once per
frozen configuration. Group-level bootstrap resamples near-clone groups.
"""

import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output

METRIC_COLUMNS = ["drug", "model", "analysis_set", "seed", "n_train", "n_validation", "n_test",
                  "test_S", "test_R", "test_prevalence", "AUROC", "average_precision",
                  "val_AUROC", "selected_lr", "selected_epoch", "selected_threshold"]


def _load_vectors(embeddings_dir, model):
    import pyarrow.parquet as pq
    name = {"bacformer": "genome_embeddings.parquet", "esm_mean": "esm_mean_embeddings.parquet"}[model]
    path = Path(embeddings_dir) / name
    if not path.exists():
        raise FileNotFoundError(f"Missing embeddings: {path}")
    table = pq.read_table(path, columns=["genome_id", "embedding"])
    return {row["genome_id"]: row["embedding"] for row in table.to_pylist()}


def _dataset(units, vectors, include_conflicts):
    rows = [u for u in units if (include_conflicts or u["label_conflict"] != "True") and u["genome_id"] in vectors]
    missing = sum(1 for u in units if (include_conflicts or u["label_conflict"] != "True") and u["genome_id"] not in vectors)
    return rows, missing


def _train_head(x_train, y_train, x_val, y_val, lr, max_epochs, patience, seed):
    import numpy as np
    import torch
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(False)
    head = torch.nn.Sequential(torch.nn.LayerNorm(x_train.shape[1]), torch.nn.Dropout(0.1), torch.nn.Linear(x_train.shape[1], 1))
    optimizer = torch.optim.Adam(head.parameters(), lr=lr)
    loss_fn = torch.nn.BCEWithLogitsLoss()
    xt = torch.tensor(x_train, dtype=torch.float32)
    yt = torch.tensor(y_train, dtype=torch.float32).unsqueeze(1)
    xv = torch.tensor(x_val, dtype=torch.float32)
    best_state, best_score, best_epoch, stale = None, -1.0, 0, 0
    val_single_class = len(set(y_val)) < 2
    for epoch in range(1, max_epochs + 1):
        head.train()
        optimizer.zero_grad()
        loss = loss_fn(head(xt), yt)
        loss.backward()
        optimizer.step()
        head.eval()
        with torch.no_grad():
            scores = head(xv).squeeze(1).numpy()
        if val_single_class:
            score = -float(loss_fn(head(xv), torch.tensor(y_val, dtype=torch.float32).unsqueeze(1)))
        else:
            from sklearn.metrics import roc_auc_score
            score = float(roc_auc_score(y_val, scores))
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
    import numpy as np
    import torch
    with torch.no_grad():
        return head(torch.tensor(x, dtype=torch.float32)).squeeze(1).numpy()


def _metrics(labels, scores):
    from sklearn.metrics import average_precision_score, roc_auc_score
    if len(set(labels)) < 2:
        return None, float(average_precision_score(labels, scores)) if labels else None
    return float(roc_auc_score(labels, scores)), float(average_precision_score(labels, scores))


def _paired_group_bootstrap(groups, labels, scores_by_model, repeats, seed):
    import numpy as np
    unique = sorted(set(groups))
    rng = np.random.default_rng(seed)
    by_group = defaultdict(list)
    for index, group in enumerate(groups):
        by_group[group].append(index)
    names = sorted(scores_by_model)
    deltas_auroc, deltas_ap = [], []
    for _ in range(repeats):
        sampled = rng.choice(len(unique), size=len(unique), replace=True)
        indices = [i for g in sampled for i in by_group[unique[g]]]
        y = [labels[i] for i in indices]
        values = {}
        for name in names:
            values[name] = _metrics(y, [scores_by_model[name][i] for i in indices])
        if values[names[0]][0] is not None and values[names[1]][0] is not None:
            deltas_auroc.append(values[names[0]][0] - values[names[1]][0])
        if values[names[0]][1] is not None and values[names[1]][1] is not None:
            deltas_ap.append(values[names[0]][1] - values[names[1]][1])
    def interval(values):
        if not values:
            return None
        return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]
    return {"n_repeats": repeats, "n_groups": len(unique),
            "delta_AUROC_ci95": interval(deltas_auroc), "delta_AP_ci95": interval(deltas_ap)}


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
    return {"n_repeats": repeats, "n_groups": len(unique), "AUROC_ci95": interval(aurocs), "AP_ci95": interval(aps)}


def run_amr_eval(embeddings_dir, units_path, output, models=("bacformer", "esm_mean"),
                 lr_grid=(0.0001, 0.001, 0.01), max_epochs=100, patience=10, seeds=(20261108, 20261109, 20261110),
                 bootstrap_repeats=1000):
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
        "lr_grid": list(lr_grid),
        "max_epochs": max_epochs,
        "patience": patience,
        "seeds": list(seeds),
        "head": "LayerNorm -> Dropout(0.1) -> Linear(960,1); BCEWithLogits; Adam",
        "selection": "learning rate and epoch selected on validation AUROC; test scored once per seed",
        "bootstrap": f"group-level resampling of near-clone groups, {bootstrap_repeats} repeats",
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
        for analysis_set, include_conflicts in (("primary", False), ("with_conflicts", True)):
            rows, missing = _dataset(drug_units, vectors[models[0]], include_conflicts)
            by_split = defaultdict(list)
            for row in rows:
                by_split[row["split"]].append(row)
            result = {"n_units": len(rows), "n_missing_embeddings": missing,
                      "partitions": {s: len(v) for s, v in by_split.items()}}
            model_scores = {}
            for model in models:
                prediction_scores, run_rows = {}, []
                for seed in seeds:
                    train = by_split.get("train", [])
                    validation = by_split.get("validation", [])
                    test = by_split.get("test", [])
                    x_train = np.asarray([vectors[model][u["genome_id"]] for u in train], dtype=np.float32)
                    y_train = np.asarray([1 if u["label"] == "R" else 0 for u in train], dtype=np.float32)
                    x_val = np.asarray([vectors[model][u["genome_id"]] for u in validation], dtype=np.float32)
                    y_val = np.asarray([1 if u["label"] == "R" else 0 for u in validation], dtype=np.float32)
                    x_test = np.asarray([vectors[model][u["genome_id"]] for u in test], dtype=np.float32)
                    y_test = np.asarray([1 if u["label"] == "R" else 0 for u in test], dtype=np.float32)
                    best = None
                    for lr in lr_grid:
                        head, val_score, epoch = _train_head(x_train, y_train, x_val, y_val, lr, max_epochs, patience, seed)
                        if best is None or val_score > best[0]:
                            best = (val_score, lr, epoch, head)
                    val_score, lr, epoch, head = best
                    scores = _predict(head, x_test)
                    auroc, ap = _metrics(y_test, scores)
                    metric_rows.append({
                        "drug": drug, "model": model, "analysis_set": analysis_set, "seed": seed,
                        "n_train": len(train), "n_validation": len(validation), "n_test": len(test),
                        "test_S": int((y_test == 0).sum()), "test_R": int((y_test == 1).sum()),
                        "test_prevalence": float(y_test.mean()) if len(y_test) else "",
                        "AUROC": "" if auroc is None else auroc, "average_precision": ap,
                        "val_AUROC": val_score, "selected_lr": lr, "selected_epoch": epoch,
                        "selected_threshold": "",
                    })
                    for unit, score, label in zip(test, scores.tolist(), y_test.tolist()):
                        run_rows.append({"drug": drug, "model": model, "analysis_set": analysis_set,
                                         "seed": seed, "genome_id": unit["genome_id"],
                                         "near_clone_group": unit["near_clone_group"],
                                         "label": label, "score": score, "split": "test"})
                    prediction_scores[seed] = scores
                prediction_rows.extend(run_rows)
                # group bootstrap on the first seed's predictions (test set)
                if run_rows and analysis_set == "primary":
                    first = [r for r in run_rows if r["seed"] == seeds[0]]
                    summaries[drug].setdefault(analysis_set, {})[model] = {
                        "bootstrap": _group_bootstrap([r["near_clone_group"] for r in first],
                                                      [r["label"] for r in first],
                                                      [r["score"] for r in first],
                                                      bootstrap_repeats, seeds[0]),
                        "runs": [r for r in metric_rows if r["drug"] == drug and r["model"] == model and r["analysis_set"] == analysis_set],
                    }
            summaries[drug][analysis_set] = dict(summaries[drug].get(analysis_set, {}), **result)
    write_tsv(output / "amr_metrics.tsv", metric_rows, METRIC_COLUMNS)
    write_tsv(output / "amr_predictions.tsv", prediction_rows,
              ["drug", "model", "analysis_set", "seed", "genome_id", "near_clone_group", "label", "score", "split"])
    # Paired model comparison on identical test units (primary analysis).
    pred_index = defaultdict(lambda: defaultdict(dict))
    label_index = {}
    for row in prediction_rows:
        pred_index[(row["drug"], row["analysis_set"], row["seed"])][row["model"]][row["genome_id"]] = row["score"]
        label_index[(row["drug"], row["analysis_set"], row["genome_id"])] = row["label"]
    comparison_rows = []
    for drug in drugs:
        for analysis_set in ("primary", "with_conflicts"):
            aligned = []
            for seed in seeds:
                scores = pred_index.get((drug, analysis_set, seed), {})
                if len(scores) < 2:
                    continue
                common = sorted(set(scores[models[0]]) & set(scores[models[1]]))
                if not common:
                    continue
                y = [label_index[(drug, analysis_set, gid)] for gid in common]
                per_model = {model: {gid: scores[model][gid] for gid in common} for model in models}
                a = _metrics(y, [per_model[models[0]][g] for g in common])
                b = _metrics(y, [per_model[models[1]][g] for g in common])
                if a[0] is None or b[0] is None:
                    continue
                aligned.append((seed, common, y, {model: [per_model[model][g] for g in common] for model in models}, a, b))
            if not aligned:
                continue
            delta_auroc = [a[0] - b[0] for _, _, _, _, a, b in aligned]
            delta_ap = [a[1] - b[1] for _, _, _, _, a, b in aligned]
            import numpy as np
            seed0, common, y, scores_lists, _, _ = aligned[0]
            groups = [next(row["near_clone_group"] for row in prediction_rows
                           if row["drug"] == drug and row["analysis_set"] == analysis_set
                           and row["seed"] == seed0 and row["model"] == models[0] and row["genome_id"] == gid)
                      for gid in common]
            paired = _paired_group_bootstrap(groups, y,
                                             {model: scores_lists[model] for model in models},
                                             bootstrap_repeats, seed0)
            comparison_rows.append({
                "drug": drug, "analysis_set": analysis_set,
                "model_a": models[0], "model_b": models[1],
                "n_test_units": len(common),
                "delta_AUROC_mean": float(np.mean(delta_auroc)),
                "delta_AUROC_sd": float(np.std(delta_auroc)),
                "delta_AUROC_ci95": json.dumps(paired["delta_AUROC_ci95"]),
                "delta_AP_mean": float(np.mean(delta_ap)),
                "delta_AP_sd": float(np.std(delta_ap)),
                "delta_AP_ci95": json.dumps(paired["delta_AP_ci95"]),
            })
    write_tsv(output / "amr_model_comparison.tsv", comparison_rows, [
        "drug", "analysis_set", "model_a", "model_b", "n_test_units",
        "delta_AUROC_mean", "delta_AUROC_sd", "delta_AUROC_ci95",
        "delta_AP_mean", "delta_AP_sd", "delta_AP_ci95"])
    summary = {
        "status": "amr_eval_complete",
        "models": list(models), "seeds": list(seeds), "analysis_sets": ["primary", "with_conflicts"],
        "drugs": summaries,
        "metrics_sha256": sha256(output / "amr_metrics.tsv"),
        "predictions_sha256": sha256(output / "amr_predictions.tsv"),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "limits": [
            "Frozen backbones with a paper-style light head; no backbone training.",
            "Label evidence is paired-accession primary; direct-only is a sensitivity subset (D0-27).",
            "Source components may cross partitions (D0-28); group bootstrap respects near-clone groups only.",
            "Single-class test partitions cannot yield AUROC and are reported as such.",
        ],
    }
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in ("summary.json", "amr_metrics.tsv", "amr_predictions.tsv", "amr_model_comparison.tsv", "contract.json")})
    (output / "README.md").write_text(
        "# AMR-02 frozen-representation evaluation\n\n"
        "Frozen Bacformer-L and matched ESM++ genome vectors, near-clone-group-isolated pre-registered "
        "split (D0-27/D0-28), paper-style light head per drug. Primary analysis excludes mixed-label "
        "clone groups; `with_conflicts` is the sensitivity analysis.\n"
    )
    return summary
