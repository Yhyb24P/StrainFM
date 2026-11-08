"""Frozen-protocol strain retrieval evaluation (protocol A train-gallery, protocol B held-out-neighborhood).

Pure-stdlib implementation. Per the frozen 2026-11-06 contract:
- k set {1, 5, 10} by default; ranking by cosine_distance = 1 - cos_sim ascending,
  ties broken by genome_id lexicographic ascending.
- Denominator of Hit/Recall/MRR/mAP is positive-supported queries only;
  open-set (no gallery positive) queries are counted separately and never
  mixed into score averages.
- Stratification: seen_in_train_neighborhood x assembly_level (complete vs draft).
- UMAP/clustering figures are diagnostics only and are not produced here.
"""

import json
import math
import random
from collections import defaultdict
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output

DEFAULT_K = (1, 5, 10)


def cosine_distance(a, b):
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    if na == 0.0 or nb == 0.0:
        raise ValueError("Zero-norm embedding")
    return 1.0 - dot / (na * nb)


def mean_average_precision(ranked_positive_flags, k, n_positives):
    """AP@k with denominator min(k, n_positives).

    Missed positives below the cut-off reduce the score; dividing by the number
    of hits (the previous implementation) would report 1.0 for a query whose
    second positive ranks last inside the top-k. Definition frozen in D0-24.
    """
    denominator = min(k, n_positives)
    if denominator <= 0:
        return 0.0
    hits = 0
    total = 0.0
    for rank, flag in enumerate(ranked_positive_flags[:k], start=1):
        if flag:
            hits += 1
            total += hits / rank
    return total / denominator


def median(values):
    if not values:
        return ""
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def _normalized(vector):
    if not vector:
        raise ValueError("Empty embedding")
    if not all(math.isfinite(x) for x in vector):
        raise ValueError("Non-finite embedding value")
    norm = math.sqrt(sum(x * x for x in vector))
    if norm == 0.0:
        raise ValueError("Zero-norm embedding")
    inv = 1.0 / norm
    return [x * inv for x in vector]


def load_embeddings(embeddings_dir, model):
    """Return ({genome_id: normalized_vector}, {genome_id: truncated_protein_count}).

    Reads the fixed artifact layout: parquet when the file exists (Bacformer
    layout carries truncated_protein_count), otherwise the TSV fallback with the
    same columns (embedding as a JSON list).
    """
    directory = Path(embeddings_dir)
    parquet = {"bacformer": directory / "genome_embeddings.parquet",
               "esm_mean": directory / "esm_mean_embeddings.parquet"}[model]
    tsv = {"bacformer": directory / "genome_embeddings.tsv",
           "esm_mean": directory / "esm_mean_embeddings.tsv"}[model]
    if parquet.exists():
        import pyarrow.parquet as pq
        columns = ["genome_id", "embedding"]
        if model == "bacformer":
            columns.append("truncated_protein_count")
        rows = pq.read_table(parquet, columns=columns).to_pylist()
    elif tsv.exists():
        raw_rows = list(read_table(tsv))
        rows = [{"genome_id": r["genome_id"],
                 "embedding": json.loads(r["embedding"]),
                 "truncated_protein_count": int(r.get("truncated_protein_count") or 0)}
                for r in raw_rows]
    else:
        raise ValueError(f"No embedding artifact for {model}: tried {parquet} and {tsv}")
    vectors, truncated = {}, {}
    dimension = None
    for row in rows:
        if row["genome_id"] in vectors:
            raise ValueError(f"Duplicate genome_id in embeddings: {row['genome_id']}")
        vector = _normalized(row["embedding"])
        if dimension is None:
            dimension = len(vector)
        elif len(vector) != dimension:
            raise ValueError(f"Inconsistent embedding dimension for {row['genome_id']}: {len(vector)} != {dimension}")
        vectors[row["genome_id"]] = vector
        truncated[row["genome_id"]] = int(row.get("truncated_protein_count") or 0)
    return vectors, truncated


def failed_genomes(embeddings_dir):
    path = Path(embeddings_dir) / "embedding_failures.tsv"
    if not path.exists():
        return set()
    return {r["genome_id"] for r in read_table(path)}


def quality_bin(assembly_level):
    return "complete" if assembly_level in ("Complete Genome", "Chromosome") else "draft"


def _select_protocol_a(manifests, neighborhood_of, train_ids, embedded, failures, split):
    """Queries = usable genomes of the requested split only; positives = train
    gallery of the same sequence neighborhood. Open-set queries (no positives)
    are kept and identified by an empty positive set. Queries present in the
    split but absent from the embeddings (or failed) are reported, not dropped
    silently."""
    usable_train = [g for g in train_ids if g in embedded and g not in failures]
    by_neighborhood = defaultdict(list)
    for gid in usable_train:
        by_neighborhood[neighborhood_of.get(gid, "")].append(gid)
    queries, missing, failed_in_split = {}, [], []
    for gid, row in sorted(manifests.items()):
        if row["split"] != split:
            continue
        if gid in failures:
            failed_in_split.append(gid)
            continue
        if gid not in embedded:
            missing.append(gid)
            continue
        positives = sorted(by_neighborhood.get(neighborhood_of.get(gid, ""), ()))
        queries[gid] = {"positives": set(positives), "gallery": usable_train, "protocol": "A"}
    return queries, missing, failed_in_split


def _select_protocol_b(manifests, neighborhood_of, split, embedded, failures, seed):
    """Held-out-neighborhood internal retrieval: neighborhoods of the target split
    with >= 2 usable members; first member after a fixed-seed shuffle is the
    query, the rest form the gallery (all positives). Nothing here trains or tunes."""
    by_neighborhood = defaultdict(list)
    for gid, row in sorted(manifests.items()):
        if row["split"] != split or gid not in embedded or gid in failures:
            continue
        by_neighborhood[neighborhood_of.get(gid, "")].append(gid)
    queries = {}
    for neighborhood, members in sorted(by_neighborhood.items()):
        if len(members) < 2:
            continue
        shuffled = members[:]
        random.Random(f"{seed}\tB\t{split}\t{neighborhood}").shuffle(shuffled)
        queries[shuffled[0]] = {"positives": set(shuffled[1:]), "gallery": sorted(shuffled[1:]), "protocol": "B"}
    return queries


def _rank_gallery(query_vector, gallery_ids, vectors):
    scored = []
    for gid in gallery_ids:
        v = vectors[gid]
        dist = 1.0 - sum(a * b for a, b in zip(query_vector, v))
        scored.append((dist, gid))
    scored.sort()
    return [gid for _, gid in scored]


def run_eval_retrieval(embeddings_dir, split_manifest, support, groups, output,
                       model="bacformer", protocol="A", split="validation",
                       k_values=DEFAULT_K, seed=20261106,
                       model_name="bacformer-large-frozen"):
    embeddings_dir, split_manifest, support, groups, output = map(Path, (embeddings_dir, split_manifest, support, groups, output))
    if output.exists() or "raw" in output.resolve().parts:
        raise ValueError("Use a new non-raw output directory")
    k_values = tuple(sorted(set(int(k) for k in k_values)))
    if not k_values or min(k_values) < 1:
        raise ValueError("k values must be positive integers")
    if protocol not in ("A", "B"):
        raise ValueError("protocol must be A or B")
    if protocol == "A" and split not in ("validation", "test"):
        raise ValueError("Protocol A supports --split validation or test")
    if protocol == "B" and split not in ("validation", "test"):
        raise ValueError("Protocol B supports --split validation or test")
    manifests = {r["genome_id"]: r for r in read_table(split_manifest)}
    neighborhood_of = {r["genome_id"]: r["sequence_neighborhood"] for r in read_table(groups)}
    train_neighborhoods = {neighborhood_of[r["genome_id"]] for r in manifests.values() if r["split"] == "train"}
    support_rows = {r["genome_id"]: r for r in read_table(support)}
    vectors, truncated = load_embeddings(embeddings_dir, model)
    failures = failed_genomes(embeddings_dir)
    missing_queries, failed_in_split, support_eligible_no_gallery_positive = [], [], 0
    support_positive_count_mismatch, support_ineligible_with_positives = 0, 0
    if protocol == "A":
        train_ids = [gid for gid, row in manifests.items() if row["split"] == "train"]
        queries, missing_queries, failed_in_split = _select_protocol_a(
            manifests, neighborhood_of, train_ids, vectors, failures, split)
        # Per-query cross-check against the frozen support file. Counts are
        # reported; queries are never silently re-labelled.
        for gid, spec in queries.items():
            row = support_rows.get(gid)
            if row is None:
                continue
            eligible = row.get("eligible_for_positive_supported_retrieval", "").strip().lower() in ("true", "1", "yes")
            expected = int(row.get("gallery_positive_count") or 0)
            actual = len(spec["positives"])
            if eligible and actual == 0:
                support_eligible_no_gallery_positive += 1
            if eligible and actual != expected:
                support_positive_count_mismatch += 1
            if not eligible and actual > 0:
                support_ineligible_with_positives += 1
    else:
        queries = _select_protocol_b(manifests, neighborhood_of, split, vectors, failures, seed)
    positive_definition = "same_Mash_sequence_neighborhood" if protocol == "A" else "held_out_neighborhood_gallery"
    contract = {
        "schema": 1,
        "purpose": "Frozen-protocol strain retrieval evaluation",
        "model": model,
        "model_name": model_name,
        "protocol": protocol,
        "split": split,
        "gallery_split": "train" if protocol == "A" else "within_neighborhood",
        "k_values": list(k_values),
        "seed": seed,
        "positive_definition": positive_definition,
        "ranking": "cosine_distance ascending; ties genome_id lexicographic ascending",
        "open_set_rule": "queries without gallery positives are counted in open_set_n and excluded from all score denominators",
        "embeddings_dir": str(embeddings_dir),
        "split_manifest_sha256": sha256(split_manifest),
        "support_sha256": sha256(support),
        "groups_sha256": sha256(groups),
        "n_embedding_failures": len(failures),
        "n_queries_missing_embeddings": len(missing_queries),
        "n_query_split_failures": len(failed_in_split),
        "missing_query_ids": missing_queries[:50],
        "failed_query_ids": failed_in_split[:50],
        "code_sha256": sha256(__file__),
    }
    if protocol == "A":
        contract["support_eligible_no_gallery_positive"] = support_eligible_no_gallery_positive
        contract["support_positive_count_mismatch"] = support_positive_count_mismatch
        contract["support_ineligible_with_positives"] = support_ineligible_with_positives
    prepare_output(output, contract)
    kmax = max(k_values)
    scored_records, strata = [], defaultdict(lambda: {"total": 0, "open": 0, "records": []})
    for gid in sorted(queries):
        spec = queries[gid]
        row = manifests[gid]
        seen = int(neighborhood_of.get(gid, "") in train_neighborhoods)
        quality = quality_bin(row.get("assembly_level", ""))
        cell = strata[(seen, quality)]
        cell["total"] += 1
        if not spec["positives"]:
            cell["open"] += 1
            continue
        ranked = _rank_gallery(vectors[gid], spec["gallery"], vectors)
        positives = spec["positives"]
        flags = [g in positives for g in ranked]
        first = next((rank for rank, flag in enumerate(flags, start=1) if flag), None)
        record = {
            "model_name": model_name,
            "query_genome_id": gid,
            "protocol": protocol,
            "k": kmax,
            "hit_at_k": int(any(flags[:kmax])),
            "recall_at_k": sum(flags[:kmax]) / len(positives),
            "ap_at_k": mean_average_precision(flags, 10, len(positives)) if kmax >= 10 else "",
            "rr": (1.0 / first) if first else 0.0,
            "first_positive_rank": first,
            "n_positive": len(positives),
            "n_gallery": len(spec["gallery"]),
            "seen_in_train_neighborhood": seen,
            "quality_bin": quality,
            "truncated": truncated.get(gid, 0),
        }
        for k in k_values:
            hits = sum(flags[:k])
            record[f"hit@{k}"] = int(hits > 0)
            record[f"recall@{k}"] = hits / len(positives)
        scored_records.append(record)
        cell["records"].append(record)
    n_total = sum(c["total"] for c in strata.values())
    n_open = sum(c["open"] for c in strata.values())
    metric_rows = []
    for (seen, quality), cell in sorted(strata.items(), key=lambda kv: (str(kv[0][0]), kv[0][1])):
        records = cell["records"]
        metric_rows.append({
            "model_name": model_name, "protocol": protocol, "split": split,
            "stratum": f"seen={seen};{quality}",
            "n_queries_total": cell["total"], "open_set_n": cell["open"], "n_queries_scored": len(records),
            **_score_block(records, k_values),
        })
    metric_rows.append({
        "model_name": model_name, "protocol": protocol, "split": split, "stratum": "ALL",
        "n_queries_total": n_total, "open_set_n": n_open, "n_queries_scored": len(scored_records),
        **_score_block(scored_records, k_values),
    })
    results_columns = ["model_name", "query_genome_id", "protocol", "k", "hit_at_k", "recall_at_k",
                       "ap_at_k", "rr", "first_positive_rank", "n_positive", "n_gallery",
                       "seen_in_train_neighborhood", "quality_bin", "truncated"]
    results_columns += [f"hit@{k}" for k in k_values] + [f"recall@{k}" for k in k_values]
    write_tsv(output / "retrieval_metrics.tsv", metric_rows, _metrics_columns(k_values))
    _write_results(output, scored_records, results_columns, k_values)
    summary = {
        "status": "retrieval_eval_complete",
        "model": model, "model_name": model_name, "protocol": protocol, "split": split,
        "n_queries_total": n_total, "open_set_n": n_open, "n_queries_scored": len(scored_records),
        "n_strata": len(strata), "k_values": list(k_values), "seed": seed,
        "n_embedding_failures": len(failures),
        "n_queries_missing_embeddings": len(missing_queries),
        "n_query_split_failures": len(failed_in_split),
        "all_hit": {f"hit@{k}": _mean(scored_records, f"hit@{k}") for k in k_values},
        "all_recall": {f"recall@{k}": _mean(scored_records, f"recall@{k}") for k in k_values},
        "all_mrr": _mean(scored_records, "rr"),
        "all_map@10": _mean(scored_records, "ap_at_k") if kmax >= 10 else "",
    }
    atomic_json(output / "summary.json", summary)
    return summary


def _score_block(records, k_values):
    block = {}
    for k in k_values:
        block[f"hit@{k}"] = _mean(records, f"hit@{k}")
        block[f"recall@{k}"] = _mean(records, f"recall@{k}")
    block["mrr"] = _mean(records, "rr")
    if max(k_values) >= 10:
        block["map@10"] = _mean_values([r["ap_at_k"] for r in records if r["ap_at_k"] != ""])
    else:
        block["map@10"] = ""
    block["median_positive_rank"] = median([r["first_positive_rank"] for r in records if r["first_positive_rank"]])
    return block


def _mean(records, field):
    return _mean_values([r[field] for r in records if r[field] != ""])


def _mean_values(values):
    if not values:
        return ""
    return sum(values) / len(values)


def _metrics_columns(k_values):
    columns = ["model_name", "protocol", "split", "stratum",
               "n_queries_total", "open_set_n", "n_queries_scored"]
    for k in k_values:
        columns.append(f"hit@{k}")
    for k in k_values:
        columns.append(f"recall@{k}")
    columns += ["mrr", "map@10", "median_positive_rank"]
    return columns


def _write_results(output, records, columns, k_values):
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
        table = pa.table({c: [r.get(c, "") for r in records] for c in columns})
        pq.write_table(table, output / "retrieval_results.parquet")
        return "parquet"
    except ImportError:
        write_tsv(output / "retrieval_results.tsv", records, columns)
        return "tsv"
