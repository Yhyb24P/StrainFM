"""Distance-consistency evaluation: model embedding distance vs Mash genetic distance.

Frozen 2026-11-06 contract:
- For each query (validation or test) sample with a fixed seed:
  <=20 same-neighborhood positives (train gallery), <=20 same-split
  different-neighborhood negatives, <=20 distance-binned candidates
  (Mash within +/-0.005 of the query's closest positive distance),
  <=20 far negatives (Mash >= 0.10 from that reference distance).
- Spearman correlation between embedding cosine_distance and mash_distance on
  the sampled pairs; no partial Mantel deconfounding is claimed.
- Strata: same_coarse_neighborhood 0/1, same_near_clone_group 0/1,
  complete vs draft (query assembly level).

Distance sources, in order:
1. The stored single-linkage neighborhoods forest (edges within the 0.005
   threshold, genome-id keyed) — authoritative and free for same-neighborhood
   and distance-binned pairs; non-adjacent pairs inside a neighborhood
   component may lack a stored edge and are reported as unresolved, never
   imputed.
2. The mash CLI run against the combined sketch for far-bin pairs. Sketch
   files are keyed by full file path, so the caller supplies the
   genome_id -> fasta path map (assemblies.tsv target_relpath under
   genome_root); the CLI is batched per query over a bounded candidate pool.
"""

import gzip
import random
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output
from .eval_retrieval import load_embeddings, failed_genomes, quality_bin

PAIR_COLUMNS = [
    "genome_id_a", "genome_id_b", "split_pair", "mash_distance",
    "same_near_clone_group", "same_coarse_neighborhood",
    "embedding_cosine_distance", "method", "sampling_bin",
]
CORRELATION_COLUMNS = ["model_name", "split", "stratum", "n_pairs", "spearman_rho", "method"]


def load_forest_distances(forest_path):
    """Undirected forest edges -> {(a, b): distance} plus per-node adjacency."""
    distances, adjacency = {}, defaultdict(set)
    with gzip.open(forest_path, "rt") as handle:
        lines = handle.readlines()
    for line in lines[1:]:
        a, b, d, _shared = line.rstrip("\n").split("\t")
        d = float(d)
        key = (a, b) if a <= b else (b, a)
        distances[key] = min(distances.get(key, d), d)
        adjacency[a].add(b)
        adjacency[b].add(a)
    return distances, adjacency


def mash_cli_batch(mash_binary, sketch, candidate_paths, workers=8):
    """Reference-less one-shot `mash dist <reference-sketch> <list> -l`:
    distances of every listed file against the sketch, keyed by listed path.
    Returns {path: distance} for resolved entries; empty dict on failure."""
    # Whole-file mode (no -i): each listed file is one query against the
    # reference sketch; these genomes are single-assembly files.
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as handle:
        handle.write("\n".join(candidate_paths) + "\n")
        list_file = handle.name
    try:
        result = subprocess.run(
            [str(mash_binary), "dist", str(sketch), list_file, "-l", "-t", "-p", str(workers)],
            capture_output=True, text=True, check=False)
        if result.returncode != 0:
            return {}
        # List-mode rows are "<query path>\t<distance>" (the reference column
        # header is skipped).
        resolved = {}
        for line in result.stdout.splitlines():
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            try:
                resolved[fields[0]] = float(fields[1])
            except ValueError:
                continue
        return resolved
    finally:
        Path(list_file).unlink(missing_ok=True)


def _ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for t in range(i, j + 1):
            ranks[order[t]] = avg
        i = j + 1
    return ranks


def spearman_rho(xs, ys):
    if len(xs) != len(ys) or len(xs) < 3:
        return None
    rx, ry = _ranks(list(xs)), _ranks(list(ys))
    mx = sum(rx) / len(rx)
    my = sum(ry) / len(ry)
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = sum((a - mx) ** 2 for a in rx)
    vy = sum((b - my) ** 2 for b in ry)
    if vx == 0 or vy == 0:
        return None
    import math
    return cov / math.sqrt(vx * vy)


def sample_distance_pairs(query_gid, positive_dists, binned_pool, rng, per_bin=20, bin_width=0.005):
    """Pick the reference distance (closest positive) and distance-binned
    partners from a pool of (partner, distance) already resolvable via forest."""
    reference = min(positive_dists) if positive_dists else None
    selected = []
    if reference is not None:
        in_width = [p for p in binned_pool if abs(p[1] - reference) <= bin_width]
        if len(in_width) > per_bin:
            in_width = rng.sample(in_width, per_bin)
        selected = sorted(in_width)
    return reference, selected


def run_eval_distance(embeddings_dir, split_manifest, groups_dir, assemblies, genome_root,
                      mash_binary, output, model="bacformer", split="validation",
                      seed=20261106, model_name="bacformer-large-frozen",
                      per_bin=20, max_far_candidates_per_query=400, cli_workers=8):
    embeddings_dir, split_manifest, groups_dir, assemblies, genome_root, mash_binary, output = map(
        Path, (embeddings_dir, split_manifest, groups_dir, assemblies, genome_root, mash_binary, output))
    if output.exists() or "raw" in output.resolve().parts:
        raise ValueError("Use a new non-raw output directory")
    if split not in ("validation", "test"):
        raise ValueError("--split must be validation or test")
    if per_bin < 1:
        raise ValueError("per_bin must be >= 1")
    if not Path(mash_binary).exists():
        raise ValueError(f"mash binary not found: {mash_binary}")
    groups_path = groups_dir / "sequence_groups.tsv"
    forest_path = groups_dir / "neighborhood_forest.tsv.gz"
    sketch = groups_dir / "genomes.msh"
    if not forest_path.exists():
        raise ValueError(f"Missing neighborhoods forest: {forest_path}")
    if not sketch.exists():
        raise ValueError(f"Missing combined sketch: {sketch}")
    group_rows = {r["genome_id"]: r for r in read_table(groups_path)}
    manifests = {r["genome_id"]: r for r in read_table(split_manifest)}
    assembly_rows = {r["genome_id"]: r for r in read_table(assemblies)}
    forest_distances, forest_adjacency = load_forest_distances(forest_path)
    vectors, _ = load_embeddings(embeddings_dir, model)
    failures = failed_genomes(embeddings_dir)
    neighborhood_of = {gid: row["sequence_neighborhood"] for gid, row in group_rows.items()}
    near_clone_of = {gid: row["near_clone_group"] for gid, row in group_rows.items()}

    def forest_dist(a, b):
        key = (a, b) if a <= b else (b, a)
        return forest_distances.get(key)

    def path_of(gid):
        rel = assembly_rows.get(gid, {}).get("target_relpath")
        return str(genome_root / rel) if rel else None

    train_ids = [gid for gid, row in manifests.items() if row["split"] == "train"]
    train_set = set(train_ids)
    by_neighborhood = defaultdict(list)
    for gid in train_ids:
        by_neighborhood[neighborhood_of.get(gid, "")].append(gid)
    usable = [gid for gid, row in manifests.items()
              if row["split"] == split and gid in vectors and gid not in failures]
    query_ids = sorted(usable)
    pair_rows = []
    n_cli_batches = 0
    n_unresolved = 0
    for query in query_ids:
        rng = random.Random(f"{seed}\tdistance\t{query}")
        q_neighborhood = neighborhood_of.get(query, "")
        positives = [g for g in sorted(by_neighborhood.get(q_neighborhood, ()))
                     if g in vectors and g not in failures]
        positive_dists = [forest_dist(query, g) for g in positives]
        positive_dists = [d for d in positive_dists if d is not None]
        same_split = [g for g in usable if g != query and neighborhood_of.get(g, "") != q_neighborhood]

        # same-split negatives: forest-resolved where possible; the rest are
        # distance-less and are dropped (counted), never imputed
        neg_pool = [(g, forest_dist(query, g)) for g in sorted(same_split)]
        neg_resolved = [p for p in neg_pool if p[1] is not None]
        neg_missing = len(neg_pool) - len(neg_resolved)
        sampled_neg = neg_resolved if len(neg_resolved) <= per_bin else rng.sample(neg_resolved, per_bin)
        sampled_neg.sort()
        reference, binned = sample_distance_pairs(query, positive_dists, neg_resolved, rng, per_bin)
        if len(binned) > per_bin:
            binned = sorted(rng.sample(binned, per_bin))

        # far negatives: bounded candidate pool resolved in one CLI batch per query
        far_resolved = {}
        if reference is not None:
            far_pool = sorted(rng.sample(same_split, min(max_far_candidates_per_query, len(same_split))))
            far_paths = [(g, path_of(g)) for g in far_pool]
            far_paths = [(g, p) for g, p in far_paths if p]
            if far_paths:
                n_cli_batches += 1
                result = mash_cli_batch(mash_binary, sketch, [p for _, p in far_paths], workers=cli_workers)
                far_resolved = {g: result[p] for g, p in far_paths if p in result}
                far_missing = len(far_paths) - len([g for g, _ in far_paths if _ in result])
                n_unresolved += far_missing
            far_selected = sorted([g for g, d in far_resolved.items() if d >= 0.10])[:per_bin]
        else:
            far_selected = []

        def emit(partner, mdist, bin_name):
            nonlocal n_unresolved
            if mdist is None:
                n_unresolved += 1
                return
            a_vec, b_vec = vectors[query], vectors[partner]
            cos_dist = 1.0 - sum(x * y for x, y in zip(a_vec, b_vec))
            pair_rows.append({
                "genome_id_a": query,
                "genome_id_b": partner,
                "split_pair": f"{split}-train" if partner in train_set else split,
                "mash_distance": mdist,
                "same_near_clone_group": int(near_clone_of.get(query, "") == near_clone_of.get(partner, "") and near_clone_of.get(query, "") != ""),
                "same_coarse_neighborhood": int(neighborhood_of.get(partner, "") == q_neighborhood),
                "embedding_cosine_distance": cos_dist,
                "method": model_name,
                "sampling_bin": bin_name,
            })

        for g in positives:
            emit(g, forest_dist(query, g), "positive_same_neighborhood")
        for g, d in sampled_neg:
            emit(g, d, "negative_same_split")
        for g, d in binned:
            emit(g, d, "negative_distance_binned")
        n_unresolved += neg_missing
        for g in far_selected:
            emit(g, far_resolved[g], "negative_far")
    if not pair_rows:
        raise ValueError("No distance pairs sampled; check embedding coverage")
    strata = defaultdict(list)
    for row in pair_rows:
        level = quality_bin(manifests[row["genome_id_a"]].get("assembly_level", ""))
        stratum = (row["same_coarse_neighborhood"], row["same_near_clone_group"], level)
        strata[stratum].append(row)
    correlation_rows = []
    for (same_nb, same_clone, level), rows in sorted(strata.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
        rho = spearman_rho([r["embedding_cosine_distance"] for r in rows], [r["mash_distance"] for r in rows])
        correlation_rows.append({
            "model_name": model_name, "split": split,
            "stratum": f"same_nb={same_nb};same_clone={same_clone};{level}",
            "n_pairs": len(rows),
            "spearman_rho": rho if rho is not None else "",
            "method": "spearman(cosine_distance, mash_distance); no partial Mantel",
        })
    rho_all = spearman_rho([r["embedding_cosine_distance"] for r in pair_rows], [r["mash_distance"] for r in pair_rows])
    correlation_rows.append({
        "model_name": model_name, "split": split, "stratum": "ALL",
        "n_pairs": len(pair_rows),
        "spearman_rho": rho_all if rho_all is not None else "",
        "method": "spearman(cosine_distance, mash_distance); no partial Mantel",
    })
    bins = defaultdict(int)
    for row in pair_rows:
        bins[row["sampling_bin"]] += 1
    contract = {
        "schema": 1,
        "purpose": "Embedding vs Mash distance consistency on fixed seeded samples",
        "model": model,
        "model_name": model_name,
        "split": split,
        "seed": seed,
        "per_bin": per_bin,
        "bin_rules": {
            "positive_same_neighborhood": "train gallery members of the query sequence_neighborhood (stored forest edges)",
            "negative_same_split": "same split, different sequence_neighborhood (stored forest edges only; non-adjacent pairs dropped and counted)",
            "negative_distance_binned": "stored-forest distance within +/-0.005 of the query's closest resolved positive distance",
            "negative_far": "Mash distance >= 0.10, resolved via mash CLI over a bounded per-query candidate pool",
        },
        "correlation": "Spearman on sampled pairs; permutation validity not claimed",
        "max_far_candidates_per_query": max_far_candidates_per_query,
        "embeddings_dir": str(embeddings_dir),
        "split_manifest_sha256": sha256(split_manifest),
        "sequence_groups_sha256": sha256(groups_path),
        "forest_sha256": sha256(forest_path),
        "genomes_msh_sha256": sha256(sketch),
        "assemblies_sha256": sha256(assemblies),
        "mash_binary": str(mash_binary),
        "n_mash_cli_batches": n_cli_batches,
        "n_unresolved_pairs_dropped": n_unresolved,
        "n_embedding_failures": len(failures),
        "code_sha256": sha256(__file__),
    }
    prepare_output(output, contract)
    write_tsv(output / "distance_correlation.tsv", correlation_rows, CORRELATION_COLUMNS)
    _write_pairs(output, pair_rows)
    summary = {
        "status": "distance_eval_complete",
        "model": model, "model_name": model_name, "split": split,
        "n_queries": len(query_ids),
        "n_pairs": len(pair_rows),
        "bin_counts": dict(sorted(bins.items())),
        "spearman_all": rho_all if rho_all is not None else "",
        "n_mash_cli_batches": n_cli_batches,
        "n_unresolved_pairs_dropped": n_unresolved,
        "n_embedding_failures": len(failures),
        "seed": seed,
    }
    atomic_json(output / "summary.json", summary)
    return summary


def _write_pairs(output, rows):
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
        table = pa.table({c: [r[c] for r in rows] for c in PAIR_COLUMNS})
        pq.write_table(table, output / "pair_distance_sample.parquet")
        return "parquet"
    except ImportError:
        write_tsv(output / "pair_distance_sample.tsv", rows, PAIR_COLUMNS)
        return "tsv"
