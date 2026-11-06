"""Score operon boundary recovery on a locally annotated genome against author labels via an explicit mapping."""

import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output


def _official_functions():
    try:
        from bacbench.tasks.operon.run_evaluation_long_read_rna_seq import (
            operon_prot_indices_to_pairwise_labels,
            predict_pairwise_operon_boundaries,
        )
    except ImportError as error:
        raise RuntimeError("BacBench pinned checkout must be on PYTHONPATH (bacbench package)") from error
    return operon_prot_indices_to_pairwise_labels, predict_pairwise_operon_boundaries


def _load_embedding_matrix(parquet_path):
    import numpy as np
    import pandas as pd
    df = pd.read_parquet(parquet_path)
    if len(df) != 1:
        raise ValueError("Expected a single-genome embeddings parquet")
    row = df.iloc[0]
    if len(row["contig_name"]) != 1:
        raise ValueError("Expected a single contig")
    embedding = np.stack(row["embeddings"][0]).astype(np.float32)
    return row, embedding


def run_operon_local_eval(embeddings_parquet, mapping_tsv, reference_parquet, output,
                          method_name="unknown", model_revision="unknown",
                          official_embeddings_parquet=None):
    import numpy as np
    import pandas as pd
    embeddings_parquet, mapping_tsv, reference_parquet, output = map(Path, (embeddings_parquet, mapping_tsv, reference_parquet, output))
    if output.exists() or "raw" in output.resolve().parts:
        raise ValueError("Use a new non-raw output directory")
    labels_fn, scores_fn = _official_functions()
    emb_row, embedding = _load_embedding_matrix(embeddings_parquet)
    n_local = embedding.shape[0]
    strand = np.asarray(emb_row["strand"][0], dtype=int)
    if len(strand) != n_local:
        raise ValueError("Local strand vector length differs from embeddings")
    mapping = pd.read_csv(mapping_tsv, sep="\t", dtype={"local_index": "Int64"})
    author_to_local = {int(r.author_index): (None if pd.isna(r.local_index) else int(r.local_index)) for r in mapping.itertuples()}
    local_to_author = {v: k for k, v in author_to_local.items() if v is not None}
    if len(local_to_author) != sum(v is not None for v in author_to_local.values()):
        raise ValueError("Duplicate local indices in mapping")
    unmapped_authors = sorted(a for a, l in author_to_local.items() if l is None)
    ref = pd.read_parquet(reference_parquet).iloc[0]
    n_author = sum(len(x) for x in ref["protein_sequence"])
    author_labels = labels_fn(ref["operon_prot_indices"][0], n_genes=n_author)
    scores = scores_fn(emb=embedding, strand=strand)
    if len(scores) != n_local - 1:
        raise ValueError("Unexpected pairwise score length")
    pairs, excluded = [], Counter()
    for l in range(n_local - 1):
        a1, a2 = local_to_author.get(l), local_to_author.get(l + 1)
        if a1 is None or a2 is None:
            excluded["local_extra_protein_pair"] += 1
            continue
        if a2 != a1 + 1:
            between = range(a1 + 1, a2) if a2 > a1 else range(a2 + 1, a1)
            if all(a in unmapped_authors for a in between):
                excluded["author_unmapped_between"] += 1
            else:
                excluded["other_non_consecutive"] += 1
            continue
        pairs.append({"local_index_a": l, "local_index_b": l + 1, "author_index_a": a1,
                      "author_index_b": a2, "label": int(author_labels[a1]), "score": float(scores[l]),
                      "strand_match": bool(strand[l] == strand[l + 1])})
    from sklearn.metrics import average_precision_score, roc_auc_score
    labels = np.array([p["label"] for p in pairs])
    pair_scores = np.array([p["score"] for p in pairs])
    metrics = {
        "method": method_name, "model_revision": model_revision,
        "AUROC": float(roc_auc_score(labels, pair_scores)) if labels.min() != labels.max() else None,
        "average_precision": float(average_precision_score(labels, pair_scores)) if labels.min() != labels.max() else None,
        "n_aligned_pairs": len(pairs), "n_positive_pairs": int(labels.sum()),
        "n_negative_pairs": int(len(labels) - labels.sum()),
        "n_local_proteins": n_local, "n_author_proteins": int(n_author),
        "aligned_pairs_of_author_pairs": len(pairs) / (n_author - 1),
        "excluded_pairs": dict(excluded),
        "unmapped_author_indices": unmapped_authors,
        "n_author_pairs": int(n_author - 1),
    }
    aligned_author_pairs = [(p["author_index_a"], p["author_index_b"]) for p in pairs]
    official_same_pairs = None
    if official_embeddings_parquet is not None:
        official_embeddings_parquet = Path(official_embeddings_parquet)
        off_row, off_embedding = _load_embedding_matrix(official_embeddings_parquet)
        off_strand = np.asarray(off_row["strand"][0], dtype=int)
        if len(off_strand) != off_embedding.shape[0] or off_embedding.shape[0] != n_author:
            raise ValueError("Official embeddings do not match author protein count")
        official_scores = scores_fn(emb=off_embedding, strand=off_strand)
        idx = np.array([a for a, _ in aligned_author_pairs])
        official_same_pairs = {
            "AUROC": float(roc_auc_score(labels, official_scores[idx])),
            "average_precision": float(average_precision_score(labels, official_scores[idx])),
            "n_pairs": int(len(idx)),
            "note": "official-input scores projected onto the identical aligned author-pair set; isolates annotation drift from pipeline differences",
            "official_embeddings_sha256": sha256(official_embeddings_parquet),
        }
    contract = {
        "schema": 1,
        "purpose": "Operon boundary scoring on local annotation restricted to pairs aligned with the author protein list",
        "embeddings_sha256": sha256(embeddings_parquet),
        "mapping_sha256": sha256(mapping_tsv),
        "reference_sha256": sha256(reference_parquet),
        "code_sha256": sha256(__file__),
        "pair_definition": "local adjacent pair (l, l+1) included iff both map to consecutive author indices; label from official author pairwise labels; score from official BacBench function",
        "coordinate_note": "embeddings column is the official flat single-row nested list; reshaped to (n_proteins, dim) before scoring",
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    prepare_output(output, contract)
    write_tsv(output / "aligned_pairs.tsv", pairs, list(pairs[0]))
    metrics["official_same_pairs"] = official_same_pairs
    atomic_json(output / "metrics.json", metrics)
    summary = dict(metrics)
    summary.update(status="operon_local_eval_complete", finished_utc=datetime.now(timezone.utc).isoformat(),
                   limitations=[
                       "Only author-aligned adjacent pairs are scored; pairs touching locally extra or unmapped proteins are excluded and counted.",
                       "Exclusions are split into local extra-protein pairs and pairs spanning unmapped author proteins; both reduce mostly negative pairs, which shifts the base rate relative to the full author pair set.",
                       "Labels are defined on the author's ordered protein list; unmapped accessions may reflect annotation-version differences, which are enumerated rather than guessed.",
                       "Single-strain local-input acceptance check; not multi-strain evidence.",
                   ])
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in ("summary.json", "metrics.json", "aligned_pairs.tsv", "contract.json")})
    return summary
