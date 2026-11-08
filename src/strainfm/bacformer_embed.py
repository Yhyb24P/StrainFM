"""Frozen ESM-C / Bacformer-L genome embeddings over one identical protein set.

Thin adapter around the pinned BacBench implementation: the ESM embeddings are
computed once per genome and reused as Bacformer inputs, so both genome vectors
(z = mean over the first `max_n_proteins` proteins in contig order) are defined
on exactly the same proteins.
"""

import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output

GENOME_COLUMNS = ["genome_id", "split", "model_name", "model_revision", "n_proteins", "truncated_protein_count", "embedding"]
ESM_COLUMNS = ["genome_id", "split", "model_name", "model_revision", "n_proteins", "embedding"]
FAILURE_COLUMNS = ["genome_id", "split", "stage", "error"]


def _genome_ids(dataset, split_by_genome):
    seen = set()
    for batch in dataset.to_batches(columns=["genome_id"]):
        seen.update(batch.column("genome_id").to_pylist())
    outside = sorted(gid for gid in seen if gid not in split_by_genome)
    if outside:
        raise ValueError(f"Protein manifest genome outside split manifest: {outside[:5]}")
    return sorted(seen)


def _load_chunk(dataset, genome_ids, split_by_genome):
    import pyarrow.dataset as ds
    table = dataset.to_table(
        columns=["genome_id", "contig_id", "contig_rank", "protein_sequence"],
        filter=ds.field("genome_id").isin(list(genome_ids)),
    )
    groups = {}
    for row in table.to_pylist():
        genome_id = row["genome_id"]
        groups.setdefault(genome_id, []).append((int(row["contig_rank"]), row["contig_id"], row["protein_sequence"]))
    for items in groups.values():
        items.sort(key=lambda item: item[0])
    return groups


def _nested_by_contig(items):
    contig_ids, nested = [], []
    for _, contig_id, sequence in items:
        if not contig_ids or contig_ids[-1] != contig_id:
            contig_ids.append(contig_id)
            nested.append([])
        nested[-1].append(sequence)
    return contig_ids, nested


def _flatten_shared(esm_nested, max_n_proteins):
    flat = [embedding for contig in esm_nested for embedding in contig]
    truncated = max(0, len(flat) - max_n_proteins)
    return flat[:max_n_proteins], truncated


class _Models:
    def __init__(self, bacformer_path, esm_path, bacformer_device, esm_device, dtype, esm_dtype="bfloat16"):
        import torch
        from transformers import AutoModel
        from bacbench.modeling.embedder import load_seq_embedder
        self.torch = torch
        self.device = bacformer_device
        self.esm = load_seq_embedder(str(esm_path), device=esm_device)
        if esm_dtype != "bfloat16":
            self.esm.model.to(esm_device, dtype=getattr(torch, esm_dtype))
            self.esm.dtype = getattr(torch, esm_dtype)
        self.bacformer = AutoModel.from_pretrained(str(bacformer_path), trust_remote_code=True).eval()
        self.bacformer.to(getattr(torch, dtype)).to(bacformer_device)

    def genome_vectors(self, contig_ids, nested_sequences, batch_size, max_seq_len, max_n_proteins, max_n_contigs):
        import numpy as np
        from bacbench.modeling.embed_prot_seqs import compute_bacformer_embeddings, compute_genome_protein_embeddings
        if len(contig_ids) > max_n_contigs:
            raise ValueError(f"Genome has {len(contig_ids)} contigs (> {max_n_contigs}); Bacformer would truncate contigs and break the shared protein set")
        esm_nested = compute_genome_protein_embeddings(
            embedder=self.esm, protein_sequences=nested_sequences, contig_ids=contig_ids,
            batch_size=batch_size, max_prot_seq_len=max_seq_len, genome_pooling_method=None,
        )
        flat, truncated = _flatten_shared(esm_nested, max_n_proteins)
        esm_vector = np.mean(np.stack(flat), axis=0).astype(np.float32)
        bacformer_nested = compute_bacformer_embeddings(
            model=self.bacformer, protein_embeddings=esm_nested, contig_ids=contig_ids,
            max_n_proteins=max_n_proteins, max_n_contigs=max_n_contigs,
            bacformer_model_type="large", genome_pooling_method=None,
        )
        bacformer_flat = [embedding for contig in bacformer_nested for embedding in contig if embedding is not None]
        if len(bacformer_flat) != len(flat):
            raise ValueError(f"Bacformer returned {len(bacformer_flat)} protein vectors for {len(flat)} shared proteins")
        bacformer_vector = np.mean(np.stack(bacformer_flat), axis=0).astype(np.float32)
        if not np.isfinite(esm_vector).all() or not np.isfinite(bacformer_vector).all():
            raise ValueError("Non-finite genome embedding")
        return esm_vector, bacformer_vector, len(flat), truncated


def _cuda_free_mib(device):
    if not device.startswith("cuda"):
        return None
    import torch
    if not torch.cuda.is_available():
        return 0
    index = int(device.split(":")[1]) if ":" in device else 0
    return torch.cuda.mem_get_info(index)[0] // (1024 * 1024)


def _bacbench_commit():
    try:
        import bacbench
        root = Path(bacbench.__file__).resolve().parents[1]
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    except Exception:
        return "unknown"


def run_bacformer_embed(protein_manifest_path, split_manifest_path, bacformer_path, esm_path, output,
                        batch_size=8, max_seq_len=1024, max_n_proteins=6000, max_n_contigs=1000,
                        device="cuda:0", dtype="bfloat16", min_free_mib=4000, chunk_size=16,
                        bacformer_revision=None, esm_revision=None,
                        esm_device="auto", esm_dtype="auto", esm_min_free_mib=1500):
    protein_manifest_path, split_manifest_path, bacformer_path, esm_path, output = map(
        Path, (protein_manifest_path, split_manifest_path, bacformer_path, esm_path, output))
    if "raw" in output.resolve().parts:
        raise ValueError("Non-raw output directory required")
    manifests = {r["genome_id"]: r for r in read_table(split_manifest_path)}
    import pyarrow.dataset as pyarrow_dataset
    dataset = pyarrow_dataset.dataset(str(protein_manifest_path), format="parquet")
    genome_ids = _genome_ids(dataset, manifests)
    if not genome_ids:
        raise ValueError("No genomes in protein manifest")
    import torch
    if device.startswith("cuda"):
        free = _cuda_free_mib(device)
        if free is not None and free < min_free_mib:
            raise RuntimeError(f"GPU {device} has {free} MiB free (< {min_free_mib} required); wait or rerun with --device cpu --dtype float32")
    if esm_device == "auto":
        free = _cuda_free_mib("cuda:0") if torch.cuda.is_available() else 0
        esm_device = "cuda:0" if free is not None and free >= esm_min_free_mib else "cpu"
    elif esm_device.startswith("cuda"):
        free = _cuda_free_mib(esm_device)
        if free is not None and free < esm_min_free_mib:
            raise RuntimeError(f"GPU {esm_device} has {free} MiB free for the ESM embedder (< {esm_min_free_mib} required)")
    if esm_dtype == "auto":
        esm_dtype = "bfloat16"
    contract = {
        "schema": 1,
        "purpose": "Frozen genome embeddings (ESM-C mean and Bacformer-L contextual mean) on one identical protein set",
        "bacformer_path": str(bacformer_path), "bacformer_revision": bacformer_revision,
        "esm_path": str(esm_path), "esm_revision": esm_revision,
        "max_n_proteins": max_n_proteins, "max_n_contigs": max_n_contigs,
        "max_seq_len": max_seq_len, "batch_size": batch_size,
        "truncation": "first max_n_proteins proteins in contig order, identical for both models",
        "pooling": "mean over per-protein vectors (contextual for Bacformer, mean-token for ESM)",
        "device": device, "esm_device": esm_device, "dtype": dtype, "esm_dtype": esm_dtype, "chunk_size": chunk_size,
        "protein_manifest_sha256": sha256(protein_manifest_path),
        "split_manifest_sha256": sha256(split_manifest_path),
        "bacbench_commit": _bacbench_commit(),
        "torch_version": torch.__version__, "cuda_version": torch.version.cuda,
        "python": platform.python_version(),
        "code_sha256": sha256(__file__),
        "training": False,
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    identity_fields = ("schema", "purpose", "bacformer_path", "bacformer_revision", "esm_path", "esm_revision",
                       "max_n_proteins", "max_n_contigs", "max_seq_len", "batch_size", "truncation", "pooling",
                       "chunk_size", "protein_manifest_sha256", "split_manifest_sha256", "bacbench_commit",
                       "torch_version", "cuda_version", "python", "training")
    if output.exists():
        contract_path = output / "contract.json"
        if not contract_path.exists():
            raise ValueError("Refusing nonempty output directory without contract")
        existing = json.loads(contract_path.read_text())
        mismatched = [field for field in identity_fields if existing.get(field) != contract.get(field)]
        if mismatched:
            raise ValueError(f"Resume contract changed: {mismatched}")
        changes = {}
        if (existing.get("device"), existing.get("esm_device"), existing.get("dtype"), existing.get("esm_dtype")) != (device, esm_device, dtype, esm_dtype):
            changes["device"] = {"from": [existing.get("device"), existing.get("esm_device"), existing.get("dtype"), existing.get("esm_dtype")],
                                 "to": [device, esm_device, dtype, esm_dtype]}
            existing["device"], existing["esm_device"], existing["dtype"], existing["esm_dtype"] = device, esm_device, dtype, esm_dtype
        if existing.get("code_sha256") != contract["code_sha256"]:
            existing.setdefault("code_sha256_history", []).append(existing["code_sha256"])
            changes["code_sha256"] = {"from": existing["code_sha256"], "to": contract["code_sha256"]}
            existing["code_sha256"] = contract["code_sha256"]
            existing["resume_code_change"] = "code changed between runs; equivalence of embedding math requires manual verification"
        if changes:
            existing.setdefault("resume_runs", []).append({"started_utc": contract["started_utc"], **changes})
            atomic_json(contract_path, existing)
        contract = existing
    else:
        prepare_output(output, contract)
    completion = output / "completion.json"
    if completion.exists():
        for name, digest in json.loads(completion.read_text()).items():
            if sha256(output / name) != digest:
                raise ValueError(f"Completed artifact changed: {name}")
        return json.loads((output / "summary.json").read_text())
    models = _Models(bacformer_path, esm_path, device, esm_device, dtype, esm_dtype)
    for start in range(0, len(genome_ids), chunk_size):
        chunk = genome_ids[start:start + chunk_size]
        prefix = output / f"embed_{start:06d}"
        done = prefix.with_suffix(".done.json")
        shard = prefix.with_suffix(".rows.jsonl")
        if done.exists() and shard.exists():
            recorded = json.loads(done.read_text()).get("rows_sha256")
            if recorded and sha256(shard) != recorded:
                raise ValueError(f"Chunk rows changed since completion: {shard.name}")
        else:
            rows, chunk_failures = [], []
            groups = _load_chunk(dataset, chunk, manifests)
            for genome_id in chunk:
                split = manifests[genome_id]["split"]
                if genome_id not in groups:
                    chunk_failures.append({"genome_id": genome_id, "split": split, "stage": "missing_proteins",
                                           "error": "No protein rows for genome in manifest"})
                    continue
                contig_ids, nested = _nested_by_contig(groups[genome_id])
                try:
                    esm_vector, bacformer_vector, used, truncated = models.genome_vectors(
                        contig_ids, nested, batch_size, max_seq_len, max_n_proteins, max_n_contigs)
                except Exception as error:
                    chunk_failures.append({"genome_id": genome_id, "split": split, "stage": "embedding", "error": f"{type(error).__name__}: {error}"})
                    continue
                rows.append({"genome_id": genome_id, "split": split, "n_proteins": used,
                             "truncated_protein_count": truncated,
                             "bacformer": [float(x) for x in bacformer_vector],
                             "esm": [float(x) for x in esm_vector]})
            with shard.open("w") as handle:
                for row in rows:
                    handle.write(json.dumps(row) + "\n")
            atomic_json(done, {"genomes": len(chunk), "embedded": len(rows), "failures": chunk_failures,
                               "rows_sha256": sha256(shard), "device": device, "esm_device": esm_device,
                               "dtype": dtype, "esm_dtype": esm_dtype})
        atomic_json(output / "progress.json", {"status": "embedding", "completed": min(start + chunk_size, len(genome_ids)), "total": len(genome_ids)})
        print(f"embed {min(start + chunk_size, len(genome_ids))}/{len(genome_ids)}", flush=True)
    _finalize(output, len(genome_ids), contract, max_n_proteins)
    return json.loads((output / "summary.json").read_text())


def _finalize(output, n_genomes_input, contract, max_n_proteins):
    import pyarrow as pa
    import pyarrow.parquet as pq
    bacformer_rows, esm_rows, summaries = [], [], []
    failures = []
    for done in sorted(output.glob("embed_*.done.json")):
        failures.extend(json.loads(done.read_text()).get("failures", []))
    failures.sort(key=lambda row: row["genome_id"])
    for shard in sorted(output.glob("embed_*.rows.jsonl")):
        with shard.open() as handle:
            for line in handle:
                if not line:
                    continue
                row = json.loads(line)
                summaries.append({"genome_id": row["genome_id"], "split": row["split"],
                                  "n_proteins": row["n_proteins"],
                                  "truncated_protein_count": row["truncated_protein_count"]})
                bacformer_rows.append({
                    "genome_id": row["genome_id"], "split": row["split"],
                    "model_name": "bacformer-large-frozen",
                    "model_revision": contract["bacformer_revision"] or "recorded-in-contract",
                    "n_proteins": row["n_proteins"], "truncated_protein_count": row["truncated_protein_count"],
                    "embedding": row["bacformer"],
                })
                esm_rows.append({
                    "genome_id": row["genome_id"], "split": row["split"],
                    "model_name": "esmplusplus-small-mean",
                    "model_revision": contract["esm_revision"] or "recorded-in-contract",
                    "n_proteins": row["n_proteins"], "embedding": row["esm"],
                })
    pq.write_table(pa.table({c: [r.get(c) for r in bacformer_rows] for c in GENOME_COLUMNS}), output / "genome_embeddings.parquet")
    pq.write_table(pa.table({c: [r.get(c) for r in esm_rows] for c in ESM_COLUMNS}), output / "esm_mean_embeddings.parquet")
    write_tsv(output / "embedding_failures.tsv", failures, FAILURE_COLUMNS)
    summary = {
        "status": "embedding_complete" if not failures else "embedding_complete_with_failures",
        "n_genomes_input": n_genomes_input,
        "n_genomes_embedded": len(summaries),
        "n_failed": len(failures),
        "n_truncated_at_max_proteins": sum(s["truncated_protein_count"] > 0 for s in summaries),
        "max_n_proteins": max_n_proteins,
        "split_counts": {s: sum(r["split"] == s for r in summaries) for s in ("train", "validation", "test")},
        "device": contract["device"], "esm_device": contract["esm_device"], "dtype": contract["dtype"], "esm_dtype": contract["esm_dtype"],
        "genome_embeddings_sha256": sha256(output / "genome_embeddings.parquet"),
        "esm_mean_embeddings_sha256": sha256(output / "esm_mean_embeddings.parquet"),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "training_performed": False,
        "limitations": [
            "Frozen pretrained inference only; no finetuning or continued pretraining was performed.",
            "Both genome vectors are means over the same first-N proteins in contig order; truncated counts are recorded per genome.",
            "Bacformer pretraining exposure to these genomes remains unknown; results diagnose within-species structure recovery, not unseen-corpus generalisation.",
        ],
    }
    atomic_json(output / "summary.json", summary)
    contract["finished_utc"] = summary["finished_utc"]
    atomic_json(output / "contract.json", contract)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in ("summary.json", "genome_embeddings.parquet", "esm_mean_embeddings.parquet", "embedding_failures.tsv", "contract.json")})
    (output / "README.md").write_text(
        "# Frozen genome embeddings\n\n"
        "ESM-C per-protein vectors are computed once and reused as Bacformer inputs, so both genome-level "
        "vectors share exactly the same protein set (first max_n_proteins proteins in contig order).\n\n"
        f"- Bacformer snapshot: {contract['bacformer_path']}\n"
        f"- ESM snapshot: {contract['esm_path']}\n"
        f"- Device/dtype: {contract['device']} / {contract['dtype']}; torch {contract['torch_version']}\n"
        f"- Protein manifest SHA-256: {contract['protein_manifest_sha256']}\n"
        f"- Split manifest SHA-256: {contract['split_manifest_sha256']}\n"
    )
    return summary
