import argparse
import json
from pathlib import Path

from .audit import run_audit
from .labels import run_label_audit
from .metadata import run_metadata_audit
from .snapshot import run_snapshot
from .metadata_bulk import run_bulk_metadata
from .cohort import run_cohort
from .diversity import run_diversity
from .split import run_split
from .download import run_download
from .virus import split_virus
from .amr_intersection import run_amr_intersection
from .operon_local_eval import run_operon_local_eval
from .bacformer_annotate import run_bacformer_annotate
from .bacformer_embed import run_bacformer_embed
from .eval_retrieval import run_eval_retrieval
from .eval_distance import run_eval_distance


def main():
    parser = argparse.ArgumentParser(description="StrainFM research CLI")
    sub = parser.add_subparsers(dest="command", required=True)
    audit = sub.add_parser("audit", help="Read-only manifest and sampled FASTA audit")
    audit.add_argument("--config", type=Path, required=True)
    audit.add_argument("--output", type=Path, required=True)
    labels = sub.add_parser("labels", help="Audit binary-label coverage by exact versioned accession")
    labels.add_argument("--labels", type=Path, required=True)
    labels.add_argument("--download-manifest", type=Path, required=True)
    labels.add_argument("--output", type=Path, required=True)
    metadata = sub.add_parser("metadata", help="Fetch NCBI assembly metadata for an existing audit sample")
    metadata.add_argument("--sample", type=Path, required=True)
    metadata.add_argument("--output", type=Path, required=True)
    snapshot = sub.add_parser("snapshot", help="Freeze and validate all mapped files for one species; resumable")
    snapshot.add_argument("--config", type=Path, required=True)
    snapshot.add_argument("--output", type=Path, required=True)
    snapshot.add_argument("--workers", type=int, default=8)
    bulk = sub.add_parser("metadata-bulk", help="Fetch versioned NCBI metadata for a species inventory; resumable")
    bulk.add_argument("--inventory", type=Path, required=True)
    bulk.add_argument("--output", type=Path, required=True)
    cohort = sub.add_parser("cohort", help="Filter and deduplicate a completed snapshot using exact identity evidence")
    cohort.add_argument("--snapshot", type=Path, required=True)
    cohort.add_argument("--metadata", type=Path, required=True)
    cohort.add_argument("--config", type=Path, required=True)
    cohort.add_argument("--output", type=Path, required=True)
    diversity = sub.add_parser("diversity", help="Sketch and group a species cohort with Mash; resumable sketch batches")
    diversity.add_argument("--representatives", type=Path, required=True)
    diversity.add_argument("--config", type=Path, required=True)
    diversity.add_argument("--output", type=Path, required=True)
    diversity.add_argument("--threads", type=int, default=8)
    split = sub.add_parser("split", help="Freeze candidate source-held-out train/validation/test after sequence grouping")
    split.add_argument("--representatives", type=Path, required=True)
    split.add_argument("--groups", type=Path, required=True)
    split.add_argument("--cohort-manifest", type=Path, required=True)
    split.add_argument("--output", type=Path, required=True)
    split.add_argument("--seed", type=int, default=20261105)
    download = sub.add_parser("download", help="Download/retry manifest FASTA files, validating existing and new files")
    download.add_argument("--manifest", type=Path, required=True)
    download.add_argument("--root", type=Path, required=True)
    virus = sub.add_parser("split-virus", help="Split UHGV into a new output directory after successful decompression")
    virus.add_argument("--manifest", type=Path, required=True)
    virus.add_argument("--input", type=Path, required=True)
    virus.add_argument("--output", type=Path, required=True)
    amr = sub.add_parser("amr-intersection", help="Intersect BacBench AMR binary labels with the frozen E. coli cohort")
    amr.add_argument("--labels", type=Path, required=True)
    amr.add_argument("--cohort-manifest", type=Path, required=True)
    amr.add_argument("--split-manifest", type=Path, required=True)
    amr.add_argument("--output", type=Path, required=True)
    local = sub.add_parser("operon-local-eval", help="Score operon boundaries on a local annotation via explicit author-list mapping")
    local.add_argument("--embeddings", type=Path, required=True)
    local.add_argument("--mapping", type=Path, required=True)
    local.add_argument("--reference", type=Path, required=True)
    local.add_argument("--output", type=Path, required=True)
    local.add_argument("--method-name", default="unknown")
    local.add_argument("--model-revision", default="unknown")
    local.add_argument("--official-embeddings", type=Path, default=None)
    annotate = sub.add_parser("bacformer-annotate", help="Build ordered protein inputs from companion GBFF or pyrodigal fallback")
    annotate.add_argument("--split-manifest", type=Path, required=True)
    annotate.add_argument("--assemblies", type=Path, required=True)
    annotate.add_argument("--genome-root", type=Path, default=Path("/mnt/e/StrainFM/data/raw/genomes"))
    annotate.add_argument("--annotation-root", type=Path, default=None)
    annotate.add_argument("--annotation-mode", choices=["auto", "gbff", "prodigal"], default="auto")
    annotate.add_argument("--tool", default="prodigal")
    annotate.add_argument("--translation-table", type=int, default=11)
    annotate.add_argument("--min-aa-len", type=int, default=30)
    annotate.add_argument("--max-proteins", type=int, default=6000)
    annotate.add_argument("--splits", default=None)
    annotate.add_argument("--genome-ids", default=None)
    annotate.add_argument("--genome-list", type=Path, default=None)
    annotate.add_argument("--limit", type=int, default=None)
    annotate.add_argument("--chunk-size", type=int, default=64)
    annotate.add_argument("--workers", type=int, default=1)
    annotate.add_argument("--output", type=Path, required=True)
    embed = sub.add_parser("bacformer-embed", help="Frozen ESM-C and Bacformer-L genome vectors on one identical protein set")
    embed.add_argument("--protein-manifest", type=Path, required=True)
    embed.add_argument("--split-manifest", type=Path, required=True)
    embed.add_argument("--bacformer-path", type=Path, required=True)
    embed.add_argument("--esm-path", type=Path, required=True)
    embed.add_argument("--bacformer-revision", default=None)
    embed.add_argument("--esm-revision", default=None)
    embed.add_argument("--batch-size", type=int, default=8)
    embed.add_argument("--max-seq-len", type=int, default=1024)
    embed.add_argument("--max-n-proteins", type=int, default=6000)
    embed.add_argument("--max-n-contigs", type=int, default=1000)
    embed.add_argument("--device", default="cuda:0")
    embed.add_argument("--esm-device", default="auto")
    embed.add_argument("--esm-dtype", default="auto")
    embed.add_argument("--esm-min-free-mib", type=int, default=1500)
    embed.add_argument("--dtype", default="bfloat16")
    embed.add_argument("--min-free-mib", type=int, default=4000)
    embed.add_argument("--chunk-size", type=int, default=16)
    embed.add_argument("--output", type=Path, required=True)
    retrieval = sub.add_parser("eval-retrieval", help="Frozen-protocol retrieval evaluation (protocol A train gallery / B held-out)")
    retrieval.add_argument("--embeddings-dir", type=Path, required=True)
    retrieval.add_argument("--split-manifest", type=Path, required=True)
    retrieval.add_argument("--support", type=Path, required=True)
    retrieval.add_argument("--groups", type=Path, required=True)
    retrieval.add_argument("--output", type=Path, required=True)
    retrieval.add_argument("--model", choices=["bacformer", "esm_mean"], default="bacformer")
    retrieval.add_argument("--protocol", choices=["A", "B"], default="A")
    retrieval.add_argument("--split", choices=["validation", "test"], default="validation")
    retrieval.add_argument("--k-values", default="1,5,10")
    retrieval.add_argument("--seed", type=int, default=20261106)
    retrieval.add_argument("--model-name", default="bacformer-large-frozen")
    distance = sub.add_parser("eval-distance", help="Embedding distance vs Mash distance consistency")
    distance.add_argument("--embeddings-dir", type=Path, required=True)
    distance.add_argument("--split-manifest", type=Path, required=True)
    distance.add_argument("--groups-dir", type=Path, required=True)
    distance.add_argument("--assemblies", type=Path, required=True)
    distance.add_argument("--genome-root", type=Path, default=Path("/mnt/e/StrainFM/data/raw/genomes"))
    distance.add_argument("--mash-binary", type=Path, required=True)
    distance.add_argument("--output", type=Path, required=True)
    distance.add_argument("--model", choices=["bacformer", "esm_mean"], default="bacformer")
    distance.add_argument("--split", choices=["validation", "test"], default="validation")
    distance.add_argument("--seed", type=int, default=20261106)
    distance.add_argument("--model-name", default="bacformer-large-frozen")
    distance.add_argument("--per-bin", type=int, default=20)
    distance.add_argument("--max-far-candidates-per-query", type=int, default=400)
    distance.add_argument("--cli-workers", type=int, default=8)
    args = parser.parse_args()
    if args.command == "audit":
        result = run_audit(args.config, args.output)
    elif args.command == "labels":
        result = run_label_audit(args.labels, args.download_manifest, args.output)
    elif args.command == "metadata":
        result = run_metadata_audit(args.sample, args.output)
    elif args.command == "snapshot":
        result = run_snapshot(args.config, args.output, args.workers)
    elif args.command == "metadata-bulk":
        result = run_bulk_metadata(args.inventory, args.output)
    elif args.command == "cohort":
        result = run_cohort(args.snapshot, args.metadata, args.config, args.output)
    elif args.command == "diversity":
        result = run_diversity(args.representatives, args.config, args.output, args.threads)
    elif args.command == "split":
        result = run_split(args.representatives, args.groups, args.cohort_manifest, args.output, args.seed)
    elif args.command == "download":
        result = run_download(args.manifest, args.root)
    elif args.command == "amr-intersection":
        result = run_amr_intersection(args.labels, args.cohort_manifest, args.split_manifest, args.output)
    elif args.command == "operon-local-eval":
        result = run_operon_local_eval(args.embeddings, args.mapping, args.reference, args.output, args.method_name, args.model_revision, args.official_embeddings)
    elif args.command == "bacformer-annotate":
        result = run_bacformer_annotate(
            args.split_manifest, args.assemblies, args.genome_root, args.output,
            annotation_root=args.annotation_root, annotation_mode=args.annotation_mode,
            tool=args.tool, translation_table=args.translation_table, min_aa_len=args.min_aa_len,
            max_proteins=args.max_proteins,
            splits=args.splits.split(",") if args.splits else None,
            genome_ids=(args.genome_ids.split(",") if args.genome_ids else
                        [line.strip() for line in args.genome_list.read_text().splitlines() if line.strip()] if args.genome_list else None),
            limit=args.limit, chunk_size=args.chunk_size, workers=args.workers)
    elif args.command == "bacformer-embed":
        result = run_bacformer_embed(
            args.protein_manifest, args.split_manifest, args.bacformer_path, args.esm_path, args.output,
            batch_size=args.batch_size, max_seq_len=args.max_seq_len, max_n_proteins=args.max_n_proteins,
            max_n_contigs=args.max_n_contigs, device=args.device, dtype=args.dtype,
            min_free_mib=args.min_free_mib, chunk_size=args.chunk_size,
            bacformer_revision=args.bacformer_revision, esm_revision=args.esm_revision,
            esm_device=args.esm_device, esm_dtype=args.esm_dtype, esm_min_free_mib=args.esm_min_free_mib)
    elif args.command == "eval-retrieval":
        result = run_eval_retrieval(
            args.embeddings_dir, args.split_manifest, args.support, args.groups, args.output,
            model=args.model, protocol=args.protocol, split=args.split,
            k_values=tuple(int(k) for k in args.k_values.split(",")), seed=args.seed,
            model_name=args.model_name)
    elif args.command == "eval-distance":
        result = run_eval_distance(
            args.embeddings_dir, args.split_manifest, args.groups_dir, args.assemblies,
            args.genome_root, args.mash_binary, args.output, model=args.model, split=args.split,
            seed=args.seed, model_name=args.model_name, per_bin=args.per_bin,
            max_far_candidates_per_query=args.max_far_candidates_per_query, cli_workers=args.cli_workers)
    else:
        result = split_virus(args.manifest, args.input, args.output)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if args.command == "download" and result["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
