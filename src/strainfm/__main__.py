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
    else:
        result = split_virus(args.manifest, args.input, args.output)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if args.command == "download" and result["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
