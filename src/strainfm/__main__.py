import argparse
import json
from pathlib import Path

from .audit import run_audit
from .labels import run_label_audit
from .metadata import run_metadata_audit


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
    args = parser.parse_args()
    if args.command == "audit":
        result = run_audit(args.config, args.output)
    elif args.command == "labels":
        result = run_label_audit(args.labels, args.download_manifest, args.output)
    else:
        result = run_metadata_audit(args.sample, args.output)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
