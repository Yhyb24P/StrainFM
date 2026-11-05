import argparse
import json
from pathlib import Path

from .audit import run_audit


def main():
    parser = argparse.ArgumentParser(description="StrainFM research CLI")
    sub = parser.add_subparsers(dest="command", required=True)
    audit = sub.add_parser("audit", help="Read-only manifest and sampled FASTA audit")
    audit.add_argument("--config", type=Path, required=True)
    audit.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run_audit(args.config, args.output)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
