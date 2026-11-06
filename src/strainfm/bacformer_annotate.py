"""Fixed protein-input annotation for frozen Bacformer/ESM embedding.

Companion GenBank (GBFF) annotations are preferred and parsed with the pinned
BacBench extractor; only genomes without usable companion annotation fall back
to pyrodigal on the nucleotide assembly. One row per protein, ordered by contig
order then CDS start. Work is chunked with per-chunk done files so an
interrupted run resumes without recomputing finished genomes.
"""

import gzip
import json
import platform
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from .audit import safe_path, sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output

PROTEIN_COLUMNS = [
    "genome_id", "assembly_accession", "contig_id", "contig_rank",
    "protein_index_in_contig", "start", "end", "strand", "protein_sequence",
    "aa_len", "annotation_tool", "annotation_version", "source_file_sha256",
]
SUMMARY_COLUMNS = [
    "genome_id", "split", "assembly_level", "status", "annotation_tool", "fallback_from_gbff", "source_file",
    "n_contigs", "n_cds", "n_cds_without_translation", "n_pseudo",
    "n_proteins", "n_proteins_over_1024aa", "n_over_max_proteins", "error",
]
FAILURE_COLUMNS = ["genome_id", "split", "error"]


def _find_gbff(annotation_root, accessions):
    if annotation_root is None:
        return None
    root = Path(annotation_root)
    candidates = set()
    for accession in dict.fromkeys(a for a in accessions if a and a != "missing"):
        for pattern in (f"{accession}*.gbff.gz", f"{accession}*.gbff",
                        f"{accession}/{accession}*.gbff.gz", f"{accession}/{accession}*.gbff",
                        f"*/{accession}*.gbff.gz", f"*/{accession}*.gbff"):
            candidates.update(root.glob(pattern))
    candidates = sorted(candidates)
    if len(candidates) > 1:
        raise ValueError(f"Multiple GBFF files for {accessions[0]}: {candidates}")
    return candidates[0] if candidates else None


def _gbff_coverage(path):
    from Bio import SeqIO
    opener = gzip.open if str(path).endswith(".gz") else open
    n_cds = n_no_translation = n_pseudo = 0
    contig_order = []
    with opener(path, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contig_order.append(record.id)
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                n_cds += 1
                if "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers:
                    n_pseudo += 1
                if "translation" not in feature.qualifiers:
                    n_no_translation += 1
    return {"n_cds": n_cds, "n_cds_without_translation": n_no_translation, "n_pseudo": n_pseudo, "contig_order": contig_order}


def _annotate_gbff(genome_id, gbff_path, max_proteins):
    try:
        from bacbench.pp.preprocess import extract_protein_info_from_gbff
    except ImportError as error:
        raise RuntimeError("GBFF mode requires the pinned BacBench checkout on PYTHONPATH") from error
    frame = extract_protein_info_from_gbff(str(gbff_path))
    coverage = _gbff_coverage(gbff_path)
    rank = {name: i for i, name in enumerate(dict.fromkeys(coverage["contig_order"]))}
    if frame.empty:
        return {"status": "no_cds", "rows": [], "coverage": coverage, "n_proteins": 0, "tool": "gbff"}
    frame = frame.copy()
    frame["contig_rank"] = frame["contig_name"].map(rank)
    if frame["contig_rank"].isna().any():
        raise ValueError(f"GBFF contig missing from record order: {genome_id}")
    frame = frame.sort_values(["contig_rank", "start", "end"], kind="stable").reset_index(drop=True)
    frame["protein_index_in_contig"] = frame.groupby("contig_rank").cumcount()
    rows = []
    for record in frame.itertuples():
        rows.append({
            "genome_id": genome_id,
            "assembly_accession": genome_id,
            "contig_id": record.contig_name,
            "contig_rank": int(record.contig_rank),
            "protein_index_in_contig": int(record.protein_index_in_contig),
            "start": int(record.start),
            "end": int(record.end),
            "strand": int(record.strand),
            "protein_sequence": record.protein_sequence,
            "aa_len": len(record.protein_sequence),
            "annotation_tool": "gbff",
            "annotation_version": "BacBench extract_protein_info_from_gbff@pinned",
            "source_file_sha256": sha256(gbff_path),
        })
    return {"status": "annotated" if rows else "no_cds", "rows": rows, "coverage": coverage,
            "n_proteins": len(rows), "tool": "gbff"}


def _annotate_prodigal(genome_id, fasta_path, min_aa_len, max_proteins):
    import pyrodigal
    from Bio import SeqIO
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    records = list(SeqIO.parse(opener(fasta_path, "rt"), "fasta"))
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train("N".join(str(record.seq) for record in records))
    version = f"pyrodigal={pyrodigal.__version__}"
    source_sha = sha256(fasta_path)
    rows, n_cds = [], 0
    for rank, record in enumerate(records):
        genes = sorted(finder.find_genes(str(record.seq)), key=lambda gene: (gene.begin, gene.end))
        index = 0
        for gene in genes:
            n_cds += 1
            translation = gene.translate()
            if not translation or len(translation) < min_aa_len:
                continue
            rows.append({
                "genome_id": genome_id,
                "assembly_accession": genome_id,
                "contig_id": record.id,
                "contig_rank": rank,
                "protein_index_in_contig": index,
                "start": int(gene.begin) - 1,
                "end": int(gene.end),
                "strand": 1 if gene.strand > 0 else -1,
                "protein_sequence": translation,
                "aa_len": len(translation),
                "annotation_tool": "pyrodigal",
                "annotation_version": version,
                "source_file_sha256": source_sha,
            })
            index += 1
    return {"status": "annotated" if rows else "no_cds", "rows": rows,
            "coverage": {"n_cds": n_cds, "n_cds_without_translation": 0, "n_pseudo": 0,
                         "contig_order": [record.id for record in records]},
            "n_proteins": len(rows), "tool": "pyrodigal"}


def _annotate_genome(genome_id, row, assemblies, genome_root, annotation_root, annotation_mode,
                     min_aa_len, max_proteins):
    accessions = [genome_id, row.get("file_accession", ""), row.get("paired_accession", "")]
    gbff = _find_gbff(annotation_root, accessions) if annotation_mode in ("auto", "gbff") else None
    if annotation_mode == "gbff" and gbff is None:
        raise FileNotFoundError(f"No GBFF found for {genome_id} under {annotation_root}")
    fasta = safe_path(genome_root, assemblies[genome_id]["target_relpath"])
    source = gbff if gbff is not None else fasta
    if gbff is None and not fasta.exists():
        raise FileNotFoundError(f"Missing FASTA and GBFF for {genome_id}")
    with tempfile.TemporaryDirectory() as staging:
        staged = Path(staging) / source.name
        shutil.copyfile(source, staged)
        if gbff is not None:
            result = _annotate_gbff(genome_id, staged, max_proteins)
            if result["status"] == "no_cds" and fasta.exists():
                fallback = Path(staging) / fasta.name
                shutil.copyfile(fasta, fallback)
                result = _annotate_prodigal(genome_id, fallback, min_aa_len, max_proteins)
                result["fallback_from_gbff"] = True
        else:
            result = _annotate_prodigal(genome_id, staged, min_aa_len, max_proteins)
    return result, source


def run_bacformer_annotate(split_manifest_path, assemblies_path, genome_root, output,
                           annotation_root=None, annotation_mode="auto", tool="prodigal",
                           translation_table=11, min_aa_len=30, max_proteins=6000,
                           splits=None, genome_ids=None, limit=None, chunk_size=64):
    split_manifest_path, assemblies_path, genome_root, output = map(Path, (split_manifest_path, assemblies_path, genome_root, output))
    if "raw" in output.resolve().parts:
        raise ValueError("Non-raw output directory required")
    if annotation_mode not in ("auto", "gbff", "prodigal"):
        raise ValueError("annotation_mode must be auto, gbff or prodigal")
    manifests = read_table(split_manifest_path)
    assemblies = {r["genome_id"]: r for r in read_table(assemblies_path)}
    if splits is not None:
        manifests = [r for r in manifests if r["split"] in set(splits)]
    if genome_ids is not None:
        wanted = set(genome_ids)
        manifests = [r for r in manifests if r["genome_id"] in wanted]
        missing_ids = wanted - {r["genome_id"] for r in manifests}
        if missing_ids:
            raise ValueError(f"Requested genomes absent from split manifest: {sorted(missing_ids)}")
    if limit is not None:
        manifests = sorted(manifests, key=lambda r: r["genome_id"])[:limit]
    if not manifests:
        raise ValueError("No genomes selected")
    missing = [gid for gid in {r["genome_id"] for r in manifests} if gid not in assemblies]
    if missing:
        raise ValueError(f"Genomes absent from assembly inventory: {missing[:5]}")
    contract = {
        "schema": 1,
        "purpose": "Fixed protein annotation of frozen split genomes for frozen Bacformer/ESM embedding",
        "annotation_mode": annotation_mode,
        "annotation_root": str(annotation_root) if annotation_root else None,
        "fallback_tool": tool,
        "translation_table": translation_table,
        "min_aa_len": min_aa_len,
        "max_proteins": max_proteins,
        "selected_splits": sorted(set(splits)) if splits is not None else ["train", "validation", "test"],
        "limit": limit, "chunk_size": chunk_size,
        "genome_root": str(genome_root),
        "split_manifest_sha256": sha256(split_manifest_path),
        "assemblies_sha256": sha256(assemblies_path),
        "code_sha256": sha256(__file__),
        "python": platform.python_version(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    identity_fields = tuple(key for key in contract if key not in ("started_utc",))
    if output.exists():
        contract_path = output / "contract.json"
        if not contract_path.exists():
            raise ValueError("Refusing nonempty output directory without contract")
        existing = json.loads(contract_path.read_text())
        mismatched = [field for field in identity_fields if field != "code_sha256" and existing.get(field) != contract.get(field)]
        if mismatched:
            raise ValueError(f"Resume contract changed: {mismatched}")
        if existing.get("code_sha256") != contract["code_sha256"]:
            existing.setdefault("code_sha256_history", []).append(existing["code_sha256"])
            existing["code_sha256"] = contract["code_sha256"]
            atomic_json(contract_path, existing)
        contract = existing
    else:
        prepare_output(output, contract)
    completion = output / "completion.json"
    if completion.exists():
        for name, digest in json.loads(completion.read_text()).items():
            if sha256(safe_path(output, name)) != digest:
                raise ValueError(f"Completed artifact changed: {name}")
        return json.loads((output / "summary.json").read_text())
    ordered = sorted(manifests, key=lambda r: r["genome_id"])
    for start in range(0, len(ordered), chunk_size):
        chunk = ordered[start:start + chunk_size]
        prefix = output / f"annotate_{start:06d}"
        done = prefix.with_suffix(".done.json")
        shard = prefix.with_suffix(".rows.jsonl")
        if not (done.exists() and shard.exists()):
            chunk_rows, chunk_summaries, chunk_failures = [], [], []
            for row in chunk:
                genome_id = row["genome_id"]
                try:
                    result, source = _annotate_genome(genome_id, row, assemblies, genome_root, annotation_root,
                                                      annotation_mode, min_aa_len, max_proteins)
                except Exception as error:
                    chunk_failures.append({"genome_id": genome_id, "split": row["split"], "error": f"{type(error).__name__}: {error}"})
                    continue
                genome_rows = result["rows"]
                coverage = result["coverage"]
                chunk_rows.extend(genome_rows)
                chunk_summaries.append({
                    "genome_id": genome_id, "split": row["split"],
                    "assembly_level": row.get("assembly_level", ""), "status": result["status"],
                    "annotation_tool": result["tool"], "fallback_from_gbff": bool(result.get("fallback_from_gbff")),
                    "source_file": str(source),
                    "n_contigs": len(coverage.get("contig_order", [])) or len({r["contig_id"] for r in genome_rows}),
                    "n_cds": coverage["n_cds"], "n_cds_without_translation": coverage["n_cds_without_translation"],
                    "n_pseudo": coverage["n_pseudo"], "n_proteins": len(genome_rows),
                    "n_proteins_over_1024aa": sum(r["aa_len"] > 1024 for r in genome_rows),
                    "n_over_max_proteins": len(genome_rows) > max_proteins, "error": "",
                })
            chunk_rows.sort(key=lambda r: (r["genome_id"], r["contig_rank"], r["start"], r["end"]))
            with shard.open("w") as handle:
                for protein_row in chunk_rows:
                    handle.write(json.dumps(protein_row, ensure_ascii=False) + "\n")
            atomic_json(done, {"genomes": len(chunk), "rows": len(chunk_rows),
                               "summaries": chunk_summaries, "failures": chunk_failures})
        completed = min(start + chunk_size, len(ordered))
        atomic_json(output / "progress.json", {"status": "annotating", "completed": completed, "total": len(ordered)})
        print(f"annotate {completed}/{len(ordered)}", flush=True)
    return _finalize_annotate(output, len(ordered), contract)


def _finalize_annotate(output, n_requested, contract):
    import pyarrow as pa
    import pyarrow.parquet as pq
    rows, summaries, failures = [], [], []
    for shard in sorted(output.glob("annotate_*.rows.jsonl")):
        with shard.open() as handle:
            for line in handle:
                if line:
                    rows.append(json.loads(line))
    for done in sorted(output.glob("annotate_*.done.json")):
        payload = json.loads(done.read_text())
        summaries.extend(payload.get("summaries", []))
        failures.extend(payload.get("failures", []))
    rows.sort(key=lambda r: (r["genome_id"], r["contig_rank"], r["start"], r["end"]))
    summaries.sort(key=lambda r: r["genome_id"])
    failures.sort(key=lambda r: r["genome_id"])
    table = pa.table({column: [row[column] for row in rows] for column in PROTEIN_COLUMNS})
    pq.write_table(table, output / "protein_manifest.parquet")
    write_tsv(output / "genome_input_summary.tsv", summaries, SUMMARY_COLUMNS)
    write_tsv(output / "annotation_failures.tsv", failures, FAILURE_COLUMNS)
    annotated = [s for s in summaries if s["status"] == "annotated"]
    fallbacks = sum(bool(s.get("fallback_from_gbff")) for s in summaries)
    summary = {
        "status": "annotation_complete" if not failures else "annotation_complete_with_failures",
        "n_genomes_requested": n_requested,
        "n_genomes_annotated": len(annotated),
        "n_failed": len(failures),
        "n_no_cds": sum(s["status"] == "no_cds" for s in summaries),
        "annotation_tool_counts": {tool: sum(s["annotation_tool"] == tool for s in summaries) for tool in ("gbff", "pyrodigal")},
        "n_gbff_fallbacks": fallbacks,
        "split_counts": {s: sum(r["split"] == s for r in summaries) for s in ("train", "validation", "test")},
        "n_proteins_total": sum(s["n_proteins"] for s in summaries),
        "n_proteins_median_per_genome": sorted(s["n_proteins"] for s in summaries)[len(summaries) // 2] if summaries else 0,
        "n_genomes_over_max_proteins": sum(s["n_over_max_proteins"] for s in summaries),
        "n_proteins_over_1024aa": sum(s["n_proteins_over_1024aa"] for s in summaries),
        "n_cds_without_translation": sum(s["n_cds_without_translation"] for s in summaries),
        "protein_manifest_sha256": sha256(output / "protein_manifest.parquet"),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "limitations": [
            "GBFF mode reuses NCBI CDS translations and never re-annotates genomes that already carry companion annotations.",
            "A GBFF without CDS features (WGS master-style record) falls back to pyrodigal and is counted in n_gbff_fallbacks.",
            "Prodigal fallback is an operational gene call, not a curated gene model.",
            "Protein order follows contig order then CDS start; circular genomes are not rotated.",
            "Max protein count enforcement happens at embedding time; over-limit genomes are counted here but not truncated here.",
        ],
    }
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in ("summary.json", "protein_manifest.parquet", "genome_input_summary.tsv", "annotation_failures.tsv", "contract.json")})
    (output / "README.md").write_text(
        "# Bacformer protein inputs\n\n"
        "One row per protein (ordered by genome, contig rank, CDS start) consumed by `bacformer-embed`. "
        "GBFF companion annotations are preferred; pyrodigal is a fallback for genomes without usable annotation.\n\n"
        f"- Requested: {summary['n_genomes_requested']}; annotated: {summary['n_genomes_annotated']}; failed: {summary['n_failed']}; no CDS: {summary['n_no_cds']}\n"
        f"- Tool counts: {summary['annotation_tool_counts']}; GBFF fallbacks: {summary['n_gbff_fallbacks']}\n"
        f"- Proteins: {summary['n_proteins_total']} total; genomes over {contract['max_proteins']}: {summary['n_genomes_over_max_proteins']}\n"
        f"- Split manifest SHA-256: {contract['split_manifest_sha256']}\n"
        f"- Protein manifest SHA-256: {summary['protein_manifest_sha256']}\n\n"
        "Consumers must truncate in contig order at the fixed `max_n_proteins` and apply the identical "
        "protein set to both ESM-C and Bacformer; do not reorder or re-filter proteins after this stage.\n"
    )
    return summary
