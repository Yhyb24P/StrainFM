"""Build a measured-AST cohort from verified local and ATB genome files."""

import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .audit import safe_path, sha256, write_tsv
from .cohort import UnionFind, rank_quality, read_table
from .metadata import known
from .snapshot import atomic_json, inspect_file


def run_amr_cohort(samples_path, local_cohort_path, config_path, output, workers=8):
    samples_path, local_cohort_path, config_path, output = map(
        Path, (samples_path, local_cohort_path, config_path, output))
    if output.exists() or "raw" in output.resolve().parts:
        raise ValueError("Use a new non-raw output directory")
    if not 1 <= workers <= 24:
        raise ValueError("workers must be between 1 and 24")
    config = json.loads(config_path.read_text())
    root = Path(config["genome_root"]).resolve()
    atb_root = Path(config["atb_genome_root"])
    samples = read_table(samples_path)
    if len({r["genome_id"] for r in samples}) != len(samples):
        raise ValueError("Duplicate sample ID")
    local = {r["genome_id"]: r for r in read_table(local_cohort_path)}
    selected = [r for r in samples if r["action"] in ("reuse_local", "download_atb")]
    if not selected:
        raise ValueError("No candidate genomes")
    items = []
    for row in selected:
        row = dict(row)
        if row["genome_id"] != row["biosample"]:
            raise ValueError("AST cohort requires one genome per BioSample")
        if row.get("ena_taxonomy_conflict") == "True":
            raise ValueError(f"Unresolved taxonomy conflict: {row['genome_id']}")
        old = local.get(row.get("local_genome_id"), {}) if row["action"] == "reuse_local" else {}
        if row["action"] == "reuse_local":
            if (not old or old["biosample"] != row["biosample"] or old["eligible"] != "True"
                    or old["file_sha256"] != row["expected_local_sha256"]):
                raise ValueError(f"Local identity evidence changed: {row['genome_id']}")
            path = Path(row["local_fasta_path"]).resolve()
            row["assembly_level"] = old.get("assembly_level", "")
        else:
            if row.get("atb_species") != "Escherichia coli" or row.get("atb_exclusion_reason"):
                raise ValueError(f"ATB metadata not eligible: {row['genome_id']}")
            path = safe_path(atb_root, row["biosample"] + ".fa.gz")
            row["file_accession"] = row["atb_assembly_accession"]
            row["assembly_level"] = ""
        relative = str(path.relative_to(root))
        safe_path(root, relative)
        # Fail before making a partial cohort if the download is incomplete.
        stat = path.stat()
        sources = set(row.get("ena_study_accessions", "").split(";"))
        sources.update(old.get("bioproject", "").split(";"))
        sources = sorted(s for s in sources if known(s))
        row.update(target_relpath=relative, initial_size=stat.st_size,
                   initial_mtime_ns=stat.st_mtime_ns, species_id="Escherichia_coli",
                   source_projects=json.dumps(sources, separators=(",", ":")),
                   bioproject=old.get("bioproject", ""))
        items.append((row, path))
    if workers == 1:
        records = [inspect_file(item) for item in items]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            records = list(pool.map(inspect_file, items))
    for row in records:
        if row["status"] != "valid":
            raise ValueError(f"Invalid genome {row['genome_id']}: {row['error']}")
        if row["action"] == "reuse_local" and row["file_sha256"] != row["expected_local_sha256"]:
            raise ValueError(f"Local genome hash changed: {row['genome_id']}")
        reasons = []
        if not config["min_bases"] <= row["total_bases"] <= config["max_bases"]:
            reasons.append("length_outside_range")
        if row["sequence_count"] > config["max_contigs"]:
            reasons.append("contig_count_above_maximum")
        if row["n50_bases"] < config["min_contig_n50"]:
            reasons.append("contig_n50_below_minimum")
        if row["ambiguous_fraction"] > config["max_ambiguous_fraction"]:
            reasons.append("ambiguous_fraction_above_maximum")
        row.update(contig_count=row["sequence_count"], contig_n50=row["n50_bases"],
                   eligible=not reasons, exclusion_reason=";".join(reasons))
    eligible = {r["genome_id"]: r for r in records if r["eligible"]}
    union, seen = UnionFind(eligible), {}
    for gid, row in eligible.items():
        digest = row["contig_multiset_sha256"]
        if digest in seen:
            union.union(seen[digest], gid)
        seen[digest] = gid
    chosen = {}
    for group, ids in union.groups().items():
        representative = min(ids, key=lambda gid: rank_quality(eligible[gid]))
        chosen.update({gid: (group, representative) for gid in ids})
    for row in records:
        group, representative = chosen.get(row["genome_id"], ("", ""))
        row.update(identity_group=group, identity_representative=representative,
                   retained_identity_representative=row["genome_id"] == representative)
    columns = list(dict.fromkeys(k for r in records for k in r))
    output.mkdir(parents=True)
    write_tsv(output / "cohort_manifest.tsv", records, columns)
    write_tsv(output / "identity_representatives.tsv",
              [r for r in records if r["retained_identity_representative"]], columns)
    summary = {
        "status": "quality_cohort_ready_for_sequence_grouping" if eligible else "no_genomes_pass_quality", "config": config,
        "source_actions": dict(Counter(r["action"] for r in samples)),
        "candidate_files": len(records), "eligible_files": len(eligible),
        "identity_groups": len(union.groups()),
        "exclusion_reasons_nonexclusive": dict(Counter(
            reason for r in records for reason in r["exclusion_reason"].split(";") if reason)),
        "input_sha256": {str(p): sha256(p) for p in (samples_path, local_cohort_path, config_path)},
        "code_sha256": {str(p): sha256(p) for p in (Path(__file__), Path(__file__).with_name("snapshot.py"), Path(__file__).with_name("audit.py"))},
        "cohort_sha256": sha256(output / "cohort_manifest.tsv"),
        "limits": ["BioSample and exact sequence groups are not independent strain counts.",
                   "All eligible BioSamples retain their own labels; representatives are for sequence grouping only.",
                   "Near-clone/source isolation and task splits are still required before model processing."]}
    atomic_json(output / "summary.json", summary)
    return summary
