import json
import tempfile
import unittest
from pathlib import Path

from strainfm.audit import write_tsv
from strainfm.cohort import read_table, run_cohort
from strainfm.split import run_split


class CohortSplitTest(unittest.TestCase):
    def test_biosample_merges_and_species_mismatch_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            samples = [{"genome_id": name, "file_accession": name, "status": "valid", "target_relpath": name, "total_bases": "1000", "ambiguous_fraction": "0", "contig_multiset_sha256": name} for name in "abc"]
            metadata = [{"file_accession": name, "report_status": "found", "organism_name": "Species a strain", "submitted_species": "Species a", "assembly_level": "Contig", "biosample": "same" if name in "ab" else "other", "bioproject": "project", "contig_count": "1", "contig_n50": "1000", "total_sequence_length": "1000"} for name in "abc"]
            metadata[-1]["submitted_species"] = "Species b"
            write_tsv(root / "samples.tsv", samples, list(samples[0]))
            write_tsv(root / "metadata.tsv", metadata, list(metadata[0]))
            config = {"scientific_name": "Species a", "genome_root": str(root), "min_bases": 100, "max_bases": 2000, "min_contig_n50": 100, "max_contigs": 10, "max_ambiguous_fraction": 0.01}
            (root / "config.json").write_text(json.dumps(config))
            summary = run_cohort(root / "samples.tsv", root / "metadata.tsv", root / "config.json", root / "output")
            self.assertEqual(summary["eligible_files"], 2)
            self.assertEqual(summary["identity_groups"], 1)
            self.assertEqual(summary["exclusion_reasons_nonexclusive"], {"metadata_species_mismatch": 1})

    def test_split_source_and_nearclone_isolation_and_query_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            samples, groups = [], []
            for i in range(20):
                gid = f"g{i:02d}"
                samples.append({"genome_id": gid, "file_accession": gid, "biosample": gid, "identity_group": gid, "bioproject": f"project{i // 2}", "assembly_level": "Contig", "contig_count": "1", "contig_n50": "1000"})
                groups.append({"genome_id": gid, "near_clone_group": "clone" if i < 2 else gid, "sequence_neighborhood": "common"})
            write_tsv(root / "samples.tsv", samples, list(samples[0]))
            write_tsv(root / "groups.tsv", groups, list(groups[0]))
            members = [{**row, "eligible": "True", "identity_representative": row["genome_id"]} for row in samples]
            write_tsv(root / "members.tsv", members, list(members[0]))
            summary = run_split(root / "samples.tsv", root / "groups.tsv", root / "members.tsv", root / "output")
            self.assertEqual(summary["near_clone_representatives"], 19)
            self.assertTrue(summary["verified_no_cross_split_identity_nearclone_or_known_project"])
            self.assertGreater(summary["retrieval_positive_supported_queries"]["test"], 0)
            emitted = read_table(root / "output/split_manifest.tsv")
            owners = {}
            for row in emitted:
                owners.setdefault(row["bioproject"], set()).add(row["split"])
            self.assertTrue(all(len(splits) == 1 for splits in owners.values()))

    def test_nonrepresentative_identity_and_clone_sources_remain_isolated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            samples, groups = [], []
            for i in range(24):
                gid = f"g{i:02d}"
                samples.append({"genome_id": gid, "file_accession": gid, "biosample": gid, "identity_group": gid, "bioproject": f"P{i}", "assembly_level": "Contig", "contig_count": "1", "contig_n50": "1000"})
                groups.append({"genome_id": gid, "near_clone_group": "clone" if i < 2 else gid, "sequence_neighborhood": "common"})
            members = [{**row, "eligible": "True", "identity_representative": row["genome_id"]} for row in samples]
            # Removed identity member links P0 to P2. Removed clone member g01
            # links that same component to P1; both relations must survive.
            members.append({**members[0], "genome_id": "duplicate", "bioproject": "P2"})
            for name, table in [("samples", samples), ("groups", groups), ("members", members)]:
                write_tsv(root / f"{name}.tsv", table, list(table[0]))
            result = run_split(root / "samples.tsv", root / "groups.tsv", root / "members.tsv", root / "output")
            emitted = {r["genome_id"]: r for r in read_table(root / "output/split_manifest.tsv")}
            self.assertNotIn("g01", emitted)
            self.assertEqual(set(json.loads(emitted["g00"]["source_projects"])), {"P0", "P1", "P2"})
            self.assertEqual(emitted["g00"]["source_group"], emitted["g02"]["source_group"])
            self.assertEqual(emitted["g00"]["split"], emitted["g02"]["split"])
            self.assertTrue(result["full_member_source_constraints_verified"])
