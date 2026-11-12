import json
import tempfile
import unittest
from pathlib import Path

from strainfm.amr_sources import _aliases, _read_pages, apply_atb_metadata, apply_ena_runs, link_samples, phenotype_state
from strainfm.audit import sha256, write_tsv


class AmrSourcesTest(unittest.TestCase):
    def test_missing_call_does_not_vote_and_intermediate_is_not_susceptible(self):
        self.assertEqual(phenotype_state([{}, {"resistant_phenotype": "Resistant"}]), "R")
        self.assertEqual(phenotype_state([{}]), "unknown")
        self.assertEqual(phenotype_state([{"resistant_phenotype": "Intermediate"}]), "I")
        self.assertEqual(phenotype_state([{"resistant_phenotype": x} for x in ("Resistant", "Susceptible")]), "conflict")

    def test_alias_merge_and_disagreement_quarantine_do_not_transfer_labels(self):
        metadata = {
            "g1": {"biosample_accession": "ERS1"},
            "g2": {"biosample_accession": "SAMN1"},
            "g3": {"biosample_accession": "SAMN2", "assembly_accession": "GCA_1.1"},
            "g4": {"biosample_accession": "", "assembly_accession": "GCA_1.1"},
        }
        records = [{"id": str(i), "genome_id": gid, "antibiotic": "ciprofloxacin",
                    "resistant_phenotype": call}
                   for i, (gid, call) in enumerate((("g1", "Resistant"), ("g2", ""),
                                                    ("g3", "Susceptible"), ("g4", "Susceptible")))]
        local = {"genome_id": "local", "biosample": "SAMN1", "file_accession": "GCA_1.1",
                 "eligible": "True", "status": "valid", "target_relpath": "one.fna.gz",
                 "file_sha256": "digest", "contig_n50": "200000", "contig_count": "20"}
        units, samples = link_samples(records, metadata, {"ERS1": "SAMN1"}, [local])
        bysample = {r["biosample"]: r for r in samples}
        self.assertEqual(bysample["SAMN1"]["ciprofloxacin"], "R")
        self.assertEqual(bysample["SAMN1"]["action"], "reuse_local")
        self.assertEqual(bysample["SAMN2"]["action"], "quarantine_identity_conflict")
        self.assertEqual(bysample["SAMN2"]["local_genome_id"], "")
        self.assertEqual(len(samples), 2)
        self.assertEqual(sum(not r["biosample"] for r in units), 1)
        self.assertEqual(next(r for r in units if r["biosample"] == "SAMN1")["record_count"], 2)

    def test_response_tampering_and_path_escape_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "snapshot"
            root.mkdir()
            data = root / "page.json"
            data.write_text('[{"id":"one"}]')
            manifest = root / "manifest.json"
            page = {"path": str(data), "sha256": sha256(data), "rows": 1}
            manifest.write_text(json.dumps([page]))
            self.assertEqual(_read_pages(manifest, root, {}), [{"id": "one"}])
            data.write_text('[{"id":"changed"}]')
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                _read_pages(manifest, root, {})
            outside = Path(tmp) / "outside.json"
            outside.write_text("[]")
            page.update(path=str(outside), sha256=sha256(outside), rows=0)
            manifest.write_text(json.dumps([page]))
            with self.assertRaisesRegex(ValueError, "outside source"):
                _read_pages(manifest, root, {})

    def test_ena_mapping_uses_external_biosample_not_primary_ers_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "ena_aliases").mkdir()
            xml = root / "ena_aliases/batch_00.xml"
            xml.write_text('<SAMPLE_SET><SAMPLE accession="ERS1"><IDENTIFIERS>'
                           '<PRIMARY_ID>ERS1</PRIMARY_ID><EXTERNAL_ID namespace="BioSample">SAMN1</EXTERNAL_ID>'
                           '</IDENTIFIERS></SAMPLE></SAMPLE_SET>')
            xml.with_suffix(".request.json").write_text(json.dumps({"sha256": sha256(xml)}))
            self.assertEqual(_aliases(root, {}), {"ERS1": "SAMN1"})

    def test_remote_fasta_presence_does_not_override_failed_species_or_qc(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "assembly.tsv"
            rows = []
            for sample in ("SAMN1", "SAMN2", "SAMN3"):
                rows.append({"sample_accession": sample, "run_accession": "SRR1", "assembly_accession": "ERZ1",
                             "dataset": "2025-05", "asm_pipe_filter": "PASS", "asm_fasta_on_osf": "1",
                             "sylph_species": "Escherichia coli", "sylph_filter": "PASS", "hq_filter": "PASS",
                             "aws_url": f"https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/{sample}.fa.gz"})
            rows[1]["sylph_species"] = "Klebsiella pneumoniae"
            rows[2]["hq_filter"] = "CONTAMINATION"
            write_tsv(path, rows, list(rows[0]))
            samples = [{"biosample": x, "action": "download_atb"} for x in ("SAMN1", "SAMN2", "SAMN3", "SAMN4")]
            apply_atb_metadata(samples, path)
            self.assertEqual([r["action"] for r in samples],
                             ["download_atb", "atb_metadata_excluded", "atb_metadata_excluded", "atb_metadata_missing"])

    def test_ena_preserves_multiple_studies_and_quarantines_taxonomy_conflict(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "runs.tsv"
            records = [{"sample_accession": "SAMN1", "run_accession": "SRR1", "study_accession": "PRJ1",
                        "scientific_name": "Escherichia coli"},
                       {"sample_accession": "SAMN1", "run_accession": "SRR2", "study_accession": "PRJ2",
                        "scientific_name": "Escherichia coli"},
                       {"sample_accession": "SAMN2", "run_accession": "SRR3", "study_accession": "PRJ3",
                        "scientific_name": "Klebsiella pneumoniae"}]
            write_tsv(path, records, list(records[0]))
            (root / "manifest.json").write_text(json.dumps([{
                "path": str(path), "sha256": sha256(path), "returned_rows": 3,
                "requested_samples": ["SAMN1", "SAMN2", "SAMN3"]}]))
            samples = [{"biosample": s, "action": "download_atb"} for s in ("SAMN1", "SAMN2", "SAMN3", "SAMN4")]
            apply_ena_runs(samples, root, {})
            self.assertEqual(samples[0]["ena_study_accessions"], "PRJ1;PRJ2")
            self.assertEqual(samples[1]["action"], "quarantine_source_taxonomy_conflict")
            self.assertEqual(samples[2]["action"], "download_atb")
            self.assertEqual(samples[2]["ena_query_status"], "no_records")
            self.assertEqual(samples[3]["ena_query_status"], "not_queried")
            path.write_text("changed")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                apply_ena_runs(samples, root, {})


if __name__ == "__main__":
    unittest.main()
