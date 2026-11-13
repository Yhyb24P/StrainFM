import json
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

from strainfm.metadata import run_metadata_audit
from strainfm.cohort import read_table


class MetadataAuditTest(unittest.TestCase):
    def test_required_input_fields_fail_before_network_or_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample = root / "sample.tsv"
            for content in ("file_accession\nGCA_000000001.1\n",
                            "genome_id\tfile_accession\tspecies_id\nid\tGCA_000000001.1\t\n"):
                sample.write_text(content)
                with patch("strainfm.metadata.urllib.request.urlopen") as fetch:
                    with self.assertRaises(ValueError):
                        run_metadata_audit(sample, root / "output")
                    fetch.assert_not_called()
                self.assertFalse((root / "output").exists())

    def test_shared_fields_preserve_audit_aliases_and_missing_values(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample = root / "sample.tsv"
            sample.write_text("genome_id\tfile_accession\tspecies_id\ttotal_bases\nGCF_000000001.1\tGCA_000000001.1\tExample_species\t4\n")
            payload = {"reports": [{"accession": "GCA_000000001.1",
                                   "organism": {"organism_name": "Example species", "tax_id": 123},
                                   "assembly_info": {"paired_assembly": {"accession": "GCF_000000001.1"},
                                                     "biosample": {"accession": "SAMN1", "attributes": [
                                                         {"name": "geographic location", "value": "France"},
                                                         {"name": "collection_date", "value": "not collected"}]}},
                                   "assembly_stats": {"total_sequence_length": 4, "contig_n50": 4}}]}
            with patch("strainfm.metadata.urllib.request.urlopen", return_value=BytesIO(json.dumps(payload).encode())):
                result = run_metadata_audit(sample, root / "output")
            row = read_table(root / "output/sample_metadata.tsv")[0]
            self.assertEqual(result["status"], "completed")
            self.assertEqual(row["reported_organism"], "Example species")
            self.assertEqual(row["biosample"], "SAMN1")
            self.assertEqual(row["ncbi_contig_n50"], "4")
            for field in ("location_known", "original_id_is_returned_or_paired", "fasta_total_length_matches"):
                self.assertEqual(row[field], "True")
            self.assertEqual(row["collection_date_known"], "False")
            self.assertEqual(row["annotated_protein_coding_genes"], "")

    def test_empty_response_is_partial_not_complete(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample = root / "sample.tsv"
            sample.write_text("genome_id\tfile_accession\tspecies_id\ttotal_bases\nGCF_000000001.1\tGCA_000000001.1\tExample_species\t4\n")
            with patch("strainfm.metadata.urllib.request.urlopen", return_value=BytesIO(json.dumps({"reports": []}).encode())):
                result = run_metadata_audit(sample, root / "output")
            self.assertEqual(result["status"], "partial")
            self.assertEqual(result["species"][0]["requested"], 1)
            self.assertEqual(result["species"][0]["found"], 0)

    def test_wrong_returned_version_remains_partial(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample = root / "sample.tsv"
            sample.write_text("genome_id\tfile_accession\tspecies_id\ttotal_bases\nGCF_000000001.1\tGCA_000000001.1\tExample_species\t4\n")
            payload = {"reports": [{"accession": "GCA_000000001.2"}]}
            with patch("strainfm.metadata.urllib.request.urlopen", return_value=BytesIO(json.dumps(payload).encode())):
                result = run_metadata_audit(sample, root / "output")
            self.assertEqual(result["status"], "partial")
            self.assertEqual(result["species"][0]["found"], 0)
