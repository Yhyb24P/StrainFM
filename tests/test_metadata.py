import json
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

from strainfm.metadata import run_metadata_audit


class MetadataAuditTest(unittest.TestCase):
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
