import json
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

from strainfm.metadata_bulk import run_bulk_metadata


class BulkMetadataTest(unittest.TestCase):
    def test_current_version_cannot_replace_requested_old_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inventory = root / "inventory.json"
            inventory.write_text(json.dumps([{"file_accession": "GCA_000000001.1"}]))
            response = {"reports": [{"accession": "GCA_000000001.2", "assembly_info": {}}]}
            with patch("strainfm.metadata_bulk.urllib.request.urlopen", return_value=BytesIO(json.dumps(response).encode())):
                result = run_bulk_metadata(inventory, root / "output")
            self.assertEqual(result["found_exact_version"], 0)
            self.assertEqual(result["unknown"], 1)
            with patch("strainfm.metadata_bulk.urllib.request.urlopen", side_effect=AssertionError("cache should avoid another request")):
                cached = run_bulk_metadata(inventory, root / "output")
            self.assertEqual(cached["unknown"], 1)

    def test_unowned_nonempty_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inventory = root / "inventory.json"
            inventory.write_text(json.dumps([{"file_accession": "GCA_000000001.1"}]))
            output = root / "output"
            output.mkdir()
            summary = output / "summary.json"
            summary.write_text("keep")
            with patch("strainfm.metadata_bulk.urllib.request.urlopen", side_effect=AssertionError("No requests allowed")):
                with self.assertRaisesRegex(ValueError, "unrelated nonempty"):
                    run_bulk_metadata(inventory, output)
            self.assertEqual(summary.read_text(), "keep")
