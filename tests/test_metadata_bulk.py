import json
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

from strainfm.metadata import run_bulk_metadata


class BulkMetadataTest(unittest.TestCase):
    def test_invalid_inventory_fails_before_network_or_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inventory = root / "inventory.json"
            for payload in ([], {}, [None], [{}], [{"file_accession": 1}], [{"file_accession": "GCA_000000001"}]):
                with self.subTest(payload=payload):
                    inventory.write_text(json.dumps(payload))
                    with patch("strainfm.metadata.urllib.request.urlopen") as fetch:
                        with self.assertRaises(ValueError):
                            run_bulk_metadata(inventory, root / "output")
                        fetch.assert_not_called()
                    self.assertFalse((root / "output").exists())

    def test_changed_cache_is_rejected_without_fetching_replacement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inventory = root / "inventory.json"
            inventory.write_text(json.dumps([{"file_accession": "GCA_000000001.1"}]))
            payload = {"reports": [{"accession": "GCA_000000001.1"}]}
            with patch("strainfm.metadata.urllib.request.urlopen", return_value=BytesIO(json.dumps(payload).encode())):
                result = run_bulk_metadata(inventory, root / "output")
            self.assertEqual(result["found_exact_version"], 1)
            (root / "output/batch_0000.json.gz").write_bytes(b"changed")
            with patch("strainfm.metadata.urllib.request.urlopen") as fetch:
                with self.assertRaisesRegex(ValueError, "Cached metadata batch changed"):
                    run_bulk_metadata(inventory, root / "output")
                fetch.assert_not_called()

    def test_current_version_cannot_replace_requested_old_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inventory = root / "inventory.json"
            inventory.write_text(json.dumps([{"file_accession": "GCA_000000001.1"}]))
            response = {"reports": [{"accession": "GCA_000000001.2", "assembly_info": {}}]}
            with patch("strainfm.metadata.urllib.request.urlopen", return_value=BytesIO(json.dumps(response).encode())):
                result = run_bulk_metadata(inventory, root / "output")
            self.assertEqual(result["found_exact_version"], 0)
            self.assertEqual(result["unknown"], 1)
            with patch("strainfm.metadata.urllib.request.urlopen", side_effect=AssertionError("cache should avoid another request")):
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
            with patch("strainfm.metadata.urllib.request.urlopen", side_effect=AssertionError("No requests allowed")):
                with self.assertRaisesRegex(ValueError, "unrelated nonempty"):
                    run_bulk_metadata(inventory, output)
            self.assertEqual(summary.read_text(), "keep")
