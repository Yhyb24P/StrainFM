import csv
import tempfile
import unittest
from pathlib import Path

from strainfm.labels import run_label_audit


class LabelAuditTest(unittest.TestCase):
    def test_exact_versions_and_missing_labels(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, labels = root / "manifest.tsv", root / "labels.csv"
            manifest.write_text("genome_id\tfile_accession\tspecies_id\tkingdom\nGCF_1.1\tGCA_1.1\tSpecies_a\tbacteria\n")
            labels.write_text("genome_name,drug_a,drug_b\nGCA_1.1,1.0,\nGCA_1.2,0.0,1.0\n")
            result = run_label_audit(labels, manifest, root / "output")
            self.assertEqual(result["matching"], {"file_accession_only": 1, "unmatched": 1})
            self.assertEqual(result["global"]["missing_cells"], 1)
            with (root / "output/species_drug_counts.tsv").open() as handle:
                records = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0]["n_label_0"], "0")
            with self.assertRaises(FileExistsError):
                run_label_audit(labels, manifest, root / "output")
