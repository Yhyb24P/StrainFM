import csv
import json
import tempfile
import unittest
from pathlib import Path

from strainfm.report import run_report


class ReportTest(unittest.TestCase):
    def test_report_preserves_estimators_and_links_every_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "ref").mkdir()
            references = [{"id": "example", "kind": "软件", "title": "A < B", "authors": "Authors",
                           "year": "2026", "url": "https://example.org", "explanation": "说明",
                           "project": "用途", "limit": "边界"}]
            (root / "ref/references.json").write_text(json.dumps(references))
            saved = root / "data/results/run"
            saved.mkdir(parents=True)
            columns = ["drug", "analysis_set", "model", "seed", "n_test", "test_R", "AUROC", "average_precision"]
            with (saved / "amr_metrics.tsv").open("w") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
                writer.writeheader()
                for seed, value in ((1, 0.5), (2, 0.9)):
                    writer.writerow(dict(zip(columns, ["drug", "primary", "bacformer", seed, 10, 4, value, value])))
            with (saved / "external_metrics.tsv").open("w") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
                writer.writeheader()
                writer.writerow(dict(zip(columns, ["drug", "primary_disjoint", "bacformer", "", 10, 4, 0.85, 0.86])))
            result = run_report(root)
            catalog = json.loads((root / "data/results/catalog.json").read_text())
            records = {record["analysis"]: record for record in catalog["records"]}
            self.assertAlmostEqual(records["primary"]["metrics"]["AUROC"], 0.7)
            self.assertEqual(records["primary_disjoint"]["metrics"]["AUROC"], 0.85)
            self.assertNotEqual(records["primary"]["estimator"], records["primary_disjoint"]["estimator"])
            self.assertEqual(result["references"], 1)
            self.assertIn('id="example"', (root / "ref/index.html").read_text())
            self.assertIn('A &lt; B', (root / "ref/index.html").read_text())
            self.assertIn('@misc{example', (root / "ref/references.bib").read_text())
            self.assertEqual(len(catalog["input_sha256"]), 2)
            original = (root / "data/results/catalog.json").read_bytes()
            run_report(root)
            self.assertEqual(original, (root / "data/results/catalog.json").read_bytes())
            source = saved / "amr_metrics.tsv"
            source.write_text(source.read_text().replace("0.9", "nan"))
            with self.assertRaisesRegex(ValueError, "Invalid AUROC"):
                run_report(root)
            self.assertEqual(original, (root / "data/results/catalog.json").read_bytes())
