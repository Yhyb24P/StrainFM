import tempfile
import unittest
from pathlib import Path

try:
    import numpy  # noqa: F401
    import sklearn  # noqa: F401
    HAVE = True
except ImportError:
    HAVE = False

from strainfm.amr_external_eval import run_amr_external_eval
from strainfm.cohort import read_table


def _write_tsv(path, rows, columns):
    path.write_text("\t".join(columns) + "\n" + "\n".join(
        "\t".join(str(row[c]) for c in columns) for row in rows) + "\n")


@unittest.skipUnless(HAVE, "requires numpy/scikit-learn")
class ExternalEvalTest(unittest.TestCase):
    def _fixture(self, root, drop_external_for=None):
        units, predictions = [], []
        labels = [0, 1, 0, 1, 0, 1, 0, 1]
        for i, label in enumerate(labels):
            genome = f"g{i}"
            group = f"c{i // 2}"
            units.append({"drug": "ciprofloxacin", "genome_id": genome, "label": "R" if label else "S",
                          "label_conflict": "False", "near_clone_group": group, "split": "test",
                          "biosample": genome, "source_component": "s"})
            for model, base in (("esm_mean", 0.1), ("bacformer", 0.2)):
                for seed in (1, 2):
                    predictions.append({"drug": "ciprofloxacin", "model": model, "analysis_set": "primary",
                                        "seed": seed, "genome_id": genome, "near_clone_group": group,
                                        "label": label, "score": base + 0.7 * label + 0.01 * seed,
                                        "split": "test"})
        _write_tsv(root / "amr_predictions.tsv", predictions,
                   ["drug", "model", "analysis_set", "seed", "genome_id", "near_clone_group", "label", "score", "split"])
        units_path = root / "units.tsv"
        _write_tsv(units_path, units, list(units[0]))
        external = [{"genome_id": u["genome_id"], "score": 0.15 + 0.6 * (u["label"] == "R")}
                    for u in units if u["genome_id"] != drop_external_for]
        external_path = root / "external.tsv"
        _write_tsv(external_path, external, ["genome_id", "score"])
        overlap_path = root / "overlap.txt"
        overlap_path.write_text("g0\ng1\n")
        return units_path, external_path, overlap_path

    def test_sets_coverage_and_determinism(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            units_path, external_path, overlap_path = self._fixture(root)
            first = run_amr_external_eval(root, units_path, external_path, root / "out1", "strainamr",
                                          overlap_ids_path=overlap_path, bootstrap_repeats=30)
            second = run_amr_external_eval(root, units_path, external_path, root / "out2", "strainamr",
                                           overlap_ids_path=overlap_path, bootstrap_repeats=30)
            self.assertEqual(first["drugs"]["ciprofloxacin"]["primary_all"]["n_test"], 8)
            self.assertEqual(first["drugs"]["ciprofloxacin"]["primary_disjoint"]["n_test"], 6)
            metrics = read_table(root / "out1" / "external_metrics.tsv")
            self.assertEqual({(r["analysis_set"], r["model"]) for r in metrics},
                             {("primary_all", "esm_mean"), ("primary_all", "bacformer"),
                              ("primary_all", "strainamr"), ("primary_disjoint", "esm_mean"),
                              ("primary_disjoint", "bacformer"), ("primary_disjoint", "strainamr")})
            pairs = read_table(root / "out1" / "external_paired.tsv")
            self.assertEqual(len(pairs), 6)
            self.assertEqual(read_table(root / "out1" / "external_metrics.tsv"),
                             read_table(root / "out2" / "external_metrics.tsv"))

    def test_missing_external_score_is_fatal(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            units_path, external_path, overlap_path = self._fixture(root, drop_external_for="g3")
            with self.assertRaises(ValueError):
                run_amr_external_eval(root, units_path, external_path, root / "out", "strainamr",
                                      overlap_ids_path=overlap_path, bootstrap_repeats=10)

    def test_native_strainamr_format_and_duplicates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            units_path, external_path, overlap_path = self._fixture(root)
            native = root / "native.tsv"
            native.write_text("Sample_ID\tLabel\tPred\tProb\n" + "".join(
                f"g{i}\t{i % 2}\t{i % 2}\t{0.2 + 0.6 * (i % 2)}\n" for i in range(8)))
            summary = run_amr_external_eval(root, units_path, native, root / "out", "strainamr",
                                            overlap_ids_path=overlap_path, bootstrap_repeats=10)
            self.assertEqual(summary["status"], "amr_external_eval_complete")
            duplicate = root / "dup.tsv"
            duplicate.write_text("genome_id\tscore\ng0\t0.1\ng0\t0.2\n")
            with self.assertRaises(ValueError):
                run_amr_external_eval(root, units_path, duplicate, root / "out2", "strainamr",
                                      bootstrap_repeats=10)


if __name__ == "__main__":
    unittest.main()
