import unittest

try:
    import numpy as np
    import sklearn  # noqa: F401
    HAVE = True
except ImportError:
    HAVE = False

from strainfm.amr_eval import _analysis_sets, _group_bootstrap, _metrics, _paired_group_bootstrap_multi


class AnalysisSetSelectionTest(unittest.TestCase):
    def test_direct_refit_requires_label_source(self):
        amr02 = [{"label_conflict": "False", "label_source": "file_accession"}]
        joint = [{"label_conflict": "False"}]
        self.assertEqual([name for name, _ in _analysis_sets(amr02)],
                         ["primary", "with_conflicts", "direct_refit"])
        self.assertEqual([name for name, _ in _analysis_sets(joint)],
                         ["primary", "with_conflicts"])


@unittest.skipUnless(HAVE, "requires numpy/scikit-learn")
class AmrEvalHelpersTest(unittest.TestCase):
    def test_metrics_single_class_returns_no_auroc(self):
        auroc, ap = _metrics([0, 0, 0], [0.1, 0.2, 0.3])
        self.assertIsNone(auroc)
        self.assertIsNotNone(ap)

    def test_bootstrap_is_deterministic(self):
        groups = ["a", "a", "b", "c", "c", "c"]
        labels = [0, 1, 0, 1, 1, 0]
        scores = [0.1, 0.9, 0.2, 0.8, 0.7, 0.3]
        first = _group_bootstrap(groups, labels, scores, 50, 20261108)
        second = _group_bootstrap(groups, labels, scores, 50, 20261108)
        self.assertEqual(first, second)
        self.assertEqual(first["n_groups"], 3)

    def test_paired_bootstrap_multi_seed(self):
        groups = ["a", "a", "b", "c", "c", "c", "d", "e"]
        labels = [0, 1, 0, 1, 1, 0, 1, 0]
        seed_scores = [
            {"model_a": [0.9, 0.8, 0.2, 0.7, 0.6, 0.3, 0.85, 0.4],
             "model_b": [0.1, 0.2, 0.8, 0.3, 0.4, 0.7, 0.15, 0.6]},
            {"model_a": [0.8, 0.7, 0.3, 0.6, 0.7, 0.2, 0.9, 0.5],
             "model_b": [0.2, 0.3, 0.7, 0.4, 0.3, 0.8, 0.1, 0.5]},
        ]
        first = _paired_group_bootstrap_multi(groups, labels, seed_scores, 100, 20261108,
                                              names=("model_a", "model_b"))
        second = _paired_group_bootstrap_multi(groups, labels, seed_scores, 100, 20261108,
                                               names=("model_a", "model_b"))
        self.assertEqual(first, second)
        self.assertIsNotNone(first["delta_AUROC_ci95"])
        with self.assertRaises(ValueError):
            _paired_group_bootstrap_multi(groups, labels, seed_scores, 10, 1, names=("model_a",))


if __name__ == "__main__":
    unittest.main()
