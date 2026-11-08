import unittest

try:
    import numpy as np
    import sklearn  # noqa: F401
    HAVE = True
except ImportError:
    HAVE = False

from strainfm.amr_eval import _group_bootstrap, _metrics


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


if __name__ == "__main__":
    unittest.main()
