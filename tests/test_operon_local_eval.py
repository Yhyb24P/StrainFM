import json
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_SCIENTIFIC = True
except ImportError:
    HAVE_SCIENTIFIC = False

from strainfm import operon_local_eval


@unittest.skipUnless(HAVE_SCIENTIFIC, "requires numpy/pandas/pyarrow")
class OperonLocalEvalTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.emb = self.root / "emb.parquet"
        self.mapping = self.root / "mapping.tsv"
        self.reference = self.root / "ref.parquet"
        import pandas as pd
        import pyarrow as pa
        import pyarrow.parquet as pq

        n = 6
        emb = np.arange(n * 3, dtype=np.float32).reshape(n, 3)
        row = {"contig_name": ["c1"], "embeddings": [emb.tolist()], "strand": [[1, 1, 1, 1, 1, 1]]}
        pq.write_table(pa.Table.from_pandas(pd.DataFrame([row])), self.emb)
        self.mapping.write_text(
            "author_index\tlocal_index\n0\t0\n1\t1\n2\t3\n3\t4\n4\t5\n5\t2\n"
        )
        pq.write_table(pa.Table.from_pandas(pd.DataFrame([{
            "protein_sequence": [["A", "B", "C", "D", "E", "F"]],
            "operon_prot_indices": [[[0, 1], [2]]],
        }])), self.reference)

        def labels_fn(operons, n_genes):
            labels = np.zeros(n_genes - 1, dtype=int)
            for operon in operons:
                for i in range(min(operon), max(operon)):
                    labels[i] = 1
            return labels

        def scores_fn(emb, strand):
            return np.linspace(0, 1, len(emb) - 1)

        operon_local_eval._official_functions = lambda: (labels_fn, scores_fn)

    def tearDown(self):
        self.tmp.cleanup()

    def test_aligned_pairs_and_exclusions(self):
        out = self.root / "out"
        summary = operon_local_eval.run_operon_local_eval(self.emb, self.mapping, self.reference, out, "fake", "rev1")
        self.assertEqual(summary["n_local_proteins"], 6)
        self.assertEqual(summary["n_aligned_pairs"], 3)
        self.assertEqual(summary["excluded_pairs"], {"other_non_consecutive": 2})
        self.assertEqual(summary["n_positive_pairs"], 1)
        self.assertTrue((out / "aligned_pairs.tsv").exists())
        self.assertEqual(json.loads((out / "metrics.json").read_text())["method"], "fake")


if __name__ == "__main__":
    unittest.main()
