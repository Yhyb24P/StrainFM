import unittest

from strainfm.bacformer_embed import _flatten_shared, _nested_by_contig


class BacformerAdapterHelpersTest(unittest.TestCase):
    def test_nested_by_contig_groups_contiguous_rows(self):
        items = [(0, "c1", "AAA"), (0, "c1", "BBB"), (1, "c2", "CCC")]
        contig_ids, nested = _nested_by_contig(items)
        self.assertEqual(contig_ids, ["c1", "c2"])
        self.assertEqual(nested, [["AAA", "BBB"], ["CCC"]])

    def test_nested_by_contig_single_contig_single_protein(self):
        contig_ids, nested = _nested_by_contig([(0, "c1", "AAA")])
        self.assertEqual(contig_ids, ["c1"])
        self.assertEqual(nested, [["AAA"]])

    def test_flatten_shared_truncates_in_contig_order(self):
        flat, truncated = _flatten_shared([[1, 2], [3, 4, 5]], 3)
        self.assertEqual(flat, [1, 2, 3])
        self.assertEqual(truncated, 2)

    def test_flatten_shared_no_truncation(self):
        flat, truncated = _flatten_shared([[1], [2]], 10)
        self.assertEqual(flat, [1, 2])
        self.assertEqual(truncated, 0)


if __name__ == "__main__":
    unittest.main()
