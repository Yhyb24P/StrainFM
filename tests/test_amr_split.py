import unittest

from strainfm.amr_split import assign_groups


class AmrSplitTest(unittest.TestCase):
    def test_deterministic_and_group_isolated(self):
        groups = {f"g{i}": ["S" if i % 2 else "R"] for i in range(300)}
        first = assign_groups(groups, 20261108, "gentamicin", min_groups=20)
        second = assign_groups(groups, 20261108, "gentamicin", min_groups=20)
        self.assertEqual(first, second)
        self.assertEqual(set(first), set(groups))
        counts = {split: sum(1 for name in first.values() if name == split)
                  for split in ("train", "validation", "test")}
        for split, count in counts.items():
            self.assertGreaterEqual(count, 20, split)

    def test_insufficient_support_raises(self):
        groups = {f"g{i}": ["S"] for i in range(10)}
        with self.assertRaises(ValueError):
            assign_groups(groups, 1, "gentamicin", min_groups=20)

    def test_seed_changes_assignment(self):
        groups = {f"g{i}": ["S"] for i in range(300)}
        a = assign_groups(groups, 1, "drug", min_groups=20)
        b = assign_groups(groups, 2, "drug", min_groups=20)
        self.assertNotEqual(a, b)


if __name__ == "__main__":
    unittest.main()
