import json
import tempfile
import unittest
from pathlib import Path

from strainfm import amr_task_split
from strainfm.amr_task_split import assign_components, run_amr_task_split

DRUGS = ("ciprofloxacin", "gentamicin")


def _component(i):
    stats = {"units": 2, "groups": 2, "S": 1 if i % 2 else 0, "R": 0 if i % 2 else 1}
    return {"key": f"c{i}", "units": 4, "drugs": {drug: dict(stats) for drug in DRUGS}}


class AssignComponentsTest(unittest.TestCase):
    def test_deterministic_and_constrained(self):
        comps = [_component(i) for i in range(120)]
        first = assign_components(comps, 20261110, DRUGS, min_groups=20)
        second = assign_components(comps, 20261110, DRUGS, min_groups=20)
        self.assertEqual(first, second)
        self.assertEqual(set(first.values()), {"train", "validation", "test"})

    def test_fails_when_support_insufficient(self):
        comps = [_component(i) for i in range(10)]
        with self.assertRaises(ValueError):
            assign_components(comps, 1, DRUGS, min_groups=20)


class RunTaskSplitTest(unittest.TestCase):
    def test_joint_isolation_and_shared_partition(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reps, groups = [], []
            for i in range(180):
                clone = f"g{i // 3}"
                study = f"PRJ{i // 9}"
                reps.append({"genome_id": f"u{i}", "biosample": f"u{i}",
                             "ciprofloxacin": "S" if i % 2 else "R",
                             "gentamicin": "R" if i % 3 else "S",
                             "ena_study_accessions": study,
                             "source_projects": json.dumps([study]),
                             "bioproject": ""})
                groups.append({"genome_id": f"u{i}", "near_clone_group": clone, "sequence_neighborhood": clone})
            reps_path, groups_path = root / "reps.tsv", root / "groups.tsv"
            amr_task_split.write_tsv(reps_path, reps, list(reps[0]))
            amr_task_split.write_tsv(groups_path, groups, list(groups[0]))
            summary = run_amr_task_split(reps_path, groups_path, root / "out", drugs=DRUGS,
                                         seed=7, min_groups=4)
            units = amr_task_split.read_table(root / "out" / "amr_task_units.tsv")
            by_clone, by_genome = {}, {}
            for row in units:
                by_clone.setdefault(row["near_clone_group"], set()).add(row["split"])
                by_genome.setdefault(row["genome_id"], set()).add(row["split"])
            self.assertTrue(all(len(v) == 1 for v in by_clone.values()))
            self.assertTrue(all(len(v) == 1 for v in by_genome.values()))
            self.assertEqual(summary["drugs"]["ciprofloxacin"]["n_units"], 180)
            self.assertEqual(summary["drugs"]["gentamicin"]["n_units"], 180)


if __name__ == "__main__":
    unittest.main()
