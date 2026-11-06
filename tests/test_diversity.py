import gzip
import json
import random
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from strainfm.audit import sha256, write_tsv
from strainfm.diversity import run_diversity


MASH = Path.home() / ".local/share/strainfm/tools/mash-2.3/mash-Linux64-v2.3/mash"


@unittest.skipUnless(MASH.is_file(), "optional Mash 2.3 binary not installed")
class DiversityTest(unittest.TestCase):
    def test_real_mash_exact_duplicates_and_unrelated_genome(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = random.Random(42)
            first = "".join(rng.choices("ACGT", k=20000))
            unrelated = "".join(rng.choices("ACGT", k=20000))
            samples = []
            for name, sequence in [("a", first), ("b", first), ("c", unrelated)]:
                path = root / f"{name}.fna.gz"
                path.write_bytes(gzip.compress(f">{name}\n{sequence}\n".encode()))
                stat = path.stat()
                samples.append({"genome_id": name, "target_relpath": path.name, "initial_size": stat.st_size, "initial_mtime_ns": stat.st_mtime_ns, "file_sha256": sha256(path)})
            write_tsv(root / "representatives.tsv", samples, list(samples[0]))
            config = {"mash_binary": str(MASH), "genome_root": str(root), "sketch_batch_size": 2, "k": 21, "sketch_size": 1000, "seed": 42, "near_clone_distance": 0.001, "neighborhood_distance": 0.005}
            (root / "config.json").write_text(json.dumps(config))
            (root / "unowned").mkdir()
            sentinel = root / "unowned/summary.json"
            sentinel.write_text("keep")
            with self.assertRaisesRegex(ValueError, "unrelated nonempty"):
                run_diversity(root / "representatives.tsv", root / "config.json", root / "unowned", threads=1)
            self.assertEqual(sentinel.read_text(), "keep")
            result = run_diversity(root / "representatives.tsv", root / "config.json", root / "output", threads=1)
            self.assertEqual(result["near_clone_groups"], 2)
            self.assertEqual(result["sequence_neighborhoods"], 2)
            self.assertEqual(result["near_clone_edges"], 1)
            self.assertEqual(result["near_clone_forest_edges"], result["records"] - result["near_clone_groups"])
            with patch("strainfm.diversity.subprocess.run", wraps=subprocess.run) as calls:
                repeated = run_diversity(root / "representatives.tsv", root / "config.json", root / "output", threads=1)
            self.assertEqual([call.args[0][1] for call in calls.call_args_list], ["--version"])
            self.assertEqual(result, repeated)
            # A crash after paste must also resume without manually deleting genomes.msh.
            (root / "output/completion.json").unlink()
            recovered = run_diversity(root / "representatives.tsv", root / "config.json", root / "output", threads=1)
            self.assertEqual(recovered["groups_sha256"], result["groups_sha256"])
            (root / "output/sequence_groups.tsv").write_text("corrupted")
            with self.assertRaisesRegex(ValueError, "Completed artifact changed"):
                run_diversity(root / "representatives.tsv", root / "config.json", root / "output", threads=1)

    def test_content_change_with_same_stat_is_rejected(self):
        import os
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "a.fna.gz"
            path.write_bytes(gzip.compress(b">a\nACGT\n", compresslevel=0, mtime=0))
            before = path.stat()
            record = {"genome_id": "a", "target_relpath": path.name, "initial_size": before.st_size, "initial_mtime_ns": before.st_mtime_ns, "file_sha256": sha256(path)}
            write_tsv(root / "representatives.tsv", [record], list(record))
            path.write_bytes(gzip.compress(b">a\nTTTT\n", compresslevel=0, mtime=0))
            os.utime(path, ns=(before.st_mtime_ns, before.st_mtime_ns))
            self.assertEqual(path.stat().st_size, before.st_size)
            config = {"mash_binary": str(MASH), "genome_root": str(root), "sketch_batch_size": 2, "near_clone_distance": .001, "neighborhood_distance": .005}
            (root / "config.json").write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, "Input content changed"):
                run_diversity(root / "representatives.tsv", root / "config.json", root / "output", threads=1)

    def test_invalid_thresholds_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = {"near_clone_distance": .02, "neighborhood_distance": .005}
            (root / "config.json").write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, "near_clone_distance"):
                run_diversity(root / "nonexistent.tsv", root / "config.json", root / "output")
