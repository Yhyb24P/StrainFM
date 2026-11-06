import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from strainfm.snapshot import run_snapshot


class SnapshotTest(unittest.TestCase):
    def test_snapshot_exact_duplicates_missing_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genomes = root / "genomes"
            genomes.mkdir()
            for name, sequence in [("a", ">a\nACGTNN\n>b\nGC\n"), ("b", ">c\ngc\n>d\nacgtnn\n")]:
                (genomes / f"{name}.fna.gz").write_bytes(gzip.compress(sequence.encode()))
            manifest = root / "manifest.tsv"
            manifest.write_text("genome_id\tfile_accession\tspecies_id\tkingdom\ttarget_relpath\n" + "".join(f"{x}\t{x}\tSpecies_a\tbacteria\t{x}.fna.gz\n" for x in "abc"))
            config = root / "config.json"
            config.write_text(json.dumps({"species": "Species_a", "download_manifest": str(manifest), "genome_root": str(genomes)}))
            output = root / "result"
            first = run_snapshot(config, output, workers=1)
            self.assertEqual(first["file_status"], {"valid": 2, "missing": 1})
            self.assertEqual(first["exact_record_multiset_duplicate_excess"], 1)
            second = run_snapshot(config, output, workers=1)
            self.assertEqual(first["assemblies_sha256"], second["assemblies_sha256"])
            self.assertEqual(len((output / "file_checks.jsonl").read_text().splitlines()), 3)
            journal = output / "file_checks.jsonl"
            journal.write_text(journal.read_text().splitlines()[0])
            recovered = run_snapshot(config, output, workers=1)
            again = run_snapshot(config, output, workers=1)
            self.assertEqual(recovered["assemblies_sha256"], again["assemblies_sha256"])
            self.assertEqual(len(journal.read_text().splitlines()), 3)
            for line in journal.read_text().splitlines():
                json.loads(line)
            config.write_text(config.read_text() + "\n")
            with self.assertRaises(ValueError):
                run_snapshot(config, output, workers=1)

    def test_corrupt_deflate_is_a_record_not_a_job_failure(self):
        from strainfm.snapshot import inspect_file
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.fna.gz"
            content = bytearray(gzip.compress(b">a\nACGT\n"))
            content[10] = 7
            path.write_bytes(content)
            stat = path.stat()
            result = inspect_file(({"initial_size": stat.st_size, "initial_mtime_ns": stat.st_mtime_ns}, path))
            self.assertEqual(result["status"], "invalid")
