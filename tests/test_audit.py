import gzip
import tempfile
import unittest
from pathlib import Path

from strainfm.audit import fasta_stats, rows, safe_path


class FastaAuditTest(unittest.TestCase):
    def test_stats_and_order_insensitive_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory) / "a.fna.gz", Path(directory) / "b.fna.gz"
            with gzip.open(a, "wt") as handle:
                handle.write(">one\nACGTNN\n>two\nGC\n")
            with gzip.open(b, "wt") as handle:
                handle.write(">renamed\ngc\n>another\nacgt\nnn\n")
            first, second = fasta_stats(a), fasta_stats(b)
            self.assertEqual(first["n50_bases"], 6)
            self.assertEqual(first["total_bases"], 8)
            self.assertEqual(first["sequence_count"], 2)
            self.assertEqual(first["n_fraction"], 0.25)
            self.assertAlmostEqual(first["gc_fraction_acgt"], 4 / 6)
            self.assertEqual(first["contig_multiset_sha256"], second["contig_multiset_sha256"])

    def test_truncated_gzip_and_invalid_fasta_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.fna.gz"
            path.write_bytes(gzip.compress(b">x\nACGT\n")[:-5])
            with self.assertRaises((EOFError, OSError)):
                fasta_stats(path)
            path.write_bytes(gzip.compress(b"ACGT\n"))
            with self.assertRaises(ValueError):
                fasta_stats(path)

    def test_path_traversal_rejected(self):
        with self.assertRaises(ValueError):
            safe_path(Path("/tmp/genomes"), "../outside.fna.gz")

    def test_trailing_empty_column_counted_but_nonempty_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.tsv"
            path.write_text("id\tname\nx\ty\t\n")
            self.assertEqual(list(rows(path, ("id", "name")))[0]["_trailing_empty_fields"], 1)
            path.write_text("id\tname\nx\ty\tz\n")
            with self.assertRaises(ValueError):
                list(rows(path, ("id", "name")))


if __name__ == "__main__":
    unittest.main()
