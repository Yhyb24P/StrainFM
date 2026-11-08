import gzip
import json
import os
import shutil
import struct
import subprocess
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from strainfm.audit import write_tsv
from strainfm.cohort import read_table
from strainfm.eval_distance import load_forest_distances, mash_cli_batch, run_eval_distance, sample_distance_pairs, spearman_rho
from strainfm.eval_retrieval import cosine_distance, mean_average_precision, run_eval_retrieval

MASH = Path(os.environ.get("STRAINFM_MASH", str(Path.home() / ".local/share/strainfm/tools/mash-2.3/mash-Linux64-v2.3/mash")))
HAVE_MASH = MASH.exists()

def _pseudo_dna(seed, length=20000):
    # random.Random is deterministic across runs/platforms; long enough to hold
    # thousands of distinct 21-mers so k=21 sketches (s=50) do not degenerate
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


FASTA_ALPHABET = _pseudo_dna(7)
TEST_SKETCH_SIZE = 50


def _write_fasta(path, name, repeat, count=3):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for i in range(count):
            handle.write(f">{name}_{i}\n")
            handle.write(((repeat * 2)[:10000]) + "\n")


def _setup(tmp, model="bacformer"):
    root = Path(tmp)
    manifest = [
        {"genome_id": "t1", "split": "train", "assembly_level": "Complete Genome"},
        {"genome_id": "t2", "split": "train", "assembly_level": "Chromosome"},
        {"genome_id": "t3", "split": "train", "assembly_level": "Contig"},
        {"genome_id": "t4", "split": "train", "assembly_level": "Scaffold"},
        {"genome_id": "t5", "split": "train", "assembly_level": "Contig"},
        {"genome_id": "t6", "split": "train", "assembly_level": "Contig"},
        {"genome_id": "v1", "split": "validation", "assembly_level": "Complete Genome"},
        {"genome_id": "v2", "split": "validation", "assembly_level": "Chromosome"},
        {"genome_id": "v3", "split": "validation", "assembly_level": "Contig"},
        {"genome_id": "v4", "split": "validation", "assembly_level": "Chromosome"},
        {"genome_id": "v5", "split": "validation", "assembly_level": "Contig"},
        {"genome_id": "v6", "split": "validation", "assembly_level": "Chromosome"},
        {"genome_id": "v7", "split": "validation", "assembly_level": "Chromosome"},
        {"genome_id": "x1", "split": "test", "assembly_level": "Complete Genome"},
    ]
    groups = [
        {"genome_id": "t1", "near_clone_group": "t1", "sequence_neighborhood": "N0"},
        {"genome_id": "t2", "near_clone_group": "t2", "sequence_neighborhood": "N0"},
        {"genome_id": "t3", "near_clone_group": "t3", "sequence_neighborhood": "N1"},
        {"genome_id": "t4", "near_clone_group": "t4", "sequence_neighborhood": "N2"},
        {"genome_id": "t5", "near_clone_group": "t5", "sequence_neighborhood": "N3"},
        {"genome_id": "t6", "near_clone_group": "t6", "sequence_neighborhood": "N4"},
        {"genome_id": "v1", "near_clone_group": "v1", "sequence_neighborhood": "N0"},
        {"genome_id": "v2", "near_clone_group": "v2", "sequence_neighborhood": "N1"},
        {"genome_id": "v3", "near_clone_group": "v3", "sequence_neighborhood": "N3"},
        {"genome_id": "v4", "near_clone_group": "v4", "sequence_neighborhood": "N5"},
        {"genome_id": "v5", "near_clone_group": "v5", "sequence_neighborhood": "N5"},
        {"genome_id": "v6", "near_clone_group": "v6", "sequence_neighborhood": "N5"},
        {"genome_id": "v7", "near_clone_group": "v7", "sequence_neighborhood": "N5"},
        {"genome_id": "x1", "near_clone_group": "x1", "sequence_neighborhood": "N0"},
    ]
    support = [
        {"genome_id": "v1", "split": "validation", "positive_definition": "same_Mash_sequence_neighborhood", "gallery_split": "train", "gallery_positive_count": "2", "eligible_for_positive_supported_retrieval": "True"},
        {"genome_id": "v2", "split": "validation", "positive_definition": "same_Mash_sequence_neighborhood", "gallery_split": "train", "gallery_positive_count": "1", "eligible_for_positive_supported_retrieval": "True"},
        {"genome_id": "v3", "split": "validation", "positive_definition": "same_Mash_sequence_neighborhood", "gallery_split": "train", "gallery_positive_count": "0", "eligible_for_positive_supported_retrieval": "False"},
        {"genome_id": "v4", "split": "validation", "positive_definition": "same_Mash_sequence_neighborhood", "gallery_split": "train", "gallery_positive_count": "0", "eligible_for_positive_supported_retrieval": "False"},
        {"genome_id": "v6", "split": "validation", "positive_definition": "same_Mash_sequence_neighborhood", "gallery_split": "train", "gallery_positive_count": "0", "eligible_for_positive_supported_retrieval": "False"},
        {"genome_id": "v7", "split": "validation", "positive_definition": "same_Mash_sequence_neighborhood", "gallery_split": "train", "gallery_positive_count": "0", "eligible_for_positive_supported_retrieval": "False"},
        {"genome_id": "x1", "split": "test", "positive_definition": "same_Mash_sequence_neighborhood", "gallery_split": "train", "gallery_positive_count": "2", "eligible_for_positive_supported_retrieval": "True"},
    ]
    write_tsv(root / "split_manifest.tsv", manifest, ["genome_id", "split", "assembly_level"])
    write_tsv(root / "sequence_groups.tsv", groups, ["genome_id", "near_clone_group", "sequence_neighborhood"])
    write_tsv(root / "support.tsv", support, list(support[0]))
    bases = {
        "N0": [1.0, 0.0], "N1": [0.0, 1.0], "N2": [1.0, 1.0], "N3": [1.0, -1.0],
        "N4": [-1.0, 1.0], "N5": [-1.0, -1.0],
    }
    neighborhoods = {r["genome_id"]: r["sequence_neighborhood"] for r in groups}
    emb_dir = root / "embeds"
    emb_dir.mkdir()
    rows = []
    for gid in ("t1", "t2", "t3", "t4", "t5", "t6", "v1", "v2", "v3", "v4", "v6", "v7", "x1"):
        base = bases[neighborhoods[gid]]
        vec = [base[0] + 0.01, base[1] + (0.01 if gid != "v4" else -0.01)]
        rows.append({"genome_id": gid, "embedding": json.dumps(vec), "truncated_protein_count": "0"})
    if model == "bacformer":
        write_tsv(emb_dir / "genome_embeddings.tsv", rows, ["genome_id", "embedding", "truncated_protein_count"])
        write_tsv(emb_dir / "embedding_failures.tsv", [{"genome_id": "v5", "split": "validation", "error": "test failure"}], ["genome_id", "split", "error"])
    else:
        write_tsv(emb_dir / "esm_mean_embeddings.tsv", rows, ["genome_id", "embedding"])

    # Real sketch artifacts for the distance path: near-identical genomes for
    # forest-resolved pairs, disjoint alphabets for far pairs.
    root_dir = root / "raw" / "genomes"
    repeats = {
        "t1": FASTA_ALPHABET, "t2": FASTA_ALPHABET, "v1": FASTA_ALPHABET, "x1": FASTA_ALPHABET,
        "t3": FASTA_ALPHABET, "v2": FASTA_ALPHABET,
        "t4": _pseudo_dna(81), "t5": _pseudo_dna(82), "t6": _pseudo_dna(83),
        "v3": FASTA_ALPHABET, "v4": FASTA_ALPHABET, "v6": FASTA_ALPHABET, "v7": FASTA_ALPHABET,
    }
    assemblies = []
    for gid, rep in repeats.items():
        rel = f"{gid}.fna"
        _write_fasta(root_dir / rel, gid, rep)
        assemblies.append({"genome_id": gid, "target_relpath": rel, "status": "valid"})
    write_tsv(root / "assemblies.tsv", assemblies, ["genome_id", "target_relpath", "status"])
    msh_dir = root / "groups"
    msh_dir.mkdir()
    (msh_dir / "sequence_groups.tsv").write_bytes((root / "sequence_groups.tsv").read_bytes())
    if HAVE_MASH:
        inputs = "\n".join(str(root_dir / f"{gid}.fna") for gid in repeats) + "\n"
        (root / "inputs.txt").write_text(inputs)
        list_file = msh_dir / "inputs.txt"
        list_file.write_text(inputs)
        subprocess.run([str(MASH), "sketch", "-l", "-k", "21", "-s", str(TEST_SKETCH_SIZE), "-S", "42",
                        "-o", str(msh_dir / "genomes"), str(list_file)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert (msh_dir / "genomes.msh").exists()
    # Stored neighborhoods forest: v1-t1-t2 chain (v1 not directly adjacent to t1),
    # v2-t3 edge; far genomes isolated.
    forest_lines = "genome_a\tgenome_b\tmash_distance\tshared_hashes\n"
    forest_lines += "t1\tt2\t0.001\t990\nv1\tt2\t0.0012\t988\nv2\tt3\t0.002\t980\n"
    with gzip.open(msh_dir / "neighborhood_forest.tsv.gz", "wt") as handle:
        handle.write(forest_lines)
    return root


class EvalRetrievalTest(unittest.TestCase):
    def test_protocol_a_metrics_open_set_and_tie_break(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = _setup(tmp)
            out = root / "out_a"
            summary = run_eval_retrieval(root / "embeds", root / "split_manifest.tsv", root / "support.tsv",
                                         root / "sequence_groups.tsv", out,
                                         model="bacformer", protocol="A", split="validation")
            # protocol A now honours --split: validation only = v1, v2, v3, v4, v6, v7 (v5 failed)
            self.assertEqual(summary["n_queries_total"], 6)
            # open set: v4/v6/v7 (neighborhood N5 has no train member)
            self.assertEqual(summary["open_set_n"], 3)
            self.assertEqual(summary["n_queries_scored"], 3)
            self.assertEqual(summary["all_hit"]["hit@1"], 1.0)
            self.assertAlmostEqual(summary["all_mrr"], 1.0)
            metric_rows = list(read_table(out / "retrieval_metrics.tsv"))
            strata = {r["stratum"] for r in metric_rows}
            self.assertIn("ALL", strata)
            all_row = next(r for r in metric_rows if r["stratum"] == "ALL")
            self.assertEqual(all_row["open_set_n"], "3")
            self.assertTrue((out / "retrieval_results.tsv").exists() or (out / "retrieval_results.parquet").exists())
            contract = json.loads((out / "contract.json").read_text())
            self.assertEqual(contract["positive_definition"], "same_Mash_sequence_neighborhood")
            self.assertEqual(contract["support_eligible_no_gallery_positive"], 0)
            self.assertEqual(contract["support_positive_count_mismatch"], 0)
            # Fixture support marks v3 ineligible although N3 has train member t5
            self.assertEqual(contract["support_ineligible_with_positives"], 1)

    def test_protocol_a_split_filter_isolation(self):
        """Synthetic-data check of the split filter: validation run must exclude the
        test query x1, and a test run must contain only x1. The real test split is
        not touched (fixtures only)."""
        def result_queries(out):
            tsv = out / "retrieval_results.tsv"
            if tsv.exists():
                return {r["query_genome_id"] for r in read_table(tsv)}
            import pyarrow.parquet as pq
            return set(pq.read_table(out / "retrieval_results.parquet", columns=["query_genome_id"]).column(0).to_pylist())
        with tempfile.TemporaryDirectory() as tmp:
            root = _setup(tmp)
            out_v = root / "out_v"
            run_eval_retrieval(root / "embeds", root / "split_manifest.tsv", root / "support.tsv",
                               root / "sequence_groups.tsv", out_v,
                               model="bacformer", protocol="A", split="validation")
            self.assertNotIn("x1", result_queries(out_v))
            self.assertEqual(json.loads((out_v / "contract.json").read_text())["n_queries_missing_embeddings"], 0)
            out_t = root / "out_t"
            summary_t = run_eval_retrieval(root / "embeds", root / "split_manifest.tsv", root / "support.tsv",
                                           root / "sequence_groups.tsv", out_t,
                                           model="bacformer", protocol="A", split="test")
            self.assertEqual(summary_t["n_queries_total"], 1)
            self.assertEqual(summary_t["open_set_n"], 0)
            self.assertEqual(summary_t["n_queries_scored"], 1)
            self.assertEqual(result_queries(out_t), {"x1"})

    def test_protocol_b_neighborhood_split_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = _setup(tmp)
            out = root / "out_b"
            summary = run_eval_retrieval(root / "embeds", root / "split_manifest.tsv", root / "support.tsv",
                                         root / "sequence_groups.tsv", out,
                                         model="bacformer", protocol="B", split="validation", seed=20261106)
            # only N5 has >=2 usable validation members (v4, v6, v7); one query, 2 gallery positives
            self.assertEqual(summary["n_queries_scored"], 1)
            self.assertEqual(summary["open_set_n"], 0)
            again = run_eval_retrieval(root / "embeds", root / "split_manifest.tsv", root / "support.tsv",
                                       root / "sequence_groups.tsv", root / "out_b2",
                                       model="bacformer", protocol="B", split="validation", seed=20261106)
            self.assertEqual(summary["all_hit"], again["all_hit"])

    def test_protocol_a_refuses_train_split(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = _setup(tmp)
            with self.assertRaises(ValueError):
                run_eval_retrieval(root / "embeds", root / "split_manifest.tsv", root / "support.tsv",
                                   root / "sequence_groups.tsv", root / "out_x",
                                   model="bacformer", protocol="A", split="train")


class EvalDistanceTest(unittest.TestCase):
    def test_forest_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with gzip.open(root / "f.tsv.gz", "wt") as handle:
                handle.write("genome_a\tgenome_b\tmash_distance\tshared_hashes\na\tb\t0.01\t990\nc\ta\t0.5\t100\n")
            distances, adjacency = load_forest_distances(root / "f.tsv.gz")
            self.assertAlmostEqual(distances[("a", "b")], 0.01)
            self.assertAlmostEqual(distances[("a", "c")], 0.5)
            self.assertEqual(adjacency["a"], {"b", "c"})

    def test_spearman(self):
        self.assertAlmostEqual(spearman_rho([1, 2, 3, 4], [2, 1, 4, 3]), 0.6)
        self.assertAlmostEqual(spearman_rho([1, 2, 3], [2.0, 4.0, 6.0]), 1.0)
        self.assertIsNone(spearman_rho([1, 1, 1], [1, 2, 3]))
        self.assertIsNone(spearman_rho([1, 2], [2, 1]))

    def test_distance_bin_selection(self):
        import random
        rng = random.Random(1)
        # 0.014 vs reference 0.01: 0.004 apart, inside the +/-0.005 width
        reference, selected = sample_distance_pairs(
            "q", [0.01, 0.02], [("a", 0.015), ("b", 0.05), ("c", 0.014)], rng)
        self.assertAlmostEqual(reference, 0.01)
        self.assertEqual([p for p, _ in selected], ["a", "c"])

    @unittest.skipUnless(HAVE_MASH, "mash binary not installed")
    def test_mash_cli_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = root / "a.fna"
            b = root / "b.fna"
            _write_fasta(a, "a", FASTA_ALPHABET)
            _write_fasta(b, "b", _pseudo_dna(99))
            sketch = root / "s.msh"
            list_file = root / "inputs.txt"
            list_file.write_text(f"{a}\n{b}\n")
            subprocess.run([str(MASH), "sketch", "-l", "-k", "21", "-s", str(TEST_SKETCH_SIZE),
                            "-o", str(root / "s"), str(list_file)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            assert sketch.exists()
            resolved = mash_cli_batch(MASH, sketch, [str(a), str(b)])
            self.assertAlmostEqual(resolved[str(a)], 0.0)
            self.assertGreater(resolved[str(b)], 0.1)

    @unittest.skipUnless(HAVE_MASH, "mash binary not installed")
    def test_run_eval_distance_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = _setup(tmp)
            out = root / "out_d"
            summary = run_eval_distance(root / "embeds", root / "split_manifest.tsv", root / "groups",
                                        root / "assemblies.tsv", root / "raw" / "genomes", MASH, out,
                                        model="bacformer", split="validation", per_bin=5)
            self.assertEqual(summary["n_queries"], 6)
            self.assertGreater(summary["n_pairs"], 0)
            self.assertIn("positive_same_neighborhood", summary["bin_counts"])
            # only v1, v2 have forest-resolvable positive references; v3's
            # positive t5 has no stored edge, so its far bin is skipped
            self.assertEqual(summary["n_mash_cli_batches"], 2)
            # v1's positive t1 has no direct forest edge (v1-t2-t1 chain): dropped, not imputed
            self.assertGreaterEqual(summary["n_unresolved_pairs_dropped"], 1)
            correlation = (out / "distance_correlation.tsv").read_text()
            self.assertIn("ALL", correlation)


class MetricHelperTest(unittest.TestCase):
    def test_cosine_distance(self):
        self.assertAlmostEqual(cosine_distance([1, 0], [1, 0]), 0.0)
        self.assertAlmostEqual(cosine_distance([1, 0], [0, 1]), 1.0)
        with self.assertRaises(ValueError):
            cosine_distance([0, 0], [1, 0])

    def test_map(self):
        self.assertAlmostEqual(mean_average_precision([True, False, True], 3, 2), (1 + 2 / 3) / 2)
        self.assertEqual(mean_average_precision([False, False], 2, 0), 0.0)
        # Denominator min(k, n_positives): second positive outside top-10 halves AP
        self.assertAlmostEqual(mean_average_precision([True] + [False] * 9 + [True], 10, 2), 0.5)


if __name__ == "__main__":
    unittest.main()
