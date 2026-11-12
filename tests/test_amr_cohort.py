import gzip
import json
import tempfile
import unittest
from pathlib import Path

from strainfm.amr_cohort import run_amr_cohort
from strainfm.audit import sha256, write_tsv
from strainfm.cohort import read_table


class AmrCohortTest(unittest.TestCase):
    def test_actual_quality_exact_groups_and_all_sources_without_label_transfer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genomes = root / "genomes"
            genomes.mkdir()
            (genomes / "local.gz").write_bytes(gzip.compress(b">local\nACGTACGT\n"))
            (genomes / "SAMN2.fa.gz").write_bytes(gzip.compress(b">renamed\nACGTACGT\n"))
            (genomes / "SAMN3.fa.gz").write_bytes(gzip.compress(b">short\nAC\n"))
            digest = sha256(genomes / "local.gz")
            local = root / "local.tsv"
            write_tsv(local, [{"genome_id": "GCA_1.1", "biosample": "SAMN1", "eligible": "True",
                               "file_sha256": digest, "bioproject": "PRJNA1", "assembly_level": "Complete Genome"}],
                      ["genome_id", "biosample", "eligible", "file_sha256", "bioproject", "assembly_level"])
            samples = root / "samples.tsv"
            records = [{"genome_id": s, "biosample": s, "action": "download_atb",
                        "atb_assembly_accession": "ERZ" + str(i), "ciprofloxacin": label,
                        "atb_species": "Escherichia coli", "atb_exclusion_reason": "",
                        "ena_study_accessions": "PRJEB2;PRJEB3"}
                       for i, (s, label) in enumerate((("SAMN1", "R"), ("SAMN2", "S"), ("SAMN3", "unknown")))]
            records[0].update(action="reuse_local", local_genome_id="GCA_1.1", file_accession="GCA_1.1",
                              expected_local_sha256=digest, local_fasta_path=str(genomes / "local.gz"))
            cols = list(dict.fromkeys(k for r in records for k in r))
            write_tsv(samples, records, cols)
            config = root / "config.json"
            config.write_text(json.dumps({"genome_root": str(genomes), "atb_genome_root": str(genomes),
                                          "min_bases": 4, "max_bases": 20, "max_contigs": 2,
                                          "min_contig_n50": 4, "max_ambiguous_fraction": 0.01}))
            out = root / "cohort"
            result = run_amr_cohort(samples, local, config, out, workers=1)
            self.assertEqual((result["eligible_files"], result["identity_groups"]), (2, 1))
            rows = read_table(out / "cohort_manifest.tsv")
            self.assertEqual([r["ciprofloxacin"] for r in rows], ["R", "S", "unknown"])
            self.assertEqual(json.loads(rows[0]["source_projects"]), ["PRJEB2", "PRJEB3", "PRJNA1"])
            self.assertEqual(rows[2]["eligible"], "False")
            (genomes / "SAMN2.fa.gz").unlink()
            with self.assertRaises(FileNotFoundError):
                run_amr_cohort(samples, local, config, root / "missing", workers=1)
            self.assertFalse((root / "missing").exists())
            (genomes / "SAMN2.fa.gz").write_bytes(gzip.compress(b">x\nACGTACGT\n"))
            (genomes / "local.gz").write_bytes(gzip.compress(b">changed\nTTTTTTTT\n"))
            with self.assertRaisesRegex(ValueError, "hash changed"):
                run_amr_cohort(samples, local, config, root / "changed", workers=1)
            self.assertFalse((root / "changed").exists())
