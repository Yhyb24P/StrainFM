import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from strainfm.amr_intersection import _read_labels, run_amr_intersection
from strainfm.audit import sha256, write_tsv
from strainfm.cohort import read_table


class AmrFeasibilityTest(unittest.TestCase):
    def fixture(self, root):
        genomes = root / 'genomes'
        genomes.mkdir()
        cohort, groups, split = [], [], []
        for i in range(7):
            accession = f'GCA_{i}.1'
            gid = f'GCF_{i}.1' if i == 0 else accession
            path = genomes / f'{accession}.fna.gz'
            path.write_bytes(gzip.compress(b'>a\nACGT\n'))
            rep = 'unlabeled_identity' if i == 2 else gid
            row = {'genome_id': gid, 'file_accession': accession, 'current_accession': accession, 'paired_accession': gid if i == 0 else '',
                   'biosample': f'BS{i}', 'bioproject': f'P{i}', 'file_sha256': sha256(path), 'target_relpath': path.name,
                   'eligible': str(i < 6), 'status': 'valid', 'exclusion_reason': '' if i < 6 else 'contig_n50_below_minimum',
                   'identity_group': rep, 'identity_representative': rep, 'retained_identity_representative': str(i != 2 and i < 6)}
            cohort.append(row)
            if i < 6:
                groups.append({'genome_id': rep, 'near_clone_group': f'N{i}'})
                selected = 'unlabeled_clone' if i == 1 else rep
                split.append({**row, 'genome_id': selected, 'file_accession': selected, 'current_accession': selected,
                              'near_clone_group': f'N{i}', 'source_group': f'P{i}', 'split': ('train', 'validation', 'test')[i // 2]})
        for name, records in [('cohort.tsv', cohort), ('groups.tsv', groups), ('split.tsv', split)]:
            write_tsv(root / name, records, list(records[0]))
        with (root / 'labels.csv').open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['genome_name', 'drug', 'oneclass', 'conflict'])
            for i in range(7):
                writer.writerow([f'GCA_{i}.1', i % 2, '0' if i < 6 else '', '0' if i == 0 else ''])
            writer.writerow(['GCA_0.1', '', '', '1'])
            writer.writerow(['GCA_0.2', '1', '', ''])  # version-only, must not merge
        return cohort

    def run_audit(self, root):
        return run_amr_intersection(root / 'labels.csv', root / 'cohort.tsv', root / 'split.tsv', root / 'output',
                                    groups_path=root / 'groups.tsv', genome_root=root / 'genomes')

    def test_loss_trace_counts_conflicts_and_group_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture(root)
            result = self.run_audit(root)
            self.assertEqual(result['matched_assemblies'], 7)
            self.assertEqual(result['unique_known_biosamples'], 7)
            self.assertEqual(result['dispositions'], {'retained_split_representative': 4, 'near_clone_nonrepresentative': 1,
                                                      'identity_nonrepresentative': 1, 'cohort_quality_excluded': 1})
            self.assertEqual(result['version_only_label_accessions'], ['GCA_0.2'])
            self.assertEqual(result['legacy_loss_mechanisms']['legacy_identifier_key_mismatch'], 1)
            self.assertEqual(result['legacy_loss_mechanisms']['legacy_recognized_representative'], 3)
            drugs = {r['drug']: r for r in read_table(root / 'output/drug_feasibility.tsv')}
            self.assertEqual((drugs['drug']['S'], drugs['drug']['R']), ('3', '3'))
            self.assertEqual(drugs['drug']['inherited_splits_both_classes'], 'True')
            self.assertEqual(drugs['drug']['protocol_decision'], 'exploratory_only_protocol_pending')
            self.assertEqual(drugs['oneclass']['protocol_decision'], 'not_currently_evaluable_three_way')
            self.assertEqual(drugs['conflict']['conflicting_biosamples'], '1')
            self.assertEqual(drugs['conflict']['quality_eligible_biosamples'], '0')
            records = {r['genome_id']: r for r in read_table(root / 'output/record_disposition.tsv')}
            self.assertEqual(records['GCA_1.1']['selected_representative_has_label'], 'False')
            # No label is emitted for an unlabeled replacement assembly.
            cells = read_table(root / 'output/sample_drug_evidence.tsv')
            self.assertFalse(any('unlabeled' in r['source_genome_ids'] for r in cells))
            self.assertFalse(result['new_split_created'])
            self.assertFalse(result['training_performed'])

    def test_changed_content_unknown_identity_and_biosample_conflict(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cohort = self.fixture(root)
            (root / 'genomes/GCA_0.1.fna.gz').write_bytes(b'changed')
            cohort[1]['biosample'] = ''
            # Two assemblies on one BioSample disagree: do not double count or majority vote.
            cohort[3]['biosample'] = cohort[2]['biosample']
            write_tsv(root / 'cohort.tsv', cohort, list(cohort[0]))
            result = self.run_audit(root)
            self.assertEqual(result['current_file_status']['changed_since_snapshot'], 1)
            self.assertEqual(result['biosamples_with_multiple_assemblies'], 1)
            drugs = {r['drug']: r for r in read_table(root / 'output/drug_feasibility.tsv')}
            self.assertEqual(drugs['drug']['conflicting_biosamples'], '1')
            self.assertEqual(drugs['drug']['unknown_biosample_units'], '1')
            self.assertEqual((drugs['drug']['S'], drugs['drug']['R']), ('1', '1'))

    def test_malformed_labels_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'labels.csv'
            for data in ['genome_name,a,a\nX,0,1\n', 'genome_name,a\nX,0,1\n']:
                path.write_text(data)
                with self.assertRaises(ValueError):
                    _read_labels(path)

    def test_paired_only_mapping_and_invalid_label_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cohort = self.fixture(root)
            cohort[1]['paired_accession'] = 'GCF_1.1'
            write_tsv(root / 'cohort.tsv', cohort, list(cohort[0]))
            path = root / 'labels.csv'
            path.write_text(path.read_text().replace('GCA_1.1,1,0,', 'GCF_1.1,1,I,'))
            self.run_audit(root)
            drugs = {r['drug']: r for r in read_table(root / 'output/drug_feasibility.tsv')}
            self.assertEqual(drugs['drug']['R'], '3')
            self.assertEqual(drugs['drug']['direct_file_biosamples_R'], '2')
            self.assertEqual(drugs['drug']['alias_only_biosamples_R'], '1')
            self.assertEqual(drugs['oneclass']['invalid_label_units'], '1')

    def test_source_constraint_cannot_be_silently_relaxed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture(root)
            split = read_table(root / 'split.tsv')
            for row in split:
                row.update(source_group='single_source', split='train')
            write_tsv(root / 'split.tsv', split, list(split[0]))
            self.run_audit(root)
            drug = next(r for r in read_table(root / 'output/drug_feasibility.tsv') if r['drug'] == 'drug')
            self.assertEqual(drug['necessary_clone_group_support'], 'True')
            self.assertEqual(drug['necessary_source_group_support'], 'False')
            self.assertEqual(drug['protocol_decision'], 'not_currently_evaluable_three_way')
            self.assertEqual(drug['protocol_blocker'], 'insufficient_source_class_groups')


if __name__ == '__main__':
    unittest.main()
