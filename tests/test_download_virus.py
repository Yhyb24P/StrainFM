import gzip
import io
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from strainfm.audit import write_tsv
from strainfm.download import download_fasta, run_download
from strainfm.virus import split_virus
from strainfm.__main__ import main


class DownloadTest(unittest.TestCase):
    def test_batch_counts_and_cli_failure_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / 'manifest.tsv'
            write_tsv(manifest, [{'target_relpath': name, 'ncbi_url': 'https://example.org/a'} for name in ['a', 'b', 'c']], ['target_relpath', 'ncbi_url'])
            with patch('strainfm.download.download_fasta', side_effect=['OK', 'skip', 'FAIL c']):
                result = run_download(manifest, Path(tmp))
            self.assertEqual(result, {'records': 3, 'downloaded': 1, 'skipped_valid': 1, 'failed': 1})
            with patch('sys.argv', ['strainfm', 'download', '--manifest', str(manifest), '--root', tmp]), patch('strainfm.__main__.run_download', return_value=result), patch('builtins.print'):
                with self.assertRaises(SystemExit) as raised:
                    main()
                self.assertEqual(raised.exception.code, 1)

    def test_valid_skip_corrupt_retry_and_failed_replacement(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / 'a.fna.gz'
            target.write_bytes(b'not gzip')
            def fetch(command, **kwargs):
                Path(command[command.index('-o') + 1]).write_bytes(gzip.compress(b'>a\nACGT\n'))
            with patch('strainfm.download.subprocess.run', side_effect=fetch) as run:
                self.assertEqual(download_fasta(root, target.name, 'https://example.org/a'), 'OK')
                self.assertEqual(download_fasta(root, target.name, 'https://example.org/a'), 'skip')
                self.assertEqual(run.call_count, 1)
            for failure in [subprocess.CalledProcessError(1, 'curl'), None]:
                target.write_bytes(b'original corrupt file')
                with patch('strainfm.download.subprocess.run', side_effect=failure):
                    self.assertTrue(download_fasta(root, target.name, 'https://example.org/a').startswith('FAIL'))
                self.assertEqual(target.read_bytes(), b'original corrupt file')
                self.assertEqual(list(root.iterdir()), [target])

    def test_unsafe_paths_rejected_before_any_download(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'genomes'
            root.mkdir()
            (root / 'link').symlink_to(Path(tmp), target_is_directory=True)
            with patch('strainfm.download.subprocess.run') as run:
                for relative in ['../escaped.fna.gz', str(Path(tmp) / 'escaped.fna.gz'), 'link/escaped.fna.gz']:
                    with self.assertRaises(ValueError):
                        download_fasta(root, relative, 'https://example.org/a')
                manifest = Path(tmp) / 'manifest.tsv'
                write_tsv(manifest, [{'target_relpath': name, 'gwh_url': 'https://example.org/a'} for name in ['good.fna.gz', '../escaped.fna.gz']], ['target_relpath', 'gwh_url'])
                with self.assertRaises(ValueError):
                    run_download(manifest, root)
                run.assert_not_called()
            self.assertFalse((Path(tmp) / 'escaped.fna.gz').exists())


class FakeDecoder:
    def __init__(self, data, code=0):
        self.stdout = io.StringIO(data)
        self.code = code
    def __enter__(self):
        return self
    def __exit__(self, *args):
        self.stdout.close()
    def wait(self):
        return self.code
    def poll(self):
        return self.code


class VirusTest(unittest.TestCase):
    def test_failed_decoder_and_existing_output_preserve_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / 'manifest.tsv'
            write_tsv(manifest, [{'genome_id': 'UHGV-1', 'species_id': 'virus'}], ['genome_id', 'species_id'])
            output = root / 'genomes'
            with patch('strainfm.virus.subprocess.Popen', return_value=FakeDecoder('>UHGV-1\nACGT\n', 1)):
                with self.assertRaisesRegex(RuntimeError, 'decompression failed'):
                    split_virus(manifest, root / 'archive', output)
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob('.uhgv-*')), [])
            output.mkdir()
            original = output / 'original.gz'
            original.write_bytes(b'preserve')
            with patch('strainfm.virus.subprocess.Popen') as decoder:
                with self.assertRaises(FileExistsError):
                    split_virus(manifest, root / 'archive', output)
                decoder.assert_not_called()
            self.assertEqual(original.read_bytes(), b'preserve')

    def test_success_duplicate_and_unsafe_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / 'manifest.tsv'
            output = root / 'genomes'
            for species, data, expected in [('virus', '>UHGV-1\nACGT\n', None), ('virus', '>UHGV-1\nACGT\n>UHGV-1\nTTTT\n', ValueError), ('../escape', '>UHGV-1\nACGT\n', ValueError)]:
                write_tsv(manifest, [{'genome_id': 'UHGV-1', 'species_id': species}], ['genome_id', 'species_id'])
                destination = output if expected is None else root / 'failed'
                with patch('strainfm.virus.subprocess.Popen', return_value=FakeDecoder(data)):
                    if expected:
                        with self.assertRaises(expected):
                            split_virus(manifest, root / 'archive', destination)
                        self.assertFalse(destination.exists())
                    else:
                        self.assertEqual(split_virus(manifest, root / 'archive', destination)['records'], 1)
                        with gzip.open(output / 'virus/UHGV-1.fna.gz', 'rt') as handle:
                            self.assertEqual(handle.read(), data)
                self.assertEqual(list(root.glob('.uhgv-*')), [])


if __name__ == '__main__':
    unittest.main()
