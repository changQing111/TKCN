"""Run: python3 -m unittest discover -s test -p 'test_bowtie*.py' -v"""
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'src' / 'bowtie_shared_specific.py'
spec = importlib.util.spec_from_file_location('panel', SCRIPT)
panel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(panel)


class MappingTests(unittest.TestCase):
    def test_hit_classes_boundaries_strands_and_targets(self):
        records = {q: (t, 'ACGTA', '1') for q, t in [('a', 'T'), ('b', 'T'), ('c', 'T'), ('d', 'T'), ('e', 'U')]}
        regions = {('G', 'T', 'chr'): [(10, 20)], ('G', 'U', 'chr'): [(50, 60)]}
        def hit(q, pos, flag=0):
            return '{}\t{}\tchr\t{}\t255\t5M\t*\t0\t0\tACGTA\t*\tNM:i:0\n'.format(q, flag, pos)
        sam = '@SQ\tSN:chr\tLN:100\n' + hit('a', 11) + hit('a', 16, 16) + hit('b', 11) + hit('b', 17) + hit('c', 10) + 'd\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\t*\n' + hit('e', 11)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'a.sam'
            path.write_text(sam)
            result = panel.classify_sam(path, 'G', records, regions, Path(d) / 'hits.tsv')
            self.assertEqual(result, {'a': [2, 0], 'b': [1, 1], 'c': [0, 1], 'd': [0, 0], 'e': [0, 1]})
            path.write_text(sam.replace('5M', '4M1S'))
            with self.assertRaises(ValueError):
                panel.classify_sam(path, 'G', records, regions, Path(d) / 'hits.tsv')
            path.write_text(sam.replace(hit('e', 11), ''))
            with self.assertRaises(ValueError):
                panel.classify_sam(path, 'G', records, regions, Path(d) / 'hits.tsv')

    def test_coordinate_conversion_and_missing_regions(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'regions.tsv'
            p.write_text('target_id\tgenome_id\treference_name\tstart\tend\tcoordinate_system\nT\tA\tchr\t11\t20\t1-based-closed\nT\tB\tchr\t10\t20\t0-based-half-open\n')
            r = panel.load_regions(p, {'T'}, {'A', 'B'})
            self.assertEqual(r['A', 'T', 'chr'], r['B', 'T', 'chr'])
            with self.assertRaises(ValueError):
                panel.load_regions(p, {'T', 'U'}, {'A', 'B'})

    @unittest.skipUnless(shutil.which(os.environ.get('BOWTIE_BIN', 'bowtie')), 'Bowtie1 not available')
    def test_real_bowtie_multi_reference_pipeline(self):
        seqs = ['ACGTCAGTACGATCGTACCGATGCTAGTCCA', 'TGCATCGATAGGCTACGTCAGGATCAGTAC', 'GATCCGATGCTACCTGATCGTAGCATCGAC']
        seqs = [s[:29] for s in seqs]
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / 'shared.tsv').write_text(''.join(s + '\t2\n' for s in seqs))
            (root / 'targets.tsv').write_text('target_id\tshared_tsv\nT\tshared.tsv\n')
            (root / 'refs.tsv').write_text('genome_id\tfasta\tindex_prefix\nA\tA.fa\t\nB\tB.fa\t\n')
            # First query has two target copies, one reverse-complemented.
            rc = seqs[0].translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            for genome in ['A', 'B']:
                (root / (genome + '.fa')).write_text('>target\n' + seqs[0] + 'N'*40 + rc + 'N'*40 + seqs[1] + 'N'*40 + (seqs[2] if genome == 'A' else 'N'*len(seqs[2])) + '\n>off\n' + seqs[1] + '\n')
            (root / 'regions.tsv').write_text('target_id\tgenome_id\treference_name\tstart\tend\tcoordinate_system\n' + ''.join('T\t{}\ttarget\t1\t{}\t1-based-closed\n'.format(g, 4*len(seqs[0])+120) for g in ['A', 'B']))
            cmd = [sys.executable, str(SCRIPT), '--targets', str(root/'targets.tsv'), '--references', str(root/'refs.tsv'), '--regions', str(root/'regions.tsv'), '--out-dir', str(root/'out'), '-t', '1', '--keep-sam']
            subprocess.run(cmd, check=True)
            self.assertEqual((root/'out/T.shared-specific.tsv').read_text(), seqs[0] + '\t2\n')
            self.assertTrue((root/'out/SUCCESS').is_file())
            classified = (root/'out/classification.tsv').read_text()
            self.assertIn('mixed', classified)
            self.assertIn('unmapped', classified)
            self.assertIn('\t2\t0\ttarget_only\t1', classified)
            subprocess.run([sys.executable, str(SCRIPT.parent / 'filter_target_specific_kmers.py'), str(root/'out/T.shared-specific.tsv'), '--out-prefix', str(root/'qc')], check=True, stdout=subprocess.PIPE)
            self.assertTrue((root/'qc.fa').is_file())
            self.assertNotEqual(subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).returncode, 0)


if __name__ == '__main__':
    unittest.main()
