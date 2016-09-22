import unittest
import tempfile
import os.path
import os

import pandas as pd
import skbio
import qiime
import numpy as np

from q2_demux._demux import BarcodeSequenceIterator
from q2_demux import emp, summary
from q2_types.per_sample_sequences import (
    FastqGzFormat, FastqManifestFormat, YamlFormat)


class BarcodeSequenceIteratorTests(unittest.TestCase):

    def test_valid(self):
        barcodes = [
            skbio.DNA('AAAA', metadata={'id': 's1/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AAAA', metadata={'id': 's2/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [25, 25, 25, 25]}),
            skbio.DNA('AACC', metadata={'id': 's3/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AACC', metadata={'id': 's4/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        sequences = [
            skbio.DNA('GGG', metadata={'id': 's1/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [22, 25, 22]}),
            skbio.DNA('CCC', metadata={'id': 's2/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [25, 25, 25]}),
            skbio.DNA('AAA', metadata={'id': 's3/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [22, 25, 22]}),
            skbio.DNA('TTT', metadata={'id': 's4/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [22, 25, 22]}),
        ]

        bsi = BarcodeSequenceIterator(barcodes, sequences)
        for i, (barcode, sequence) in enumerate(bsi):
            self.assertEqual(barcode, barcodes[i])
            self.assertEqual(sequence, sequences[i])

    def test_too_few_barcodes(self):
        barcodes = [
            skbio.DNA('AAAA', metadata={'id': 's1', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AAAA', metadata={'id': 's2', 'description': 'abc'},
                      positional_metadata={'quality': [25, 25, 25, 25]}),
            skbio.DNA('AACC', metadata={'id': 's3', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        sequences = [
            skbio.DNA('GGG', metadata={'id': 's1', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22]}),
            skbio.DNA('CCC', metadata={'id': 's2', 'description': 'abc'},
                      positional_metadata={'quality': [25, 25, 25]}),
            skbio.DNA('AAA', metadata={'id': 's3', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22]}),
            skbio.DNA('TTT', metadata={'id': 's4', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22]}),
        ]

        bsi = BarcodeSequenceIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_too_few_sequences(self):
        barcodes = [
            skbio.DNA('AAAA', metadata={'id': 's1', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AAAA', metadata={'id': 's2', 'description': 'abc'},
                      positional_metadata={'quality': [25, 25, 25, 25]}),
            skbio.DNA('AACC', metadata={'id': 's3', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AACC', metadata={'id': 's4', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        sequences = [
            skbio.DNA('GGG', metadata={'id': 's1', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22]}),
        ]

        bsi = BarcodeSequenceIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_id(self):
        barcodes = [
            skbio.DNA('AAAA', metadata={'id': 's1/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AAAA', metadata={'id': 's2/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [25, 25, 25, 25]}),
            skbio.DNA('AACC', metadata={'id': 's3/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AACC', metadata={'id': 's5/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        sequences = [
            skbio.DNA('GGG', metadata={'id': 's1/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [22, 25, 22]}),
            skbio.DNA('CCC', metadata={'id': 's2/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [25, 25, 25]}),
            skbio.DNA('AAA', metadata={'id': 's3/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [22, 25, 22]}),
            skbio.DNA('TTT', metadata={'id': 's4/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [22, 25, 22]}),
        ]

        bsi = BarcodeSequenceIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_description(self):
        barcodes = [
            skbio.DNA('AAAA', metadata={'id': 's1/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AAAA', metadata={'id': 's2/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [25, 25, 25, 25]}),
            skbio.DNA('AACC', metadata={'id': 's3/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AACC', metadata={'id': 's4/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        sequences = [
            skbio.DNA('GGG', metadata={'id': 's1/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [22, 25, 22]}),
            skbio.DNA('CCC', metadata={'id': 's2/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [25, 25, 25]}),
            skbio.DNA('AAA', metadata={'id': 's3/1', 'description': 'abc/1'},
                      positional_metadata={'quality': [22, 25, 22]}),
            skbio.DNA('TTT', metadata={'id': 's4/1', 'description': 'abd/1'},
                      positional_metadata={'quality': [22, 25, 22]}),
        ]

        bsi = BarcodeSequenceIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_handles_slashes_in_id(self):
        # mismatch is detected as being before the last slash, even if there
        # is more than one slash
        barcodes = [
            skbio.DNA('AAAA', metadata={'id': 's1/2/2', 'description': 'a/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        sequences = [
            skbio.DNA('GGG', metadata={'id': 's1/1/1', 'description': 'a/1'}),
        ]

        bsi = BarcodeSequenceIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_handles_slashes_in_description(self):
        # mismatch is detected as being before the last slash, even if there
        # is more than one slash
        barcodes = [
            skbio.DNA('AAAA', metadata={'id': 's1/2', 'description': 'a/2/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        sequences = [
            skbio.DNA('GGG', metadata={'id': 's1/1', 'description': 'a/1/1'}),
        ]

        bsi = BarcodeSequenceIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)


class EmpTests(unittest.TestCase):

    def setUp(self):
        barcodes = [
            skbio.DNA('AAAA', metadata={'id': 's1/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AAAA', metadata={'id': 's2/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [25, 25, 25, 25]}),
            skbio.DNA('AACC', metadata={'id': 's3/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AACC', metadata={'id': 's4/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        # setting the dtype on these quality scores as they end up being
        # compared to the output, which gets written as uint8
        self.sequences = [
            skbio.DNA('GGG', metadata={'id': 's1/1', 'description': 'abc/1'},
                      positional_metadata={'quality':
                      np.array([22, 25, 22], dtype='uint8')}),
            skbio.DNA('CCC', metadata={'id': 's2/1', 'description': 'abc/1'},
                      positional_metadata={'quality':
                      np.array([29, 29, 29], dtype='uint8')}),
            skbio.DNA('AAA', metadata={'id': 's3/1', 'description': 'abc/1'},
                      positional_metadata={'quality':
                      np.array([22, 20, 22], dtype='uint8')}),
            skbio.DNA('TTT', metadata={'id': 's4/1', 'description': 'abc/1'},
                      positional_metadata={'quality':
                      np.array([22, 25, 22], dtype='uint8')}),
        ]
        self.bsi = BarcodeSequenceIterator(barcodes, self.sequences)

        barcodes = pd.Series(['AAAA', 'AACC'], index=['sample1', 'sample2'])
        self.barcodes = qiime.MetadataCategory(barcodes)

    def _compare_manifests(self, act_manifest, exp_manifest):
        # strip comment lines before comparing
        act_manifest = [l for l in act_manifest if not l.startswith('#')]
        self.assertEqual(act_manifest, exp_manifest)

    def test_valid(self):
        actual = emp(self.bsi, self.barcodes)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # two per-sample files were written
        self.assertEqual(len(output_fastq), 2)

        # sequences in sample1 are correct
        sample1_seqs = skbio.io.read(output_fastq[0][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample1_seqs = list(sample1_seqs)
        self.assertEqual(len(sample1_seqs), 2)
        self.assertEqual(sample1_seqs[0], self.sequences[0])
        self.assertEqual(sample1_seqs[1], self.sequences[1])

        # sequences in sample2 are correct
        sample2_seqs = skbio.io.read(output_fastq[1][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample2_seqs = list(sample2_seqs)
        self.assertEqual(len(sample2_seqs), 2)
        self.assertEqual(sample2_seqs[0], self.sequences[2])
        self.assertEqual(sample2_seqs[1], self.sequences[3])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_2_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    def test_variable_length_barcodes(self):
        barcodes = pd.Series(['AAA', 'AACC'], index=['sample1', 'sample2'])
        barcodes = qiime.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp(self.bsi, barcodes)

    def test_duplicate_barcodes(self):
        barcodes = pd.Series(['AACC', 'AACC'], index=['sample1', 'sample2'])
        barcodes = qiime.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp(self.bsi, barcodes)

    def test_no_matched_barcodes(self):
        barcodes = pd.Series(['CCCC', 'GGCC'], index=['sample1', 'sample2'])
        barcodes = qiime.MetadataCategory(barcodes)
        with self.assertRaises(ValueError):
            emp(self.bsi, barcodes)

    def test_rev_comp_mapping_barcodes(self):
        barcodes = pd.Series(['TTTT', 'GGTT'], index=['sample1', 'sample2'])
        barcodes = qiime.MetadataCategory(barcodes)
        actual = emp(self.bsi, barcodes, rev_comp_mapping_barcodes=True)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # two per-sample files were written
        self.assertEqual(len(output_fastq), 2)

        # sequences in sample1 are correct
        sample1_seqs = skbio.io.read(output_fastq[0][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample1_seqs = list(sample1_seqs)
        self.assertEqual(len(sample1_seqs), 2)
        self.assertEqual(sample1_seqs[0], self.sequences[0])
        self.assertEqual(sample1_seqs[1], self.sequences[1])

        # sequences in sample2 are correct
        sample2_seqs = skbio.io.read(output_fastq[1][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample2_seqs = list(sample2_seqs)
        self.assertEqual(len(sample2_seqs), 2)
        self.assertEqual(sample2_seqs[0], self.sequences[2])
        self.assertEqual(sample2_seqs[1], self.sequences[3])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_2_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    def test_rev_comp_barcodes(self):
        barcodes = [
            skbio.DNA('TTTT', metadata={'id': 's1', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('TTTT', metadata={'id': 's2', 'description': 'abc'},
                      positional_metadata={'quality': [25, 25, 25, 25]}),
            skbio.DNA('GGTT', metadata={'id': 's3', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('GGTT', metadata={'id': 's4', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        bsi = BarcodeSequenceIterator(barcodes, self.sequences)
        actual = emp(bsi, self.barcodes, rev_comp_barcodes=True)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # two per-sample files were written
        self.assertEqual(len(output_fastq), 2)

        # sequences in sample1 are correct
        sample1_seqs = skbio.io.read(output_fastq[0][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample1_seqs = list(sample1_seqs)
        self.assertEqual(len(sample1_seqs), 2)
        self.assertEqual(sample1_seqs[0], self.sequences[0])
        self.assertEqual(sample1_seqs[1], self.sequences[1])

        # sequences in sample2 are correct
        sample2_seqs = skbio.io.read(output_fastq[1][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample2_seqs = list(sample2_seqs)
        self.assertEqual(len(sample2_seqs), 2)
        self.assertEqual(sample2_seqs[0], self.sequences[2])
        self.assertEqual(sample2_seqs[1], self.sequences[3])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_2_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)

    def test_barcode_trimming(self):
        # these barcodes are longer then the ones in the mapping file, so
        # only the first barcode_length bases should be read
        barcodes = [
            skbio.DNA('AAAAG', metadata={'id': 's1', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18, 29]}),
            skbio.DNA('AAAAG', metadata={'id': 's2', 'description': 'abc'},
                      positional_metadata={'quality': [25, 25, 25, 25, 29]}),
            skbio.DNA('AACCG', metadata={'id': 's3', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18, 29]}),
            skbio.DNA('AACCG', metadata={'id': 's4', 'description': 'abc'},
                      positional_metadata={'quality': [22, 25, 22, 18, 29]}),
        ]
        bsi = BarcodeSequenceIterator(barcodes, self.sequences)
        actual = emp(bsi, self.barcodes)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # two per-sample files were written
        self.assertEqual(len(output_fastq), 2)

        # sequences in sample1 are correct
        sample1_seqs = skbio.io.read(output_fastq[0][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample1_seqs = list(sample1_seqs)
        self.assertEqual(len(sample1_seqs), 2)
        self.assertEqual(sample1_seqs[0], self.sequences[0])
        self.assertEqual(sample1_seqs[1], self.sequences[1])

        # sequences in sample2 are correct
        sample2_seqs = skbio.io.read(output_fastq[1][1].open(),
                                     format='fastq',
                                     phred_offset=33, compression='gzip',
                                     constructor=skbio.DNA)
        sample2_seqs = list(sample2_seqs)
        self.assertEqual(len(sample2_seqs), 2)
        self.assertEqual(sample2_seqs[0], self.sequences[2])
        self.assertEqual(sample2_seqs[1], self.sequences[3])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_2_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

        # metadata is correct
        act_metadata = list(actual.metadata.view(YamlFormat).open())
        exp_metadata = ["{phred-offset: 33}\n"]
        self.assertEqual(act_metadata, exp_metadata)


class SummaryTests(unittest.TestCase):

    def setUp(self):
        barcodes = [
            skbio.DNA('AAAA', metadata={'id': 's1/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AAAA', metadata={'id': 's2/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [25, 25, 25, 25]}),
            skbio.DNA('AAAA', metadata={'id': 's3/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
            skbio.DNA('AACC', metadata={'id': 's4/2', 'description': 'abc/2'},
                      positional_metadata={'quality': [22, 25, 22, 18]}),
        ]
        sequences = [
            skbio.DNA('GGG', metadata={'id': 's1/1', 'description': 'abc/1'},
                      positional_metadata={'quality':
                      np.array([22, 25, 22])}),
            skbio.DNA('CCC', metadata={'id': 's2/1', 'description': 'abc/1'},
                      positional_metadata={'quality':
                      np.array([29, 29, 29])}),
            skbio.DNA('AAA', metadata={'id': 's3/1', 'description': 'abc/1'},
                      positional_metadata={'quality':
                      np.array([22, 20, 22])}),
            skbio.DNA('TTT', metadata={'id': 's4/1', 'description': 'abc/1'},
                      positional_metadata={'quality':
                      np.array([22, 25, 22])}),
        ]
        bsi = BarcodeSequenceIterator(barcodes, sequences)

        barcode_map = pd.Series(['AAAA', 'AACC'], index=['sample1', 'sample2'])
        barcode_map = qiime.MetadataCategory(barcode_map)

        self.demux_data = emp(bsi, barcode_map)

    def test_basic(self):
        # test that an index.html file is created and that it has size > 0
        with tempfile.TemporaryDirectory() as output_dir:
            result = summary(output_dir, self.demux_data)

            self.assertTrue(result is None)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.getsize(index_fp) > 0)
            csv_fp = os.path.join(output_dir, 'per-sample-fastq-counts.csv')
            self.assertTrue(os.path.exists(csv_fp))
            self.assertTrue(os.path.getsize(csv_fp) > 0)
            pdf_fp = os.path.join(output_dir, 'demultiplex-summary.pdf')
            self.assertTrue(os.path.exists(pdf_fp))
            self.assertTrue(os.path.getsize(pdf_fp) > 0)
            png_fp = os.path.join(output_dir, 'demultiplex-summary.png')
            self.assertTrue(os.path.exists(png_fp))
            self.assertTrue(os.path.getsize(png_fp) > 0)
            with open(index_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<td>Minimum:</td><td>1</td>', html)
                self.assertIn('<td>Maximum:</td><td>3</td>', html)
