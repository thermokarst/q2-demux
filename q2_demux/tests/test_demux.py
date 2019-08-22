# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import unittest.mock as mock
import os.path
import tempfile
import json
import shutil
import random

import pandas as pd
import pandas.testing as pdt
import skbio
import qiime2
import numpy.testing as npt

from qiime2.plugin.testing import TestPluginBase
from q2_demux._demux import (BarcodeSequenceFastqIterator,
                             BarcodePairedSequenceFastqIterator)
from q2_demux import emp_single, emp_paired, summarize
from q2_types.per_sample_sequences import (
    FastqGzFormat, FastqManifestFormat,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_demux._summarize._visualizer import (_PlotQualView,
                                             _decode_qual_to_phred33)


class BarcodeSequenceFastqIteratorTests(unittest.TestCase):

    def test_valid(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        for i, (barcode, sequence) in enumerate(bsi):
            self.assertEqual(barcode, barcodes[i])
            self.assertEqual(sequence, sequences[i])

    def test_too_few_barcodes(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_too_few_sequences(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_id(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s5/1 abc/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_description(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abd/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_handles_slashes_in_id(self):
        # mismatch is detected as being before the last slash, even if there
        # is more than one slash
        barcodes = [('@s1/2/2 abc/2', 'AAAA', '+', 'YYYY')]
        sequences = [('@s1/1/1 abc/1', 'GGG', '+', 'YYY')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_mismatched_handles_slashes_in_description(self):
        # mismatch is detected as being before the last slash, even if there
        # is more than one slash
        barcodes = [('@s1/2 a/2/2', 'AAAA', '+', 'YYYY')]
        sequences = [('@s1/1 a/1/1', 'GGG', '+', 'YYY')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

    def test_no_description(self):
        barcodes = [('@s1/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1', 'GGG', '+', 'YYY'),
                     ('@s2/1', 'CCC', '+', 'PPP'),
                     ('@s3/1', 'AAA', '+', 'PPP'),
                     ('@s4/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        for i, (barcode, sequence) in enumerate(bsi):
            self.assertEqual(barcode, barcodes[i])
            self.assertEqual(sequence, sequences[i])

    def test_only_one_description(self):
        barcodes = [('@s1/2 abc', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc', 'AAAA', '+', 'PPPP'),
                    ('@s3/2 abc', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1', 'GGG', '+', 'YYY'),
                     ('@s2/1', 'CCC', '+', 'PPP'),
                     ('@s3/1', 'AAA', '+', 'PPP'),
                     ('@s4/1', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)

        barcodes = [('@s1/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2', 'AAAA', '+', 'PPPP'),
                    ('@s3/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2', 'AACC', '+', 'PPPP')]

        sequences = [('@s1/1 abc', 'GGG', '+', 'YYY'),
                     ('@s2/1 abc', 'CCC', '+', 'PPP'),
                     ('@s3/1 abc', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc', 'TTT', '+', 'PPP')]

        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)
        with self.assertRaises(ValueError):
            list(bsi)


class EmpTestingUtils:
    def _compare_sequence_to_record(self, sequence, fields):
        header_line = ' '.join([sequence.metadata['id'],
                                sequence.metadata['description']])
        self.assertEqual(fields[0][1:], header_line)
        self.assertEqual(fields[1], str(sequence))
        npt.assert_array_equal(_decode_qual_to_phred33(fields[3]),
                               sequence.positional_metadata['quality'])

    def _compare_manifests(self, act_manifest, exp_manifest):
        # strip comment lines before comparing
        act_manifest = [l for l in act_manifest if not l.startswith('#')]
        self.assertEqual(act_manifest, exp_manifest)

    def _validate_sample_fastq(self, fastq, sequences, indices):
        seqs = skbio.io.read(fastq, format='fastq', phred_offset=33,
                             compression='gzip', constructor=skbio.DNA)
        seqs = list(seqs)
        self.assertEqual(len(seqs), len(indices))
        for idx, i in enumerate(indices):
            self._compare_sequence_to_record(seqs[idx], sequences[i])


class EmpSingleTests(unittest.TestCase, EmpTestingUtils):

    def setUp(self):
        barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'AACC', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'AAAA', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'CGGC', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'GGAA', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'CGGC', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'CGGC', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'GGAA', '+', 'PPPP')]

        golaybarcodes = [  # ATGATGCGACCA -> ACGATGCGACCA
                         ('@s1/2 abc/2', 'ATGATGCGACCA', '+', 'YYYYYYYYYYYY'),
                         ('@s2/2 abc/2', 'AGCTATCCACGA', '+', 'PPPPPPPPPPPP'),
                         ('@s3/2 abc/2', 'ACACACTATGGC', '+', 'PPPPPPPPPPPP'),
                         ('@s4/2 abc/2', 'AGCTATCCACGA', '+', 'PPPPPPPPPPPP'),
                         ('@s5/2 abc/2', 'ACACACTATGGC', '+', 'PPPPPPPPPPPP'),
                         ('@s6/2 abc/2', 'ACGATGCGACCA', '+', 'PPPPPPPPPPPP'),
                         # CATTGTATCAAC -> CATCGTATCAAC
                         ('@s7/2 abc/2', 'CATTGTATCAAC', '+', 'PPPPPPPPPPPP'),
                         ('@s8/2 abc/2', 'CTAACGCAGGGG', '+', 'PPPPPPPPPPPP'),
                         ('@s9/2 abc/2', 'CATCGTATCAAC', '+', 'PPPPPPPPPPPP'),
                         ('@s10/2 abc/2', 'CATCGTATCAAC', '+', 'PPPPPPPPPPPP'),
                         ('@s11/2 abc/2', 'CTAACGCAGTCA', '+', 'PPPPPPPPPPPP')]

        self.sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                          ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                          ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                          ('@s4/1 abc/1', 'TTT', '+', 'PPP'),
                          ('@s5/1 abc/1', 'ATA', '+', 'PPP'),
                          ('@s6/1 abc/1', 'TAT', '+', 'PPP'),
                          ('@s7/1 abc/1', 'CGC', '+', 'PPP'),
                          ('@s8/1 abc/1', 'GCG', '+', 'PPP'),
                          ('@s9/1 abc/1', 'ACG', '+', 'PPP'),
                          ('@s10/1 abc/1', 'GCA', '+', 'PPP'),
                          ('@s11/1 abc/1', 'TGA', '+', 'PPP')]
        self.bsi = BarcodeSequenceFastqIterator(barcodes, self.sequences)
        barcode_map = pd.Series(
            ['AAAA', 'AACC', 'TTAA', 'GGAA', 'CGGC'], name='bc',
            index=pd.Index(['sample1', 'sample2', 'sample3', 'sample4',
                            'sample5'], name='id')
        )
        self.barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        self.bsi_werr = BarcodeSequenceFastqIterator(golaybarcodes,
                                                     self.sequences)

        golay_barcode_map = pd.Series(
            ['ACGATGCGACCA', 'ACACACTATGGC', 'AGCTATCCACGA',
             'CTAACGCAGTCA', 'CATCGTATCAAC'], name='bc',
            index=pd.Index(['sample1', 'sample2', 'sample3', 'sample4',
                            'sample5'], name='id')
        )
        self.golay_barcode_map = qiime2.CategoricalMetadataColumn(
            golay_barcode_map)

    def test_valid(self):
        actual, error_detail = emp_single(self.bsi, self.barcode_map,
                                          golay_error_correction=False)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # five per-sample files were written
        self.assertEqual(len(output_fastq), 5)

        # sequences in sample1 are correct
        self._validate_sample_fastq(
            output_fastq[0][1].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            output_fastq[1][1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            output_fastq[2][1].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            output_fastq[3][1].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            output_fastq[4][1].open(), self.sequences, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

    def test_valid_ecc_no_golay(self):
        _, ecc = emp_single(self.bsi, self.barcode_map,
                            golay_error_correction=False)
        exp_errors = pd.DataFrame([
            ['sample1', '@s1/2 abc/2',  'AAAA', None, None],
            ['sample3', '@s2/2 abc/2',  'TTAA', None, None],
            ['sample2', '@s3/2 abc/2',  'AACC', None, None],
            ['sample3', '@s4/2 abc/2',  'TTAA', None, None],
            ['sample2', '@s5/2 abc/2',  'AACC', None, None],
            ['sample1', '@s6/2 abc/2',  'AAAA', None, None],
            ['sample5', '@s7/2 abc/2',  'CGGC', None, None],
            ['sample4', '@s8/2 abc/2',  'GGAA', None, None],
            ['sample5', '@s9/2 abc/2',  'CGGC', None, None],
            ['sample5', '@s10/2 abc/2', 'CGGC', None, None],
            ['sample4', '@s11/2 abc/2', 'GGAA', None, None]
            ],
            columns=['sample', 'barcode-sequence-id',
                     'barcode-uncorrected', 'barcode-corrected',
                     'barcode-errors'],
            index=pd.Index(['record-%02d' % i for i in range(1, 12)],
                           name='id'))
        pdt.assert_frame_equal(ecc, exp_errors)

    def test_valid_with_barcode_errors(self):
        actual, error_detail = emp_single(self.bsi_werr,
                                          self.golay_barcode_map,
                                          golay_error_correction=True)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # five per-sample files were written
        self.assertEqual(len(output_fastq), 5)

        # sequences in sample1 are correct
        self._validate_sample_fastq(
            output_fastq[0][1].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            output_fastq[1][1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            output_fastq[2][1].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            output_fastq[3][1].open(), self.sequences, [10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            output_fastq[4][1].open(), self.sequences, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)
        exp_errors = pd.DataFrame([
            ['sample1', '@s1/2 abc/2',  'ATGATGCGACCA', 'ACGATGCGACCA', 1],
            ['sample3', '@s2/2 abc/2',  'AGCTATCCACGA', 'AGCTATCCACGA', 0],
            ['sample2', '@s3/2 abc/2',  'ACACACTATGGC', 'ACACACTATGGC', 0],
            ['sample3', '@s4/2 abc/2',  'AGCTATCCACGA', 'AGCTATCCACGA', 0],
            ['sample2', '@s5/2 abc/2',  'ACACACTATGGC', 'ACACACTATGGC', 0],
            ['sample1', '@s6/2 abc/2',  'ACGATGCGACCA', 'ACGATGCGACCA', 0],
            ['sample5', '@s7/2 abc/2',  'CATTGTATCAAC', 'CATCGTATCAAC', 1],
            [None,      '@s8/2 abc/2',  'CTAACGCAGGGG', None,           4],
            ['sample5', '@s9/2 abc/2',  'CATCGTATCAAC', 'CATCGTATCAAC', 0],
            ['sample5', '@s10/2 abc/2', 'CATCGTATCAAC', 'CATCGTATCAAC', 0],
            ['sample4', '@s11/2 abc/2', 'CTAACGCAGTCA', 'CTAACGCAGTCA', 0]
            ],
            columns=['sample', 'barcode-sequence-id',
                     'barcode-uncorrected', 'barcode-corrected',
                     'barcode-errors'],
            index=pd.Index(['record-%02d' % i for i in range(1, 12)],
                           name='id'))
        pdt.assert_frame_equal(error_detail, exp_errors)

    @mock.patch('q2_demux._demux.OPEN_FH_LIMIT', 3)
    def test_valid_small_open_fh_limit(self):
        self.test_valid()

    def test_variable_length_barcodes(self):
        barcodes = pd.Series(['AAA', 'AACC'], name='bc',
                             index=pd.Index(['sample1', 'sample2'], name='id'))
        barcodes = qiime2.CategoricalMetadataColumn(barcodes)
        with self.assertRaises(ValueError):
            emp_single(self.bsi, barcodes, golay_error_correction=False)

    def test_duplicate_barcodes(self):
        barcodes = pd.Series(['AACC', 'AACC'], name='bc',
                             index=pd.Index(['sample1', 'sample2'], name='id'))
        barcodes = qiime2.CategoricalMetadataColumn(barcodes)
        with self.assertRaises(ValueError):
            emp_single(self.bsi, barcodes, golay_error_correction=False)

    def test_no_matched_barcodes(self):
        barcodes = pd.Series(['CCCC', 'GGCC'], name='bc',
                             index=pd.Index(['sample1', 'sample2'], name='id'))
        barcodes = qiime2.CategoricalMetadataColumn(barcodes)
        with self.assertRaises(ValueError):
            emp_single(self.bsi, barcodes, golay_error_correction=False)

    def test_rev_comp_mapping_barcodes(self):
        barcodes = pd.Series(
            ['TTTT', 'GGTT', 'TTAA', 'TTCC', 'GCCG'], name='bc',
            index=pd.Index(['sample1', 'sample2', 'sample3', 'sample4',
                            'sample5'], name='id')
        )
        barcodes = qiime2.CategoricalMetadataColumn(barcodes)
        actual, ecc = emp_single(self.bsi, barcodes,
                                 rev_comp_mapping_barcodes=True,
                                 golay_error_correction=False)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # five per-sample files were written
        self.assertEqual(len(output_fastq), 5)

        # sequences in sample1 are correct
        self._validate_sample_fastq(
            output_fastq[0][1].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            output_fastq[1][1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            output_fastq[2][1].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            output_fastq[3][1].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            output_fastq[4][1].open(), self.sequences, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

    def test_rev_comp_barcodes(self):
        barcodes = [('@s1/2 abc/2', 'TTTT', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'GGTT', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'GGTT', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'TTTT', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'TTCC', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'TTCC', '+', 'PPPP')]
        bsi = BarcodeSequenceFastqIterator(barcodes, self.sequences)
        actual, ecc = emp_single(bsi, self.barcode_map,
                                 golay_error_correction=False,
                                 rev_comp_barcodes=True)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # five per-sample files were written
        self.assertEqual(len(output_fastq), 5)

        # sequences in sample1 are correct
        self._validate_sample_fastq(
            output_fastq[0][1].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            output_fastq[1][1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            output_fastq[2][1].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            output_fastq[3][1].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            output_fastq[4][1].open(), self.sequences, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)

    def test_barcode_trimming(self):
        # these barcodes are longer then the ones in the mapping file, so
        # only the first barcode_length bases should be read
        barcodes = [('@s1/2 abc/2', 'AAAAG', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAAG', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACCG', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAAG', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'AACCG', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'AAAAG', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'GGAAG', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'GGAAG', '+', 'PPPP')]
        bsi = BarcodeSequenceFastqIterator(barcodes, self.sequences)
        actual, ecc = emp_single(bsi, self.barcode_map,
                                 golay_error_correction=False)
        output_fastq = list(actual.sequences.iter_views(FastqGzFormat))
        # five per-sample files were written
        self.assertEqual(len(output_fastq), 5)

        # sequences in sample1 are correct
        self._validate_sample_fastq(
            output_fastq[0][1].open(), self.sequences, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            output_fastq[1][1].open(), self.sequences, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            output_fastq[2][1].open(), self.sequences, [1, 3])

        # sequences in sample4 are correct
        self._validate_sample_fastq(
            output_fastq[3][1].open(), self.sequences, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            output_fastq[4][1].open(), self.sequences, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n']
        self._compare_manifests(act_manifest, exp_manifest)


class EmpPairedTests(unittest.TestCase, EmpTestingUtils):
    def setUp(self):
        self.barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                         ('@s2/2 abc/2', 'TTAA', '+', 'PPPP'),
                         ('@s3/2 abc/2', 'AACC', '+', 'PPPP'),
                         ('@s4/2 abc/2', 'TTAA', '+', 'PPPP'),
                         ('@s5/2 abc/2', 'AACC', '+', 'PPPP'),
                         ('@s6/2 abc/2', 'AAAA', '+', 'PPPP'),
                         ('@s7/2 abc/2', 'CGGC', '+', 'PPPP'),
                         ('@s8/2 abc/2', 'GGAA', '+', 'PPPP'),
                         ('@s9/2 abc/2', 'CGGC', '+', 'PPPP'),
                         ('@s10/2 abc/2', 'CGGC', '+', 'PPPP'),
                         ('@s11/2 abc/2', 'GGAA', '+', 'PPPP')]

        golaybarcodes = [  # ATGATGCGACCA -> ACGATGCGACCA
                         ('@s1/2 abc/2', 'ATGATGCGACCA', '+', 'YYYYYYYYYYYY'),
                         ('@s2/2 abc/2', 'AGCTATCCACGA', '+', 'PPPPPPPPPPPP'),
                         ('@s3/2 abc/2', 'ACACACTATGGC', '+', 'PPPPPPPPPPPP'),
                         ('@s4/2 abc/2', 'AGCTATCCACGA', '+', 'PPPPPPPPPPPP'),
                         ('@s5/2 abc/2', 'ACACACTATGGC', '+', 'PPPPPPPPPPPP'),
                         ('@s6/2 abc/2', 'ACGATGCGACCA', '+', 'PPPPPPPPPPPP'),
                         # CATTGTATCAAC -> CATCGTATCAAC
                         ('@s7/2 abc/2', 'CATTGTATCAAC', '+', 'PPPPPPPPPPPP'),
                         ('@s8/2 abc/2', 'CTAACGCAGGGG', '+', 'PPPPPPPPPPPP'),
                         ('@s9/2 abc/2', 'CATCGTATCAAC', '+', 'PPPPPPPPPPPP'),
                         ('@s10/2 abc/2', 'CATCGTATCAAC', '+', 'PPPPPPPPPPPP'),
                         ('@s11/2 abc/2', 'CTAACGCAGTCA', '+', 'PPPPPPPPPPPP')]

        self.forward = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                        ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                        ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                        ('@s4/1 abc/1', 'TTT', '+', 'PPP'),
                        ('@s5/1 abc/1', 'ATA', '+', 'PPP'),
                        ('@s6/1 abc/1', 'TAT', '+', 'PPP'),
                        ('@s7/1 abc/1', 'CGC', '+', 'PPP'),
                        ('@s8/1 abc/1', 'GCG', '+', 'PPP'),
                        ('@s9/1 abc/1', 'ACG', '+', 'PPP'),
                        ('@s10/1 abc/1', 'GCA', '+', 'PPP'),
                        ('@s11/1 abc/1', 'TGA', '+', 'PPP')]

        self.reverse = [('@s1/1 abc/1', 'CCC', '+', 'YYY'),
                        ('@s2/1 abc/1', 'GGG', '+', 'PPP'),
                        ('@s3/1 abc/1', 'TTT', '+', 'PPP'),
                        ('@s4/1 abc/1', 'AAA', '+', 'PPP'),
                        ('@s5/1 abc/1', 'TAT', '+', 'PPP'),
                        ('@s6/1 abc/1', 'ATA', '+', 'PPP'),
                        ('@s7/1 abc/1', 'GCG', '+', 'PPP'),
                        ('@s8/1 abc/1', 'CGC', '+', 'PPP'),
                        ('@s9/1 abc/1', 'CGT', '+', 'PPP'),
                        ('@s10/1 abc/1', 'TGC', '+', 'PPP'),
                        ('@s11/1 abc/1', 'TCA', '+', 'PPP')]

        self.bpsi = BarcodePairedSequenceFastqIterator(
            self.barcodes, self.forward, self.reverse)

        barcode_map = pd.Series(
            ['AAAA', 'AACC', 'TTAA', 'GGAA', 'CGGC'], name='bc',
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')
        )
        self.barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        self.bpsi_werr = BarcodePairedSequenceFastqIterator(golaybarcodes,
                                                            self.forward,
                                                            self.reverse)

        golay_barcode_map = pd.Series(
            ['ACGATGCGACCA', 'ACACACTATGGC', 'AGCTATCCACGA',
             'CTAACGCAGTCA', 'CATCGTATCAAC'], name='bc',
            index=pd.Index(['sample1', 'sample2', 'sample3', 'sample4',
                            'sample5'], name='id')
        )
        self.golay_barcode_map = qiime2.CategoricalMetadataColumn(
            golay_barcode_map)

    def check_valid(self, *args, **kwargs):
        actual, ecc = emp_paired(*args, **kwargs)

        # five forward sample files
        forward_fastq = [
            view for path, view in actual.sequences.iter_views(FastqGzFormat)
            if 'R1_001.fastq' in path.name]
        self.assertEqual(len(forward_fastq), 5)

        # five reverse sample files
        reverse_fastq = [
            view for path, view in actual.sequences.iter_views(FastqGzFormat)
            if 'R2_001.fastq' in path.name]
        self.assertEqual(len(reverse_fastq), 5)

        # FORWARD:
        # sequences in sample1 are correct
        self._validate_sample_fastq(
            forward_fastq[0].open(), self.forward, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            forward_fastq[1].open(), self.forward, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            forward_fastq[2].open(), self.forward, [1, 3])

        # sequences in sample4 are correct
        if kwargs['golay_error_correction']:
            self._validate_sample_fastq(
                forward_fastq[3].open(), self.forward, [10, ])
        else:
            self._validate_sample_fastq(
                forward_fastq[3].open(), self.forward, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            forward_fastq[4].open(), self.forward, [6, 8, 9])

        # REVERSE:
        # sequences in sample1 are correct
        self._validate_sample_fastq(
            reverse_fastq[0].open(), self.reverse, [0, 5])

        # sequences in sample2 are correct
        self._validate_sample_fastq(
            reverse_fastq[1].open(), self.reverse, [2, 4])

        # sequences in sample3 are correct
        self._validate_sample_fastq(
            reverse_fastq[2].open(), self.reverse, [1, 3])

        # sequences in sample4 are correct
        if kwargs['golay_error_correction']:
            self._validate_sample_fastq(
                reverse_fastq[3].open(), self.reverse, [10, ])
        else:
            self._validate_sample_fastq(
                reverse_fastq[3].open(), self.reverse, [7, 10])

        # sequences in sample5 are correct
        self._validate_sample_fastq(
            reverse_fastq[4].open(), self.reverse, [6, 8, 9])

        # manifest is correct
        act_manifest = list(actual.manifest.view(FastqManifestFormat).open())
        exp_manifest = ['sample-id,filename,direction\n',
                        'sample1,sample1_1_L001_R1_001.fastq.gz,forward\n',
                        'sample1,sample1_1_L001_R2_001.fastq.gz,reverse\n',
                        'sample3,sample3_2_L001_R1_001.fastq.gz,forward\n',
                        'sample3,sample3_2_L001_R2_001.fastq.gz,reverse\n',
                        'sample2,sample2_3_L001_R1_001.fastq.gz,forward\n',
                        'sample2,sample2_3_L001_R2_001.fastq.gz,reverse\n',
                        'sample5,sample5_4_L001_R1_001.fastq.gz,forward\n',
                        'sample5,sample5_4_L001_R2_001.fastq.gz,reverse\n',
                        'sample4,sample4_5_L001_R1_001.fastq.gz,forward\n',
                        'sample4,sample4_5_L001_R2_001.fastq.gz,reverse\n']

        self._compare_manifests(act_manifest, exp_manifest)

        if kwargs['golay_error_correction']:
            exp_errors = pd.DataFrame([
                ['sample1', '@s1/2 abc/2',  'ATGATGCGACCA', 'ACGATGCGACCA', 1],
                ['sample3', '@s2/2 abc/2',  'AGCTATCCACGA', 'AGCTATCCACGA', 0],
                ['sample2', '@s3/2 abc/2',  'ACACACTATGGC', 'ACACACTATGGC', 0],
                ['sample3', '@s4/2 abc/2',  'AGCTATCCACGA', 'AGCTATCCACGA', 0],
                ['sample2', '@s5/2 abc/2',  'ACACACTATGGC', 'ACACACTATGGC', 0],
                ['sample1', '@s6/2 abc/2',  'ACGATGCGACCA', 'ACGATGCGACCA', 0],
                ['sample5', '@s7/2 abc/2',  'CATTGTATCAAC', 'CATCGTATCAAC', 1],
                [None,      '@s8/2 abc/2',  'CTAACGCAGGGG', None,           4],
                ['sample5', '@s9/2 abc/2',  'CATCGTATCAAC', 'CATCGTATCAAC', 0],
                ['sample5', '@s10/2 abc/2', 'CATCGTATCAAC', 'CATCGTATCAAC', 0],
                ['sample4', '@s11/2 abc/2', 'CTAACGCAGTCA', 'CTAACGCAGTCA', 0]
                ],
                columns=['sample', 'barcode-sequence-id',
                         'barcode-uncorrected', 'barcode-corrected',
                         'barcode-errors'],
                index=pd.Index(['record-%02d' % i for i in range(1, 12)],
                               name='id'))
            pdt.assert_frame_equal(ecc, exp_errors)

    def test_valid(self):
        self.check_valid(self.bpsi, self.barcode_map,
                         golay_error_correction=False)

    def test_valid_ecc_no_golay(self):
        _, ecc = emp_paired(self.bpsi, self.barcode_map,
                            golay_error_correction=False)
        exp_errors = pd.DataFrame([
            ['sample1', '@s1/2 abc/2',  'AAAA', None, None],
            ['sample3', '@s2/2 abc/2',  'TTAA', None, None],
            ['sample2', '@s3/2 abc/2',  'AACC', None, None],
            ['sample3', '@s4/2 abc/2',  'TTAA', None, None],
            ['sample2', '@s5/2 abc/2',  'AACC', None, None],
            ['sample1', '@s6/2 abc/2',  'AAAA', None, None],
            ['sample5', '@s7/2 abc/2',  'CGGC', None, None],
            ['sample4', '@s8/2 abc/2',  'GGAA', None, None],
            ['sample5', '@s9/2 abc/2',  'CGGC', None, None],
            ['sample5', '@s10/2 abc/2', 'CGGC', None, None],
            ['sample4', '@s11/2 abc/2', 'GGAA', None, None]
            ],
            columns=['sample', 'barcode-sequence-id',
                     'barcode-uncorrected', 'barcode-corrected',
                     'barcode-errors'],
            index=pd.Index(['record-%02d' % i for i in range(1, 12)],
                           name='id'))
        pdt.assert_frame_equal(ecc, exp_errors)

    def test_valid_with_barcode_errors(self):
        self.check_valid(self.bpsi_werr, self.golay_barcode_map,
                         golay_error_correction=True)

    @mock.patch('q2_demux._demux.OPEN_FH_LIMIT', 6)
    def test_valid_small_open_fh_limit(self):
        self.test_valid()

    def test_variable_length_barcodes(self):
        barcodes = pd.Series(
            ['AAA', 'AACC'], name='bc',
            index=pd.Index(['sample1', 'sample2'], name='id')
        )
        barcodes = qiime2.CategoricalMetadataColumn(barcodes)
        with self.assertRaises(ValueError):
            emp_paired(self.bpsi, barcodes, golay_error_correction=False)

    def test_duplicate_barcodes(self):
        barcodes = pd.Series(
            ['AACC', 'AACC'], name='bc',
            index=pd.Index(['sample1', 'sample2'], name='id')
        )
        barcodes = qiime2.CategoricalMetadataColumn(barcodes)
        with self.assertRaises(ValueError):
            emp_paired(self.bpsi, barcodes, golay_error_correction=False)

    def test_no_matched_barcodes(self):
        barcodes = pd.Series(
            ['CCCC', 'GGCC'], name='bc',
            index=pd.Index(['sample1', 'sample2'], name='id')
        )
        barcodes = qiime2.CategoricalMetadataColumn(barcodes)
        with self.assertRaises(ValueError):
            emp_paired(self.bpsi, barcodes, golay_error_correction=False)

    def test_rev_comp_mapping_barcodes(self):
        barcodes = pd.Series(
            ['TTTT', 'GGTT', 'TTAA', 'TTCC', 'GCCG'], name='bc',
            index=pd.Index(['sample1', 'sample2', 'sample3', 'sample4',
                            'sample5'], name='id')
        )
        barcodes = qiime2.CategoricalMetadataColumn(barcodes)
        self.check_valid(self.bpsi, barcodes, rev_comp_mapping_barcodes=True,
                         golay_error_correction=False)

    def test_rev_comp_barcodes(self):
        barcodes = [('@s1/2 abc/2', 'TTTT', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'GGTT', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAA', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'GGTT', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'TTTT', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'TTCC', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'GCCG', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'TTCC', '+', 'PPPP')]
        bpsi = BarcodePairedSequenceFastqIterator(
            barcodes, self.forward, self.reverse)
        self.check_valid(bpsi, self.barcode_map, rev_comp_barcodes=True,
                         golay_error_correction=False)

    def test_barcode_trimming(self):
        # these barcodes are longer then the ones in the mapping file, so
        # only the first barcode_length bases should be read
        barcodes = [('@s1/2 abc/2', 'AAAAG', '+', 'YYYY'),
                    ('@s2/2 abc/2', 'TTAAG', '+', 'PPPP'),
                    ('@s3/2 abc/2', 'AACCG', '+', 'PPPP'),
                    ('@s4/2 abc/2', 'TTAAG', '+', 'PPPP'),
                    ('@s5/2 abc/2', 'AACCG', '+', 'PPPP'),
                    ('@s6/2 abc/2', 'AAAAG', '+', 'PPPP'),
                    ('@s7/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s8/2 abc/2', 'GGAAG', '+', 'PPPP'),
                    ('@s9/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s10/2 abc/2', 'CGGCG', '+', 'PPPP'),
                    ('@s11/2 abc/2', 'GGAAG', '+', 'PPPP')]
        bpsi = BarcodePairedSequenceFastqIterator(
            barcodes, self.forward, self.reverse)
        self.check_valid(bpsi, self.barcode_map, golay_error_correction=False)


class SummarizeTests(TestPluginBase):
    package = 'q2_demux.tests'

    def setUp(self):
        super().setUp()
        self.barcodes = [('@s1/2 abc/2', 'AAAA', '+', 'YYYY'),
                         ('@s2/2 abc/2', 'AAAA', '+', 'PPPP'),
                         ('@s3/2 abc/2', 'AAAA', '+', 'PPPP'),
                         ('@s4/2 abc/2', 'AACC', '+', 'PPPP')]

        self.sequences = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                          ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                          ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                          ('@s4/1 abc/1', 'TTT', '+', 'PPP')]

    def test_basic(self):
        bsi = BarcodeSequenceFastqIterator(self.barcodes, self.sequences)

        barcode_map = pd.Series(
            ['AAAA', 'AACC'], name='bc',
            index=pd.Index(['sample_1', 'sample2'], name='id')
        )
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        demux_data, ecc = emp_single(bsi, barcode_map,
                                     golay_error_correction=False)
        # test that an index.html file is created and that it has size > 0
        with tempfile.TemporaryDirectory() as output_dir:
            # TODO: Remove _PlotQualView wrapper
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=False), n=2)
            self.assertTrue(result is None)
            index_fp = os.path.join(output_dir, 'overview.html')
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
            qual_forward_fp = os.path.join(
                output_dir, 'forward-seven-number-summaries.csv')
            self.assertTrue(os.path.exists(qual_forward_fp))
            self.assertTrue(os.path.getsize(qual_forward_fp) > 0)
            with open(index_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<td>Minimum:</td><td>1</td>', html)
                self.assertIn('<td>Maximum:</td><td>3</td>', html)
            with open(csv_fp, 'r') as ch:
                csv = ch.read()
                self.assertIn('sample_1', csv)

    def test_single_sample(self):
        bsi = BarcodeSequenceFastqIterator(self.barcodes[:1],
                                           self.sequences[:1])

        barcode_map = pd.Series(
            ['AAAA'], name='bc',
            index=pd.Index(['sample1'], name='id')
        )
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        demux_data, ecc = emp_single(bsi, barcode_map,
                                     golay_error_correction=False)
        # test that an index.html file is created and that it has size > 0
        with tempfile.TemporaryDirectory() as output_dir:
            # TODO: Remove _PlotQualView wrapper
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=False), n=1)
            self.assertTrue(result is None)
            index_fp = os.path.join(output_dir, 'overview.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.getsize(index_fp) > 0)
            csv_fp = os.path.join(output_dir, 'per-sample-fastq-counts.csv')
            self.assertTrue(os.path.exists(csv_fp))
            self.assertTrue(os.path.getsize(csv_fp) > 0)
            pdf_fp = os.path.join(output_dir, 'demultiplex-summary.pdf')
            self.assertFalse(os.path.exists(pdf_fp))
            png_fp = os.path.join(output_dir, 'demultiplex-summary.png')
            self.assertFalse(os.path.exists(png_fp))
            with open(index_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<td>Minimum:</td><td>1</td>', html)
                self.assertIn('<td>Maximum:</td><td>1</td>', html)

    def test_single_sample_multiple_files(self):
        # Note, this case came up on the QIIME 2 Forum, due to a user running
        # `summarize` with demuxed sequences that all had the same Sample ID
        # in the internal MANIFEST file. With versions of seaborn greater than
        # 0.8.0 this is no longer a problem, but on prior versions, the user
        # would see the following error:
        # TypeError: len() of unsized object
        files = ['sample1_S0_L001_R1_001.fastq.gz',
                 'sample2_S0_L001_R1_001.fastq.gz',
                 'sample3_S0_L001_R1_001.fastq.gz',
                 'MANIFEST', 'metadata.yml']
        for f in files:
            shutil.copy(
                self.get_data_path('single_sample_multiple_files/%s' % f),
                self.temp_dir.name)

        demux_data = SingleLanePerSampleSingleEndFastqDirFmt(
                self.temp_dir.name, mode='r')
        with tempfile.TemporaryDirectory() as output_dir:
            # TODO: Remove _PlotQualView wrapper
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=False), n=1)
            self.assertTrue(result is None)
            index_fp = os.path.join(output_dir, 'overview.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.getsize(index_fp) > 0)
            csv_fp = os.path.join(output_dir, 'per-sample-fastq-counts.csv')
            self.assertTrue(os.path.exists(csv_fp))
            self.assertTrue(os.path.getsize(csv_fp) > 0)
            pdf_fp = os.path.join(output_dir, 'demultiplex-summary.pdf')
            self.assertTrue(os.path.exists(pdf_fp))
            png_fp = os.path.join(output_dir, 'demultiplex-summary.png')
            self.assertTrue(os.path.exists(png_fp))
            with open(index_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<td>Minimum:</td><td>1</td>', html)
                self.assertIn('<td>Maximum:</td><td>1</td>', html)

    def test_paired_end(self):
        barcodes = self.barcodes[:3]

        forward = self.sequences[:3]

        reverse = [('@s1/1 abc/1', 'CCC', '+', 'YYY'),
                   ('@s2/1 abc/1', 'GGG', '+', 'PPP'),
                   ('@s3/1 abc/1', 'TTT', '+', 'PPP')]

        bpsi = BarcodePairedSequenceFastqIterator(barcodes, forward, reverse)

        barcode_map = pd.Series(
            ['AAAA', 'AACC', 'TTAA'], name='bc',
            index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')
        )
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        demux_data, ecc = emp_paired(bpsi, barcode_map,
                                     golay_error_correction=False)
        with tempfile.TemporaryDirectory() as output_dir:
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=True), n=2)
            self.assertTrue(result is None)
            plot_fp = os.path.join(output_dir, 'quality-plot.html')
            qual_forward_fp = os.path.join(
                output_dir, 'forward-seven-number-summaries.csv')
            self.assertTrue(os.path.exists(qual_forward_fp))
            self.assertTrue(os.path.getsize(qual_forward_fp) > 0)
            qual_reverse_fp = os.path.join(
                output_dir, 'reverse-seven-number-summaries.csv')
            self.assertTrue(os.path.exists(qual_reverse_fp))
            self.assertTrue(os.path.getsize(qual_reverse_fp) > 0)
            with open(plot_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<h5 class="text-center">Forward Reads</h5>',
                              html)
                self.assertIn('<h5 class="text-center">Reverse Reads</h5>',
                              html)

    def test_subsample_higher_than_seqs_count(self):
        barcodes = self.barcodes[:1]

        sequences = self.sequences[:1]
        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)

        barcode_map = pd.Series(['AAAA'], name='bc',
                                index=pd.Index(['sample1'], name='id'))
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        demux_data, ecc = emp_single(bsi, barcode_map,
                                     golay_error_correction=False)
        with tempfile.TemporaryDirectory() as output_dir:
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=False), n=50)
            self.assertTrue(result is None)
            plot_fp = os.path.join(output_dir, 'quality-plot.html')
            with open(plot_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<strong>Warning:</strong>', html)

    def test_phred_score_out_of_range(self):
        barcodes = self.barcodes[:3]

        sequences = [('@s1/1 abc/1', 'GGG', '+', 'jjj'),
                     ('@s2/1 abc/1', 'CCC', '+', 'iii'),
                     ('@s3/1 abc/1', 'AAA', '+', 'hhh')]
        bsi = BarcodeSequenceFastqIterator(barcodes, sequences)

        barcode_map = pd.Series(
            ['AAAA', 'AACC', 'TTAA'], name='bc',
            index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')
        )
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        demux_data, ecc = emp_single(bsi, barcode_map,
                                     golay_error_correction=False)
        with tempfile.TemporaryDirectory() as output_dir:
            result = summarize(output_dir, _PlotQualView(demux_data,
                                                         paired=False), n=50)
            self.assertTrue(result is None)
            plot_fp = os.path.join(output_dir, 'quality-plot.html')
            with open(plot_fp, 'r') as fh:
                html = fh.read()
                self.assertIn('<strong>Danger:</strong>', html)

    def test_inconsistent_sequence_length_single(self):
        sequences = [('@s1/1 abc/1', 'GGGGGGG', '+', 'YYYYYYY'),
                     ('@s2/1 abc/1', 'CCCCC', '+', 'PPPPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc/1', 'T', '+', 'P')]
        bsi = BarcodeSequenceFastqIterator(self.barcodes, sequences)

        barcode_map = pd.Series(
            ['AAAA', 'AACC'], name='bc',
            index=pd.Index(['sample1', 'sample2'], name='id')
        )
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        demux_data, ecc = emp_single(bsi, barcode_map,
                                     golay_error_correction=False)
        lengths = [1, 3, 5, 7]
        for n in range(1, 6):
            with tempfile.TemporaryDirectory() as output_dir:
                lengths_ = lengths[0:5-n] if n < 4 else [1]
                # TODO: Remove _PlotQualView wrapper
                summarize(output_dir, _PlotQualView(demux_data,
                                                    paired=False), n=n)
                plot_fp = os.path.join(output_dir, 'data.jsonp')
                with open(plot_fp, 'r') as fh:
                    jsonp = fh.read()
                    json_ = jsonp.replace('app.init(',
                                          '[').replace(');', ']')
                    payload = json.loads(json_)[0]
                    self.assertEqual(payload["totalSeqCount"], 4)
                    self.assertIn(payload["minSeqLen"]["forward"], lengths_)
                    self.assertEqual(payload["minSeqLen"]["reverse"], None)
                    self.assertEqual(payload["n"], min(n, 4))

    def test_inconsistent_sequence_length_paired(self):
        forward = [('@s1/1 abc/1', 'G', '+', 'Y'),
                   ('@s2/1 abc/1', 'CCC', '+', 'PPP'),
                   ('@s3/1 abc/1', 'AAAAA', '+', 'PPPPP'),
                   ('@s4/1 abc/1', 'TTTTTTT', '+', 'PPPPPPP')]
        reverse = [('@s1/1 abc/1', 'AAAAAAA', '+', 'YYYYYYY'),
                   ('@s2/1 abc/1', 'TTTTT', '+', 'PPPPP'),
                   ('@s3/1 abc/1', 'GGG', '+', 'PPP'),
                   ('@s4/1 abc/1', 'C', '+', 'P')]
        bpsi = BarcodePairedSequenceFastqIterator(self.barcodes, forward,
                                                  reverse)

        barcode_map = pd.Series(
            ['AAAA', 'AACC'], name='bc',
            index=pd.Index(['sample1', 'sample2'], name='id')
        )
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        demux_data, ecc = emp_paired(bpsi, barcode_map,
                                     golay_error_correction=False)
        lengths = [1, 3, 5, 7]
        for n in range(1, 6):
            with tempfile.TemporaryDirectory() as output_dir:
                lengths_ = lengths[0:5-n] if n < 4 else [1]
                # TODO: Remove _PlotQualView wrapper
                summarize(output_dir, _PlotQualView(demux_data,
                                                    paired=True), n=n)
                plot_fp = os.path.join(output_dir, 'data.jsonp')
                with open(plot_fp, 'r') as fh:
                    jsonp = fh.read()
                    json_ = jsonp.replace('app.init(',
                                          '[').replace(');', ']')
                    payload = json.loads(json_)[0]
                    self.assertEqual(payload["totalSeqCount"], 4)
                    self.assertIn(payload["minSeqLen"]["forward"], lengths_)
                    self.assertIn(payload["minSeqLen"]["reverse"], lengths_)
                    self.assertEqual(payload["n"], min(n, 4))

    def test_sequence_length_uses_subsample_single(self):
        random.seed(6)  # Will select s1 and s2 which aren't the shortest ones

        sequences = [('@s1/1 abc/1', 'GGGGGGG', '+', 'YYYYYYY'),
                     ('@s2/1 abc/1', 'CCCCC', '+', 'PPPPP'),
                     ('@s3/1 abc/1', 'AAA', '+', 'PPP'),
                     ('@s4/1 abc/1', 'T', '+', 'P')]
        bsi = BarcodeSequenceFastqIterator(self.barcodes, sequences)

        barcode_map = pd.Series(
            ['AAAA', 'AACC'], name='bc',
            index=pd.Index(['sample1', 'sample2'], name='id')
        )
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        demux_data, ecc = emp_single(bsi, barcode_map,
                                     golay_error_correction=False)
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, _PlotQualView(demux_data, paired=False), n=2)
            plot_fp = os.path.join(output_dir, 'data.jsonp')
            with open(plot_fp, 'r') as fh:
                jsonp = fh.read()
                json_ = jsonp.replace('app.init(',
                                      '[').replace(');', ']')
                payload = json.loads(json_)[0]
                self.assertEqual(payload["minSeqLen"]["forward"], 5)
                self.assertEqual(payload["minSeqLen"]["reverse"], None)

    def test_sequence_length_uses_subsample_paired(self):
        random.seed(6)  # Will select s1 and s2 which aren't the shortest pairs

        forward = [('@s1/1 abc/1', 'GGG', '+', 'YYY'),
                   ('@s2/1 abc/1', 'CCCCC', '+', 'PPPPP'),
                   ('@s3/1 abc/1', 'A', '+', 'P'),
                   ('@s4/1 abc/1', 'TTTTTTT', '+', 'PPPPPPP')]
        reverse = [('@s1/1 abc/1', 'AAAAA', '+', 'YYYYY'),
                   ('@s2/1 abc/1', 'TTTTTTT', '+', 'PPPPPPP'),
                   ('@s3/1 abc/1', 'GGG', '+', 'PPP'),
                   ('@s4/1 abc/1', 'C', '+', 'P')]
        bpsi = BarcodePairedSequenceFastqIterator(self.barcodes, forward,
                                                  reverse)

        barcode_map = pd.Series(
            ['AAAA', 'AACC'], name='bc',
            index=pd.Index(['sample1', 'sample2'], name='id')
        )
        barcode_map = qiime2.CategoricalMetadataColumn(barcode_map)

        demux_data, ecc = emp_paired(bpsi, barcode_map,
                                     golay_error_correction=False)
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, _PlotQualView(demux_data, paired=True), n=2)
            plot_fp = os.path.join(output_dir, 'data.jsonp')
            with open(plot_fp, 'r') as fh:
                jsonp = fh.read()
                json_ = jsonp.replace('app.init(',
                                      '[').replace(');', ']')
                payload = json.loads(json_)[0]
                self.assertEqual(payload["minSeqLen"]["forward"], 3)
                self.assertEqual(payload["minSeqLen"]["reverse"], 5)

    def test_empty_single_end(self):
        empty = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('summarize_empty/empty_single_end'), mode='r')
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, _PlotQualView(empty, paired=False), n=1)
        # Checkpoint assertion
        self.assertTrue(True)

    def test_empty_paired_end(self):
        empty = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path('summarize_empty/empty_paired_end'), mode='r')
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, _PlotQualView(empty, paired=True), n=1)
        # Checkpoint assertion
        self.assertTrue(True)

    def test_empty_paired_end_forward(self):
        empty = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path(
                'summarize_empty/empty_forward_in_paired_end'), mode='r')
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, _PlotQualView(empty, paired=True), n=1)
        # Checkpoint assertion
        self.assertTrue(True)

    def test_empty_paired_end_reverse(self):
        empty = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path(
                'summarize_empty/empty_reverse_in_paired_end'), mode='r')
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, _PlotQualView(empty, paired=True), n=1)
        # Checkpoint assertion
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
