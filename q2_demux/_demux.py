import os.path
import gzip
import yaml
import itertools
import collections

import skbio
import pandas as pd

import qiime
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt, FastqManifestFormat,
    YamlFormat, FastqGzFormat)


class BarcodeSequenceIterator(collections.Iterator):
    def __init__(self, barcode_generator, sequence_generator):
        self.barcode_generator = barcode_generator
        self.sequence_generator = sequence_generator

    def __iter__(self):
        # Adapted from q2-types
        for barcode, sequence in itertools.zip_longest(
                self.barcode_generator, self.sequence_generator):
            if barcode is None:
                raise ValueError('More sequences than barcodes.')
            if sequence is None:
                raise ValueError('More barcodes than sequences.')
            if barcode.metadata['id'] != sequence.metadata['id']:
                raise ValueError(
                    'Mismatched sequence ids: %s and %s' %
                    (barcode.metadata['id'], sequence.metadata['id']))
            if barcode.metadata['description'] != \
               sequence.metadata['description']:
                raise ValueError(
                    'Mismatched sequence descriptions: %s and %s' %
                    (barcode.metadata['description'],
                     sequence.metadata['description']))
            yield barcode, sequence

    def __next__(self):
        return next(self.barcode_generator), next(self.sequence_generator)


def summary(output_dir: str, data: SingleLanePerSampleSingleEndFastqDirFmt) \
        -> None:
    per_sample_fastqs = list(data.sequences.iter_views(FastqGzFormat))
    per_sample_fastq_counts = {}
    for per_sample_fastq in per_sample_fastqs:
        seqs = skbio.io.read(per_sample_fastq[1].open(),
                             format='fastq',
                             phred_offset=33, compression='gzip',
                             constructor=skbio.DNA)
        count = 0
        for seq in seqs:
            count += 1
        per_sample_fastq_counts[per_sample_fastq[0]] = count
    result = pd.Series(per_sample_fastq_counts)
    result.sort_values(inplace=True, ascending=False)
    result.to_csv(os.path.join(output_dir,
                  'per-sample-fastq-counts.csv'))
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        fh.write(' <h1>Demultiplexed sequence counts summary</h1>\n')
        fh.write(' Minimum: %d<br>\n' % result.min())
        fh.write(' Median: %d<br>\n' % result.max())
        fh.write(' Mean: %1.3f<br>\n' % result.mean())
        fh.write(' Maximum: %1.3f<p>\n\n' % result.median())
        fh.write(' <h1>Demultiplexed sequence counts detail</h1>\n')
        fh.write(result.to_frame('Counts').to_html())
        fh.write(' <a href="./per-sample-fastq-counts.csv">csv</a>\n')
        fh.write('</body></html>')


def emp(seqs: BarcodeSequenceIterator,
        barcodes: qiime.MetadataCategory,
        rev_comp_barcodes: bool=False,
        rev_comp_mapping_barcodes: bool=False) \
        -> SingleLanePerSampleSingleEndFastqDirFmt:
    barcode_map = {}
    barcode_len = None
    for sample_id, barcode in barcodes.to_series().iteritems():
        if barcode_len is None:
            barcode_len = len(barcode)
        elif len(barcode) != barcode_len:
            raise ValueError('Barcodes of different lengths were detected: '
                             '%d != %d. Variable length barcodes are not '
                             'supported.' % (len(barcode), barcode_len))
        if rev_comp_mapping_barcodes:
            barcode = str(skbio.DNA(barcode).reverse_complement())
        if barcode in barcode_map:
            raise ValueError('A duplicate barcode was detected. The barcode '
                             '%s was observed for samples %s and %s.'
                             % (barcode, sample_id, barcode_map[barcode]))
        barcode_map[barcode] = sample_id

    result = SingleLanePerSampleSingleEndFastqDirFmt()

    manifest = FastqManifestFormat()
    manifest_fh = manifest.open()
    manifest_fh.write('sample-id,filename,direction\n')

    per_sample_fastqs = {}

    for barcode, sequence in seqs:
        if rev_comp_barcodes:
            barcode = barcode.reverse_complement()
        barcode = str(barcode)
        barcode = barcode[:barcode_len]

        try:
            sample_id = barcode_map[barcode]
        except KeyError:
            # TODO: this should ultimately be logged, but we don't currently
            # have support for that.
            continue

        if sample_id not in per_sample_fastqs:
            # The barcode id, lane number and read number are not relevant
            # here. We might ultimately want to use a dir format other than
            # SingleLanePerSampleSingleEndFastqDirFmt which doesn't care
            # about this information. Similarly, the direction of the read
            # isn't relevant here anymore.
            barcode_id = len(per_sample_fastqs) + 1
            path = result.sequences.path_maker(sample_id=sample_id,
                                               barcode_id=barcode_id,
                                               lane_number=1,
                                               read_number=1)
            per_sample_fastqs[sample_id] = gzip.open(str(path), mode='w')
            manifest_fh.write('%s,%s,%s\n' % (sample_id, path.name, 'forward'))

        # skbio can't currently write to a gzip file handle, so this is a
        # hack to do that.
        fastq_lines = ''.join(
             sequence.write([], format='fastq', phred_offset=33))
        fastq_lines = fastq_lines.encode('utf-8')
        per_sample_fastqs[sample_id].write(fastq_lines)

    if len(per_sample_fastqs) == 0:
        raise ValueError('No valid barcodes were identified.')

    for fh in per_sample_fastqs.values():
        fh.close()

    manifest_fh.close()
    result.manifest.write_data(manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    result.metadata.write_data(metadata, YamlFormat)

    return result
