import os.path
import gzip
import yaml
import itertools
import collections

import skbio
import pandas as pd
import seaborn as sns

import qiime
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt, FastqManifestFormat, YamlFormat,
    FastqGzFormat)


def _read_fastq_seqs(filepath):
    # This function is adapted from @jairideout's SO post:
    # http://stackoverflow.com/a/39302117/3424666
    fh = gzip.open(filepath, 'rt')
    for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
        yield (seq_header.strip(), seq.strip(), qual_header.strip(),
               qual.strip())


def _extract_common_header(header):
    return header.rsplit('/', 1)[0]


class BarcodeSequenceFastqIterator(collections.Iterator):
    def __init__(self, barcode_generator, sequence_generator):
        self.barcode_generator = barcode_generator
        self.sequence_generator = sequence_generator

    def __iter__(self):
        # Adapted from q2-types
        for barcode_fields, sequence_fields in itertools.zip_longest(
                self.barcode_generator, self.sequence_generator):
            if barcode_fields is None:
                raise ValueError('More sequences were provided than barcodes.')
            if sequence_fields is None:
                raise ValueError('More barcodes were provided than sequences.')
            # The id or description fields may end with "/read-number", which
            # and read number will differ between the sequence and barcode
            # reads. Confirm that they are identical up until the last /
            barcode_header_fields = \
                barcode_fields[0][1:].split(' ', maxsplit=1)
            sequence_header_fields = \
                sequence_fields[0][1:].split(' ', maxsplit=1)

            # confirm that both barcode and sequene have header lines with
            # an equal number of space-separated fields
            if len(barcode_header_fields) != len(sequence_header_fields):
                raise ValueError('Mismatched header lines: %s and %s' %
                                 (barcode_fields[0], sequence_fields[0]))

            # confirm that the id fields are equal
            if _extract_common_header(barcode_header_fields[0]) != \
               _extract_common_header(sequence_header_fields[0]):
                raise ValueError(
                    'Mismatched sequence ids: %s and %s' %
                    (_extract_common_header(barcode_header_fields[0]),
                     _extract_common_header(sequence_header_fields[0])))

            # if a description field is present, confirm that they're equal
            # note that the number of fields can only be 1 or 2, due to
            # maxsplit
            if len(barcode_header_fields) == 2:
                if _extract_common_header(barcode_header_fields[1]) != \
                   _extract_common_header(sequence_header_fields[1]):
                    raise ValueError(
                        'Mismatched sequence descriptions: %s and %s' %
                        (_extract_common_header(barcode_header_fields[1]),
                         _extract_common_header(sequence_header_fields[1])))

            yield barcode_fields, sequence_fields

    def __next__(self):
        return next(self.barcode_generator), next(self.sequence_generator)


def summary(output_dir: str, data: SingleLanePerSampleSingleEndFastqDirFmt) \
        -> None:
    per_sample_fastqs = list(data.sequences.iter_views(FastqGzFormat))
    per_sample_fastq_counts = {}
    for per_sample_fastq in per_sample_fastqs:
        seqs = _read_fastq_seqs(str(per_sample_fastq[1]))
        count = 0
        for seq in seqs:
            count += 1
        per_sample_fastq_counts[per_sample_fastq[0]] = count

    result = pd.Series(per_sample_fastq_counts)
    result.sort_values(inplace=True, ascending=False)
    result.to_csv(os.path.join(output_dir,
                  'per-sample-fastq-counts.csv'))
    ax = sns.distplot(result, kde=False)
    ax.set_xlabel('Number of sequences')
    ax.set_xlabel('Frequency')
    fig = ax.get_figure()
    fig.savefig(os.path.join(output_dir, 'demultiplex-summary.png'))
    fig.savefig(os.path.join(output_dir, 'demultiplex-summary.pdf'))
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        fh.write(' <h1>Demultiplexed sequence counts summary</h1>\n')
        fh.write(' Minimum: %d<br>\n' % result.min())
        fh.write(' Median: %d<br>\n' % result.max())
        fh.write(' Mean: %1.3f<br>\n' % result.mean())
        fh.write(' Maximum: %1.3f<br>\n' % result.median())
        fh.write(' Total number of sequences demultiplexed: %d<p>\n\n' %
                 result.sum())
        fh.write(' <img src="demultiplex-summary.png"><br>\n')
        fh.write(' <a href="demultiplex-summary.pdf">PDF</a><p>\n\n')
        fh.write(' <h1>Demultiplexed sequence counts detail</h1>\n')
        fh.write(result.to_frame('Counts').to_html())
        fh.write(' <a href="./per-sample-fastq-counts.csv">csv</a>\n')
        fh.write('</body></html>')


def emp(seqs: BarcodeSequenceFastqIterator,
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
    manifest_fh.write('# direction is not meaningful in this file as these\n')
    manifest_fh.write('# data may be derived from forward, reverse, or \n')
    manifest_fh.write('# joined reads\n')

    per_sample_fastqs = {}

    for barcode_fields, sequence_fields in seqs:
        barcode_read = barcode_fields[1]
        if rev_comp_barcodes:
            barcode_read = str(skbio.DNA(barcode_read).reverse_complement())
        barcode_read = barcode_read[:barcode_len]

        try:
            sample_id = barcode_map[barcode_read]
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

        fastq_lines = '\n'.join(list(sequence_fields)) + '\n'
        fastq_lines = fastq_lines.encode('utf-8')
        per_sample_fastqs[sample_id].write(fastq_lines)

    if len(per_sample_fastqs) == 0:
        raise ValueError('No sequences were mapped to samples. Check that '
                         'your barcodes are in the correct orientation (see '
                         'rev_comp_barcodes and/or rev_comp_mapping_barcodes '
                         'options).')

    for sid, fh in per_sample_fastqs.items():
        fh.close()

    manifest_fh.close()
    result.manifest.write_data(manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    result.metadata.write_data(metadata, YamlFormat)

    return result
