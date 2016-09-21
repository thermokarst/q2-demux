import tempfile
import shutil
import os.path
import itertools
import gzip

import skbio
from q2_types.per_sample_sequences import FastqGzFormat

from .plugin_setup import plugin
from ._demux import BarcodeSequenceIterator
from ._format import EMPMultiplexedDirFmt, EMPMultiplexedSingleEndDirFmt

def _read_fastq_seqs(filepath):
    # This function is adapted from @jairideout's SO post:
    # http://stackoverflow.com/a/39302117/3424666
    fh = gzip.open(filepath, 'rt')
    for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
        yield (seq_header.strip(), seq.strip(), qual_header.strip(),
               qual.strip())

@plugin.register_transformer
def _1(dirfmt: EMPMultiplexedDirFmt) -> BarcodeSequenceIterator:
    barcode_generator =  _read_fastq_seqs(
        str(dirfmt.barcodes.view(FastqGzFormat).path))
    sequence_generator = _read_fastq_seqs(
        str(dirfmt.sequences.view(FastqGzFormat).path))
    result = BarcodeSequenceIterator(barcode_generator, sequence_generator)
    # ensure that dirfmt stays in scope as long as result does so these
    # generators will work.
    result.__dirfmt = dirfmt
    return result


@plugin.register_transformer
def _2(dirfmt: EMPMultiplexedSingleEndDirFmt) -> EMPMultiplexedDirFmt:
    result = tempfile.mkdtemp()
    sequences_fp = os.path.join(result,
                                'sequences.fastq.gz')
    barcodes_fp = os.path.join(result,
                               'barcodes.fastq.gz')
    shutil.copyfile(str(dirfmt.sequences.view(FastqGzFormat).path),
                    sequences_fp)
    shutil.copyfile(str(dirfmt.barcodes.view(FastqGzFormat).path),
                    barcodes_fp)
    return result
