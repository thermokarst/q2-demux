import tempfile
import shutil
import os.path

import skbio
from q2_types.per_sample_sequences import FastqGzFormat

from .plugin_setup import plugin
from ._demux import BarcodeSequenceIterator
from ._format import EMPMultiplexedDirFmt, EMPMultiplexedSingleEndDirFmt


@plugin.register_transformer
def _1(dirfmt: EMPMultiplexedDirFmt) -> BarcodeSequenceIterator:
    barcode_generator = skbio.io.read(
        str(dirfmt.barcodes.view(FastqGzFormat).path), format='fastq',
        constructor=skbio.DNA, phred_offset=33)
    sequence_generator = skbio.io.read(
        str(dirfmt.sequences.view(FastqGzFormat).path), format='fastq',
        constructor=skbio.DNA, phred_offset=33)
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
