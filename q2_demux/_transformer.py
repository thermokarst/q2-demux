import shutil

import skbio
from q2_types.per_sample_sequences import FastqGzFormat

from .plugin_setup import plugin
from ._demux import BarcodeSequenceIterator
from ._format import EMPMultiplexedDirFmt, EMPMultiplexedSingleEndDirFmt


@plugin.register_transformer
def _1(dirfmt: EMPMultiplexedDirFmt) -> BarcodeSequenceIterator:
    barcode_generator = skbio.io.read(
        str(dirfmt.barcodes.view(FastqGzFormat)), format='fastq',
        constructor=skbio.DNA, phred_offset=33, verify=False)
    sequence_generator = skbio.io.read(
        str(dirfmt.sequences.view(FastqGzFormat)), format='fastq',
        constructor=skbio.DNA, phred_offset=33, verify=False)
    result = BarcodeSequenceIterator(barcode_generator, sequence_generator)
    # ensure that dirfmt stays in scope as long as result does so these
    # generators will work.
    result.__dirfmt = dirfmt
    return result


@plugin.register_transformer
def _2(dirfmt: EMPMultiplexedSingleEndDirFmt) -> EMPMultiplexedDirFmt:
    # TODO: revisit this API to simpify defining transformers
    result = EMPMultiplexedDirFmt().path

    sequences_fp = str(result / 'sequences.fastq.gz')
    barcodes_fp = str(result / 'barcodes.fastq.gz')
    shutil.copyfile(str(dirfmt.sequences.view(FastqGzFormat)), sequences_fp)
    shutil.copyfile(str(dirfmt.barcodes.view(FastqGzFormat)), barcodes_fp)

    return result
