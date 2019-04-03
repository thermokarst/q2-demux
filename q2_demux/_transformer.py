# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import shutil

import pandas as pd
import numpy as np
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqGzFormat)
from qiime2 import Metadata

from .plugin_setup import plugin
from ._demux import (BarcodeSequenceFastqIterator,
                     BarcodePairedSequenceFastqIterator, _read_fastq_seqs)
from ._format import (EMPMultiplexedDirFmt,
                      EMPSingleEndDirFmt, EMPSingleEndCasavaDirFmt,
                      EMPPairedEndDirFmt, EMPPairedEndCasavaDirFmt,
                      ErrorCorrectionDetailsFmt)
from ._summarize import _PlotQualView


@plugin.register_transformer
def _1(dirfmt: EMPSingleEndDirFmt) -> BarcodeSequenceFastqIterator:
    barcode_generator = _read_fastq_seqs(
        str(dirfmt.barcodes.view(FastqGzFormat)))
    sequence_generator = _read_fastq_seqs(
        str(dirfmt.sequences.view(FastqGzFormat)))
    result = BarcodeSequenceFastqIterator(barcode_generator,
                                          sequence_generator)
    # ensure that dirfmt stays in scope as long as result does so these
    # generators will work.
    result.__dirfmt = dirfmt
    return result


# TODO: remove this when names are aliased
@plugin.register_transformer
def _1_legacy(dirfmt: EMPMultiplexedDirFmt) -> BarcodeSequenceFastqIterator:
    return _1(dirfmt)


# NOTE: a legacy transformer isn't needed for EMPMultiplexedSingleEndDirFmt
# as no artifacts exist in this form, it is used for import only.
@plugin.register_transformer
def _2(dirfmt: EMPSingleEndCasavaDirFmt) -> EMPSingleEndDirFmt:
    # TODO: revisit this API to simpify defining transformers
    result = EMPMultiplexedDirFmt().path

    sequences_fp = str(result / 'sequences.fastq.gz')
    barcodes_fp = str(result / 'barcodes.fastq.gz')
    shutil.copyfile(str(dirfmt.sequences.view(FastqGzFormat)), sequences_fp)
    shutil.copyfile(str(dirfmt.barcodes.view(FastqGzFormat)), barcodes_fp)

    return result


@plugin.register_transformer
def _3(dirfmt: EMPPairedEndCasavaDirFmt) -> EMPPairedEndDirFmt:
    result = EMPMultiplexedDirFmt()
    root = result.path

    forward_fp = str(root / 'forward.fastq.gz')
    reverse_fp = str(root / 'reverse.fastq.gz')
    barcodes_fp = str(root / 'barcodes.fastq.gz')
    shutil.copyfile(str(dirfmt.forward.view(FastqGzFormat)), forward_fp)
    shutil.copyfile(str(dirfmt.reverse.view(FastqGzFormat)), reverse_fp)
    shutil.copyfile(str(dirfmt.barcodes.view(FastqGzFormat)), barcodes_fp)

    return result


@plugin.register_transformer
def _4(dirfmt: EMPPairedEndDirFmt) -> BarcodePairedSequenceFastqIterator:
    barcode_generator = _read_fastq_seqs(
        str(dirfmt.barcodes.view(FastqGzFormat)))
    forward_generator = _read_fastq_seqs(
        str(dirfmt.forward.view(FastqGzFormat)))
    reverse_generator = _read_fastq_seqs(
        str(dirfmt.reverse.view(FastqGzFormat)))
    result = BarcodePairedSequenceFastqIterator(barcode_generator,
                                                forward_generator,
                                                reverse_generator)
    # ensure that dirfmt stays in scope as long as result does so these
    # generators will work.
    result.__dirfmt = dirfmt
    return result


# TODO: Remove _PlotQualView once QIIME 2 #220 completed
@plugin.register_transformer
def _5(dirfmt: SingleLanePerSampleSingleEndFastqDirFmt) -> _PlotQualView:
    return _PlotQualView(dirfmt, paired=False)


@plugin.register_transformer
def _6(dirfmt: SingleLanePerSamplePairedEndFastqDirFmt) -> _PlotQualView:
    return _PlotQualView(dirfmt, paired=True)


@plugin.register_transformer
def _7(dirfmt: EMPPairedEndDirFmt) -> BarcodeSequenceFastqIterator:
    barcode_generator = _read_fastq_seqs(
        str(dirfmt.barcodes.view(FastqGzFormat)))
    sequence_generator = _read_fastq_seqs(
        str(dirfmt.forward.view(FastqGzFormat)))
    result = BarcodeSequenceFastqIterator(barcode_generator,
                                          sequence_generator)
    # ensure that dirfmt stays in scope as long as result does so these
    # generators will work.
    result.__dirfmt = dirfmt
    return result


_ec_details_column_dtypes = {
    'sample-id': np.str,
    'barcode-sequence-id': np.str,
    'barcode-uncorrected': np.str,
    'barcode-corrected': np.str,
    'barcode-errors': np.number
}


def _ec_details_to_df(ff):
    # https://github.com/pandas-dev/pandas/issues/9435
    df = pd.read_csv(str(ff), dtype=_ec_details_column_dtypes)
    df.set_index('sample-id', inplace=True)
    return df


@plugin.register_transformer
def _8(data: pd.DataFrame) -> ErrorCorrectionDetailsFmt:
    ff = ErrorCorrectionDetailsFmt()
    data.to_csv(str(ff))
    return ff


@plugin.register_transformer
def _9(ff: ErrorCorrectionDetailsFmt) -> pd.DataFrame:
    return _ec_details_to_df(ff)


@plugin.register_transformer
def _10(ff: ErrorCorrectionDetailsFmt) -> Metadata:
    return Metadata(_ec_details_to_df(ff))
