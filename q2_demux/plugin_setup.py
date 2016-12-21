# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

import qiime2
import qiime2.plugin
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import SequencesWithQuality

import q2_demux
from ._type import RawSequences
from ._format import EMPMultiplexedDirFmt, EMPMultiplexedSingleEndDirFmt


plugin = qiime2.plugin.Plugin(
    name='demux',
    version=q2_demux.__version__,
    website='https://github.com/qiime2/q2-demux',
    package='q2_demux',
    user_support_text=None,
    citation_text=None
)

plugin.register_semantic_types(RawSequences)

plugin.register_formats(EMPMultiplexedDirFmt, EMPMultiplexedSingleEndDirFmt)

plugin.register_semantic_type_to_format(
    RawSequences,
    artifact_format=EMPMultiplexedDirFmt
)

plugin.methods.register_function(
    function=q2_demux.emp,
    inputs={'seqs': RawSequences},
    parameters={'barcodes': qiime2.plugin.MetadataCategory,
                'rev_comp_barcodes': qiime2.plugin.Bool,
                'rev_comp_mapping_barcodes': qiime2.plugin.Bool},
    outputs=[('per_sample_sequences', SampleData[SequencesWithQuality])],
    name='Demultiplex sequence data generated with the EMP protocol.',
    description=('Demultiplex sequence data (i.e., map barcode reads to '
                 'sample ids) for data generated with the Earth Microbiome '
                 'Project (EMP) amplicon sequencing protocol. Details about '
                 'this protocol can be found at '
                 'http://www.earthmicrobiome.org/emp-standard-protocols/')
)

plugin.visualizers.register_function(
    function=q2_demux.summarize,
    inputs={'data': SampleData[SequencesWithQuality]},
    parameters={},
    name='Summarize counts per sample.',
    description=('Generate a summary of counts per sample from sequence data '
                 'that has already been demultiplexed')
)

importlib.import_module('q2_demux._transformer')
