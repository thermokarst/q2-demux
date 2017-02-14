# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

import pandas as pd
import seaborn as sns

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt, FastqGzFormat
)
from q2_demux._demux import _read_fastq_seqs
import q2templates


TEMPLATES = pkg_resources.resource_filename('q2_demux', '_summarize')


def summarize(output_dir: str, data: SingleLanePerSampleSingleEndFastqDirFmt) \
        -> None:
    per_sample_fastqs = list(data.sequences.iter_views(FastqGzFormat))
    per_sample_fastq_counts = {}
    for relpath, view in per_sample_fastqs:
        seqs = _read_fastq_seqs(str(view))
        count = 0
        for seq in seqs:
            count += 1
        sample_name = relpath.name.split('_', 1)[0]
        per_sample_fastq_counts[sample_name] = count

    result = pd.Series(per_sample_fastq_counts)
    result.name = 'Sequence count'
    result.index.name = 'Sample name'
    result.sort_values(inplace=True, ascending=False)
    result.to_csv(os.path.join(output_dir, 'per-sample-fastq-counts.csv'),
                  header=True, index=True)

    show_plot = len(per_sample_fastqs) > 1
    if show_plot:
        ax = sns.distplot(result, kde=False)
        ax.set_xlabel('Number of sequences')
        ax.set_ylabel('Frequency')
        fig = ax.get_figure()
        fig.savefig(os.path.join(output_dir, 'demultiplex-summary.png'))
        fig.savefig(os.path.join(output_dir, 'demultiplex-summary.pdf'))

    html = result.to_frame().to_html(classes='table table-striped table-hover')
    html = html.replace('border="1"', 'border="0"')
    index = os.path.join(TEMPLATES, 'assets', 'index.html')
    overview_template = os.path.join(TEMPLATES, 'assets', 'overview.html')
    quality_template = os.path.join(TEMPLATES, 'assets', 'quality-plot.html')
    context = {
        'result_data': {
            'min': result.min(),
            'median': result.median(),
            'mean': result.mean(),
            'max': result.max(),
            'sum': result.sum()
        },
        'result': html,
        'show_plot': show_plot,
        'tabs': [{'title': 'Overview',
                  'url': 'overview.html'},
                 {'title': 'Interactive Quality Plot',
                  'url': 'quality-plot.html'}]
    }
    templates = [index, overview_template, quality_template]
    q2templates.render(templates, output_dir, context=context)
