# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import shutil
import random

import pandas as pd
import seaborn as sns
import numpy as np

from q2_types.per_sample_sequences import FastqGzFormat
from q2_demux._demux import _read_fastq_seqs
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_demux', '_summarize')


def _decode_qual_to_phred(qual_str):
    # this function is adapted from scikit-bio
    qual = np.fromstring(qual_str, dtype=np.uint8) - 33
    return qual


# TODO: Remove _PlotQualView once QIIME 2 #220 completed
class _PlotQualView:
    """
    A very simple pass-through view which is made up of a single-end or
    paired-end directory format with a bool indicating if single or paired.
    """
    def __init__(self, directory_format, paired):
        self._directory_format = directory_format
        self.directory = str(directory_format)
        self.paired = paired


def summarize(output_dir: str, data: _PlotQualView) -> None:
    paired = data.paired
    data = data._directory_format

    per_sample_fastqs = list(data.sequences.iter_views(FastqGzFormat))
    per_sample_fastq_counts = {}
    sample_map = {}
    for relpath, view in per_sample_fastqs:
        count = 0
        for seq in _read_fastq_seqs(str(view)):
            count += 1
        sample_name = relpath.name.split('_', 1)[0]
        per_sample_fastq_counts[sample_name] = count

    result = pd.Series(per_sample_fastq_counts)
    result.name = 'Sequence count'
    result.index.name = 'Sample name'
    result.sort_values(inplace=True, ascending=False)
    result.to_csv(os.path.join(output_dir, 'per-sample-fastq-counts.csv'),
                  header=True, index=True)

    quality_scores = {}
    subsample_ns = sorted(random.sample(range(1, result.sum() + 1), 10))
    target = subsample_ns.pop(0)
    current = 1
    for relpath, view in per_sample_fastqs:
        for seq in _read_fastq_seqs(str(view)):
            if current == target:
                quality_scores[seq[0]] = _decode_qual_to_phred(seq[3])
                if subsample_ns:
                    target = subsample_ns.pop(0)
                else:
                    break
            current += 1
    subsample_seqs = pd.DataFrame.from_dict(quality_scores, orient='index')

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

    shutil.copytree(os.path.join(TEMPLATES, 'assets', 'app'),
                    os.path.join(output_dir, 'app'))

    with open(os.path.join(output_dir, 'data.jsonp'), 'w') as fh:
        fh.write("app.init(")
        subsample_seqs.to_json(fh, orient='values')
        fh.write(');')
