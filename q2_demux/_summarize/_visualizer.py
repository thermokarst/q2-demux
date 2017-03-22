# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import os
import pkg_resources
import shutil
import random

import pandas as pd
import seaborn as sns
import numpy as np

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


def summarize(output_dir: str, data: _PlotQualView, n: int=10) -> None:
    paired = data.paired
    data = data._directory_format

    manifest = pd.read_csv(os.path.join(str(data), data.manifest.name),
                           header=0, comment='#')
    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(data), x))

    fwd = manifest[manifest.direction == 'forward'].filename.tolist()
    rev = manifest[manifest.direction == 'reverse'].filename.tolist()

    per_sample_fastq_counts = {}
    for file in fwd:
        count = 0
        for seq in _read_fastq_seqs(os.path.join(str(data), file)):
            count += 1
        sample_name = os.path.basename(file).split('_', 1)[0]
        per_sample_fastq_counts[sample_name] = count

    result = pd.Series(per_sample_fastq_counts)
    result.name = 'Sequence count'
    result.index.name = 'Sample name'
    result.sort_values(inplace=True, ascending=False)
    result.to_csv(os.path.join(output_dir, 'per-sample-fastq-counts.csv'),
                  header=True, index=True)

    quality_scores = collections.defaultdict(list)
    subsample_ns = sorted(random.sample(range(1, result.sum() + 1), n))
    target = subsample_ns.pop(0)
    current = 1
    if paired:
        for f, r in zip(sorted(fwd), sorted(rev)):
            for seq1, seq2 in zip(_read_fastq_seqs(f), _read_fastq_seqs(r)):
                if current == target:
                    quality_scores['forward'].append(
                        _decode_qual_to_phred(seq1[3]))
                    quality_scores['reverse'].append(
                        _decode_qual_to_phred(seq2[3]))
                    if subsample_ns:
                        target = subsample_ns.pop(0)
                    else:
                        target = -1
                        break
                current += 1
            if target == -1:
                break
    else:
        for f in fwd:
            for seq in _read_fastq_seqs(f):
                if current == target:
                    quality_scores['forward'].append(
                        _decode_qual_to_phred(seq[3]))
                    if subsample_ns:
                        target = subsample_ns.pop(0)
                    else:
                        target = -1
                        break
                current += 1
            if target == -1:
                break
    subsample_seqs = pd.DataFrame.from_dict(quality_scores)

    show_plot = len(fwd) > 1
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
        'paired': paired,
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
        subsample_seqs.to_json(fh, orient="records")
        fh.write(');')
