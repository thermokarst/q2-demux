# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
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
import json

import pandas as pd
import seaborn as sns
import numpy as np

from q2_demux._demux import _read_fastq_seqs
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_demux', '_summarize')


def _decode_qual_to_phred33(qual_str):
    # this function is adapted from scikit-bio
    qual = np.frombuffer(qual_str.encode('ascii'), dtype=np.uint8) - 33
    return qual


# TODO: Remove _PlotQualView once QIIME 2 #220 completed
class _PlotQualView:
    """
    A very simple pass-through view which is made up of a single-end or
    paired-end directory format with a bool indicating if single or paired.
    """
    def __init__(self, directory_format, paired):
        self.directory_format = directory_format
        self.paired = paired


def _link_sample_n_to_file(file_records, counts, subsample_ns, direction):
    results = collections.defaultdict(list)
    for num in subsample_ns:
        total = 0
        for record in file_records[direction]:
            sample_id = record['sample_id']
            filename = record['filename']

            total += counts[direction][sample_id]
            if num < total:
                idx = counts[direction][sample_id] - (total - num)
                results[filename].append(idx)
                break
    return results


def _subsample_paired(fastq_map):
    qual_sample = collections.defaultdict(list)
    min_seq_len = {'forward': float('inf'), 'reverse': float('inf')}
    for fwd, rev, index in fastq_map:
        file_pair = zip(_read_fastq_seqs(fwd), _read_fastq_seqs(rev))
        for i, (fseq, rseq) in enumerate(file_pair):
            if i == index[0]:
                min_seq_len['forward'] = min(min_seq_len['forward'],
                                             len(fseq[1]))
                min_seq_len['reverse'] = min(min_seq_len['reverse'],
                                             len(rseq[1]))
                qual_sample['forward'].append(_decode_qual_to_phred33(fseq[3]))
                qual_sample['reverse'].append(_decode_qual_to_phred33(rseq[3]))
                index.pop(0)
                if len(index) == 0:
                    break
    return qual_sample, min_seq_len


def _subsample_single(fastq_map):
    qual_sample = collections.defaultdict(list)
    min_seq_len = {'forward': float('inf'), 'reverse': None}
    for file, index in fastq_map:
        for i, seq in enumerate(_read_fastq_seqs(file)):
            if i == index[0]:
                min_seq_len['forward'] = min(min_seq_len['forward'],
                                             len(seq[1]))
                qual_sample['forward'].append(_decode_qual_to_phred33(seq[3]))
                index.pop(0)
                if len(index) == 0:
                    break
    return qual_sample, min_seq_len


def _compute_stats_of_df(df):
    df_stats = df.describe(
        percentiles=[0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98])
    drop_cols = df_stats.index.isin(['std', 'mean', 'min', 'max'])
    df_stats = df_stats[~drop_cols]
    return df_stats


def _build_seq_len_table(qscores: pd.DataFrame) -> str:
    sequence_lengths = qscores.notnull().sum(axis=1).copy()
    stats = _compute_stats_of_df(sequence_lengths)

    stats[stats.index != 'count'] = \
        stats[stats.index != 'count'].astype(int).apply('{} nts'.format)

    stats.rename(index={'50%': '50% (Median)',
                        'count': 'Total Sequences Sampled'},
                 inplace=True)
    frame = stats.to_frame(name="")
    return q2templates.df_to_html(frame)


def summarize(output_dir: str, data: _PlotQualView, n: int = 10000) -> None:
    paired = data.paired
    data = data.directory_format
    summary_columns = ['Minimum', 'Median', 'Mean', 'Maximum', 'Total']
    context = {
        'result_data': pd.DataFrame([], columns=summary_columns),
        'result': {'forward': None, 'reverse': None},
        'n_samples': {'forward': None, 'reverse': None},
        'show_plot': {'forward': None, 'reverse': None},
        'paired': paired,
        'tabs': [{'title': 'Overview', 'url': 'overview.html'}],
        'dangers': [],
        'warnings': [],
    }

    manifest = data.manifest.view(pd.DataFrame)

    directions = ['forward', 'reverse'] if paired else ['forward']
    file_records = {'forward': [], 'reverse': []}
    per_sample_fastq_counts = {'forward': {}, 'reverse': {}}
    subsample_size = {'forward': n, 'reverse': n}
    links = {'forward': {}, 'reverse': {}}

    for sample_id, row in manifest.iterrows():
        for direction in directions:
            count = 0
            filename = row[direction]
            for seq in _read_fastq_seqs(filename):
                count += 1
            per_sample_fastq_counts[direction][sample_id] = count
            file_records[direction].append({
                'filename': filename,
                'sample_id': sample_id,
            })

    for direction in directions:
        # Prepare summary
        result = pd.Series(per_sample_fastq_counts[direction])
        result.name = 'Sequence count'
        result.index.name = 'Sample name'
        result.sort_values(inplace=True, ascending=False)

        # Create a TSV
        result_fn = 'per-sample-fastq-counts-%s.tsv' % (direction, )
        result_path = os.path.join(output_dir, result_fn)
        result.to_csv(result_path, header=True, index=True, sep='\t')

        sequence_count = result.sum()
        if subsample_size[direction] > sequence_count:
            subsample_size[direction] = sequence_count
            context.warnings.append('A subsample value was provided that is '
                                    'greater than the amount of sequences '
                                    'across all samples. The plot was '
                                    'generated using all available sequences.')

        subsample_ns = sorted(random.sample(range(sequence_count), subsample_size[direction]))
        links[direction] = _link_sample_n_to_file(file_records,
                                                  per_sample_fastq_counts,
                                                  subsample_ns,
                                                  direction)

        sample_map = [(k, v) for k, v in links[direction].items()]
        quality_scores, min_seq_len = _subsample_single(sample_map)

        show_plot = len(sample_map) > 0

        ax = sns.distplot(result, kde=False)
        ax.set_xlabel('Number of sequences')
        ax.set_ylabel('Frequency')
        fig = ax.get_figure()
        fig.savefig(os.path.join(output_dir,
                                 'demultiplex-summary-%s.png' % (direction, )))
        fig.savefig(os.path.join(output_dir,
                                 'demultiplex-summary-%s.pdf' % (direction, )))

        html_df = result.to_frame().reset_index(drop=False)
        html = q2templates.df_to_html(html_df, index=False)

    df = pd.DataFrame([[result.min(), result.median(), result.mean(),
                        result.max(), sequence_count]],
                      index=[direction],
                      columns=summary_columns)
    # RESUME: looks like i am not appending in place - fix that when I resume this work.
    # Also, don't forget to transpose the table before rendering as HTML
    context['result_data'].append(df)
    context['result'][direction] = html
    context['n_samples'][direction] = result.count()
    context['show_plot'][direction] = show_plot

    index = os.path.join(TEMPLATES, 'assets', 'index.html')
    overview_template = os.path.join(TEMPLATES, 'assets', 'overview.html')
    templates = [index, overview_template]

    forward_scores = pd.DataFrame(quality_scores['forward'])
    if not forward_scores.empty:
        forward_stats = _compute_stats_of_df(forward_scores)
        forward_stats.to_csv(os.path.join(output_dir,
                             'forward-seven-number-summaries.csv'),
                             header=True, index=True)
        forward_length_table = _build_seq_len_table(forward_scores)

        if (forward_stats.loc['50%'] > 45).any():
            context['dangers'].append('Some of the PHRED quality values are '
                                      'out of range. This is likely because '
                                      'an incorrect PHRED offset was chosen '
                                      'on import of your raw data. You can '
                                      'learn how to choose your PHRED offset '
                                      'during import in the importing '
                                      'tutorial.')

        templates.append(os.path.join(TEMPLATES, 'assets',
                                      'quality-plot.html'))

        context['forward_length_table'] = forward_length_table
        context['tabs'].append({'title': 'Interactive Quality Plot',
                               'url': 'quality-plot.html'})

    # Required initilization for conditional display of the table
    reverse_length_table = None
    if paired:
        reverse_scores = pd.DataFrame(quality_scores['reverse'])
        if not reverse_scores.empty:
            reverse_stats = _compute_stats_of_df(reverse_scores)
            reverse_stats.to_csv(os.path.join(output_dir,
                                 'reverse-seven-number-summaries.csv'),
                                 header=True, index=True)
            reverse_length_table = _build_seq_len_table(reverse_scores)
            context['reverse_length_table'] = reverse_length_table
            # If we have reverse reads and not forward reads (due to some sort
            # of error) the quality template will still need to be added
            quality_template = os.path.join(TEMPLATES, 'assets',
                                            'quality-plot.html')
            if quality_template not in templates:
                templates.append(quality_template)

    print(context['result_data'])
    context['result_data'] = q2templates.df_to_html(context['result_data'])
    q2templates.render(templates, output_dir, context=context)

    shutil.copytree(os.path.join(TEMPLATES, 'assets', 'dist'),
                    os.path.join(output_dir, 'dist'))

    with open(os.path.join(output_dir, 'data.jsonp'), 'w') as fh:
        fh.write("app.init(")
        json.dump({'n': int(n), 'totalSeqCount': int(sequence_count),
                   'minSeqLen': min_seq_len}, fh)
        fh.write(',')
        if not forward_scores.empty:
            forward_stats.to_json(fh)
        else:
            forward_scores.to_json(fh)
        if paired and not reverse_scores.empty:
            fh.write(',')
            reverse_stats.to_json(fh)
        fh.write(');')
