#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Mon Apr 25, 2022 at 07:24 PM -0400
#
# Description: plot fit variables w/ w/o decay-in-flight smearing

import mplhep
import uproot
import numpy as np

from argparse import ArgumentParser
from pyTuplingUtils.utils import gen_histo_stacked_baseline, gen_histo
from pyTuplingUtils.plot import (
    plot_prepare, plot_histo, ax_add_args_histo, plot_step, ax_add_args_step
)
from pyTuplingUtils.io import read_branches_dict


################
# Configurable #
################

DEFAULT_COLORS = ['#00429d', '#5585b7', '#8bd189', '#d15d5f', '#93003a']
DEFAULT_OVERALL_COLORS = ['black', 'red']
LEGEND_LOC = {
    'q2': 'upper left',
    'mm2': 'upper left',
    'el': 'upper right'
}

MISID_WTS = ['wmis_norm']
MISID_TAGS = {
    'is_misid_pi': r'$\pi$ tag',
    'is_misid_k': r'$K$ tag',
    'is_misid_p': r'$p$ tag',
    'is_misid_e': r'$e$ tag',
    'is_misid_g': r'ghost tag',
}
PLOT_VARS = {
    'q2': r'$q^2$ [GeV$^2$]',
    'mm2': r'$m_{miss}^2$ [GeV$^2$]',
    'el': r'$E_l$ [GeV]'
}
PLOT_RANGE = {
    'q2': [-0.4, 12.6],
    'mm2': [-2.0, 10.9],
    'el': [0.1, 2.65]
}
SMR_WTS = [
    ['', '_smr_pi', '_smr_k'],
    ['_no_smr', '_pi_smr', '_k_smr']
]


#######################
# Command line helper #
#######################

def parse_input():
    parser = ArgumentParser(
        description='plot fit variables w/ w/o decay-in-flight smearing.')

    parser.add_argument('-i', '--input', nargs='+',
                        help='specify main input ntuple.')

    parser.add_argument('-a', '--aux', nargs='+',
                        help='specify auxilliary ntuple containing misID weights.')

    parser.add_argument('-o', '--output', help='specify output folder.')

    parser.add_argument('-t', '--tree', help='specify tree name.',
                        default='tree')

    parser.add_argument('--bins', help='specify number of bins.', type=int,
                        default=300)

    parser.add_argument('-p', '--prefix', help='specify plot filename prefix.',
                        default='D0')

    parser.add_argument('--title', help='specify plot title.',
                        default=r'2016 misID, $D^0$')

    parser.add_argument('--ylabel', help='specify plot y label.',
                        default='Number of events')

    return parser.parse_args()


###########
# Helpers #
###########

def load_vars(ntps, tree, variables):
    variables_to_load = [v for v in variables
                         for n in [uproot.open(x) for x in ntps]
                         if v in n[tree]]
    for v in variables_to_load:
        print(f'  Loading {v}')
    return read_branches_dict(ntps, tree, variables_to_load)


def plot_comp(histos, legends, title, xlabel, ylabel, output_dir,
              output_filename,
              legend_loc='upper left', suffix='pdf', colors=DEFAULT_COLORS):
    data = [h[0] for h in histos]
    binspecs = [h[1] for h in histos]
    baselines = gen_histo_stacked_baseline(data)

    plotters = []
    for lbl, hist, bins, bot, clr in \
            zip(legends, data, binspecs, baselines, colors):
        add_args = ax_add_args_histo(lbl, clr, baseline=bot)
        plotters.append(
            lambda fig, ax, b=bins, h=hist+bot, add=add_args: plot_histo(
                b, h, add, figure=fig, axis=ax, show_legend=False))

    fig, ax, _ = plot_prepare(xlabel=xlabel, ylabel=ylabel, title=title,
                              show_legend=False)
    for p in plotters:
        p(fig, ax)

    handles, leg_lbls = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], leg_lbls[::-1], numpoints=1,
              loc=legend_loc, frameon='true')

    fig.savefig(f'{output_dir}/{output_filename}.{suffix}')


def plot_overall(histos, legends, title, xlabel, ylabel, output_dir,
                 output_filename,
                 legend_loc='upper left', suffix='pdf',
                 colors=DEFAULT_OVERALL_COLORS):
    data = [h[0] for h in histos]
    binspecs = [h[1] for h in histos]

    plotters = []
    for lbl, hist, bins, clr in \
            zip(legends, data, binspecs, colors):
        add_args = ax_add_args_step(lbl, clr)
        plotters.append(
            lambda fig, ax, b=bins, h=hist, add=add_args: plot_step(
                b, h, add, figure=fig, axis=ax, show_legend=False))

    fig, ax, _ = plot_prepare(xlabel=xlabel, ylabel=ylabel, title=title,
                              show_legend=False)
    for p in plotters:
        p(fig, ax)

    handles, leg_lbls = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], leg_lbls[::-1], numpoints=1,
              loc=legend_loc, frameon='true')

    fig.savefig(f'{output_dir}/{output_filename}.{suffix}')


########
# Main #
########

if __name__ == '__main__':
    mplhep.style.use('LHCb2')
    args = parse_input()
    ntp_main = args.input
    ntp_aux = args.aux

    all_vars = list(PLOT_VARS.keys()) + MISID_WTS + list(MISID_TAGS.keys()) + \
        [v + s for v in PLOT_VARS for s in SMR_WTS[0]] + \
        [w + s for w in MISID_WTS for s in SMR_WTS[1]]

    brs = load_vars(ntp_main, args.tree, all_vars)
    brs.update(load_vars(ntp_aux, args.tree, all_vars))

    for misid_wt in MISID_WTS:
        wt_tags = []
        wt_misid = []
        for tag_cut in MISID_TAGS:
            wt_tags.append(brs[tag_cut].astype(float))
            wt_misid.append(brs[misid_wt]*brs[tag_cut])

        wt_smr = []
        for smr_suf in SMR_WTS[1]:
            wt_smr.append(brs[misid_wt+smr_suf])

        for v in PLOT_VARS:
            histos_tags = []
            histos_misid = []
            histos_smr = []
            data_range = PLOT_RANGE[v]

            for w in wt_tags:
                histos_tags.append(gen_histo(
                    brs[v], args.bins, data_range=data_range, weights=w))

            for w in wt_misid:
                histos_misid.append(gen_histo(
                    brs[v], args.bins, data_range=data_range, weights=w))

            for w in wt_misid:
                tmp_histos = []
                for v_suffix, w_smr in zip(SMR_WTS[0], wt_smr):
                    tmp_histos.append(gen_histo(
                        brs[v+v_suffix], args.bins, data_range=data_range,
                        weights=w_smr*w))
                merged_tmp_histo = np.add.reduce([h[0] for h in tmp_histos])
                histos_smr.append((merged_tmp_histo, tmp_histos[0][1]))

            xlabel = PLOT_VARS[v]
            # raw data, just taggged
            plot_comp(histos_tags, MISID_TAGS.values(),
                 f'{args.title} (tags)',
                 xlabel, args.ylabel, args.output,
                 f'{args.prefix}_{misid_wt}_{v}_tags',
                 legend_loc=LEGEND_LOC[v])

            # unsmeared, w/ misID weights
            plot_comp(histos_misid, MISID_TAGS.values(),
                 f'{args.title} (misID weighted)',
                 xlabel, args.ylabel, args.output,
                 f'{args.prefix}_{misid_wt}_{v}_misid',
                 legend_loc=LEGEND_LOC[v])

            # smeared, w/ misID weights
            plot_comp(histos_smr, MISID_TAGS.values(),
                 f'{args.title} (misID weighted smeared)',
                 xlabel, args.ylabel, args.output,
                 f'{args.prefix}_{misid_wt}_{v}_smr',
                 legend_loc=LEGEND_LOC[v])

            # shape comparison
            # NOTE: The overall number of event is NOT guaranteed to conserve,
            #       as after smearing, some of the events might fall out of the
            #       fit variable acceptance cut
            merged_histo_misid = (
                np.add.reduce([h[0] for h in histos_misid]),
                histos_misid[0][1]
            )
            merged_histo_smr = (
                np.add.reduce([h[0] for h in histos_smr]),
                histos_smr[0][1]
            )
            plot_overall([merged_histo_misid, merged_histo_smr],
                         ['unsmeared', 'smeared'],
                         f'{args.title} (misID unsmeared/smeared comparison)',
                         xlabel, args.ylabel, args.output,
                         f'{args.prefix}_{misid_wt}_{v}_comp',
                         legend_loc=LEGEND_LOC[v])
