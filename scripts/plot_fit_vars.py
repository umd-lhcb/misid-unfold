#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Tue Apr 19, 2022 at 09:40 PM -0400
#
# Description: plot fit variables w/ w/o decay-in-flight smearing

import mplhep
import uproot

from argparse import ArgumentParser
from pyTuplingUtils.utils import gen_histo_stacked_baseline, gen_histo
from pyTuplingUtils.plot import (
    plot_prepare, plot_histo, ax_add_args_histo
)


################
# Configurable #
################

DEFAULT_COLORS = ['#00429d', '#5585b7', '#8bd189', '#d15d5f', '#93003a']
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
    'mm2': r'$m_{miss}^2 [GeV$^2$]',
    'el': r'$E_l$ [GeV]'
}
PLOT_RANGE = {
    'q2': [-0.4, 12.6],
    'mm2': [-2.0, 10.9],
    'el': [0.1, 2.65]
}


#######################
# Command line helper #
#######################

def parse_input():
    parser = ArgumentParser(
        description='plot fit variables w/ w/o decay-in-flight smearing.')

    parser.add_argument('-i', '--input', help='specify main input ntuple.')

    parser.add_argument('-a', '--aux', help='specify auxilliary ntuple containing misID weights.')

    parser.add_argument('-o', '--output', help='specify output folder.')

    parser.add_argument('-t', '--tree', help='specify tree name.',
                        default='tree')

    parser.add_argument('--bins', help='specify number of bins.', type=int,
                        default=40)

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

def load_vars(ntp, tree, variables):
    variables_to_load = [v for v in variables if v in ntp[tree]]
    for v in variables_to_load:
        print(f'  Loading {v}')
    return ntp[tree].arrays(variables_to_load, library='np')


def plot(histos, legends, title, xlabel, ylabel, output_dir, output_filename,
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


########
# Main #
########

if __name__ == '__main__':
    mplhep.style.use('LHCb2')
    args = parse_input()
    ntp_main = uproot.open(args.input)
    ntp_aux = uproot.open(args.aux)

    all_vars = list(PLOT_VARS.keys()) + MISID_WTS + list(MISID_TAGS.keys())

    brs = load_vars(ntp_main, args.tree, all_vars)
    brs.update(load_vars(ntp_aux, args.tree, all_vars))

    for misid_wt in MISID_WTS:
        wt_tags = []
        wt_misid = []

        for tag_cut in MISID_TAGS:
            wt_tags.append(brs[tag_cut].astype(float))
            wt_misid.append(brs[misid_wt]*brs[tag_cut])

        for v in PLOT_VARS:
            histos_tags = []
            histos_misid = []
            data_range = PLOT_RANGE[v]

            for w in wt_tags:
                histos_tags.append(gen_histo(
                    brs[v], args.bins, data_range=data_range, weights=w))

            for w in wt_misid:
                histos_misid.append(gen_histo(
                    brs[v], args.bins, data_range=data_range, weights=w))

            xlabel = PLOT_VARS[v]
            plot(histos_tags, MISID_TAGS.values(),
                 f'{args.title} (tags)',
                 xlabel, args.ylabel, args.output,
                 f'{args.prefix}_{misid_wt}_{v}_tags',
                 legend_loc=LEGEND_LOC[v])
            plot(histos_misid, MISID_TAGS.values(),
                 f'{args.title} (misID weighted)',
                 xlabel, args.ylabel, args.output,
                 f'{args.prefix}_{misid_wt}_{v}_misid',
                 legend_loc=LEGEND_LOC[v])
