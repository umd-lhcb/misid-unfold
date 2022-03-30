#!/usr/bin/env python
#  Author: Yipeng Sun
#  License: GPLv3
#  Last Change: Wed Mar 30, 2022 at 02:52 PM -0400
#
# Description: histogram plotter (for this project)

import re
import uproot
import mplhep
import numpy as np

from argparse import ArgumentParser
from pyTuplingUtils.utils import gen_histo_stacked_baseline
from pyTuplingUtils.plot import (
    plot_top, plot_histo, ax_add_args_histo
)


################
# Configurable #
################

DEFAULT_COLORS = [
    'mediumseagreen',
    'lightsteelblue',
    'mediumpurple',
    'salmon',
    'palegoldenrod',
]

KNOWN_PARTICLES = {
    'd0': r'$D^0$',
    'e': r'$e$',
    'pi': r'$\pi$',
    'k': r'$K$',
    'g': 'ghost',
    'mu': r'$\mu$',
}

GET_PARTICLE = lambda x: KNOWN_PARTICLES[x.lower()] \
    if x.lower() in KNOWN_PARTICLES else x


#######################
# Command line helper #
#######################

def parse_input():
    parser = ArgumentParser(description='histogram plotter (for this project).')

    parser.add_argument('-i', '--input', help='specify input ntuple.')

    parser.add_argument('-o', '--output', help='specify output folder.')

    parser.add_argument('-p', '--prefix', action='append', default=['D0'],
                        help='specify histo prefixes to plot.')

    parser.add_argument('-s', '--suffix', default='True',
                        help='specify histo suffix to plot.')

    parser.add_argument('-v', '--vars', nargs='+',
                        default=['p', 'eta', 'ntracks'],
                        help='specify histo binning variables.')

    parser.add_argument('-l', '--labels', nargs='+',
                        default=[r'$p$ [MeV]', r'$\eta$', r'nTracks'],
                        help='specify histo binning variables in displayed format.')

    return parser.parse_args()


###########
# Helpers #
###########

def find_key(name, list_of_keys, sep='__'):
    for k in list_of_keys:
        if k + sep in name:
            return True
    return False


def name_cleanup(name):
    return name.split(';')[0]


def prefix_gen(histo_spec):
    raw_name = list(histo_spec)[0]
    ptcl, descr = raw_name.split('__')
    mod = ''
    if 'true' in descr.lower():
        mod = 'true'
    elif 'tag' in descr.lower():
        mod = 'tag'
    return f'{ptcl}_{mod}'


def label_gen(name, known=KNOWN_PARTICLES):
    name = name_cleanup(name).split('__')[1]
    regex = re.search(r'(\w+)(True|Tag)', name)
    return f'{GET_PARTICLE(regex.group(1))} {regex.group(2)}'


def title_gen(name, known=KNOWN_PARTICLES):
    ptcl, descr = name.split('_')
    return f'{GET_PARTICLE(ptcl)} {descr.lower()}'


########
# Plot #
########

def plot(histo_spec, bin_vars, bin_names, output_dir,
         colors=DEFAULT_COLORS, suffix='png'):
    prefix = prefix_gen(histo_spec)

    for idx, v in enumerate(bin_vars):
        ordering = sorted(histo_spec)
        labels = [label_gen(i) for i in ordering]
        histos = [np.sum(histo_spec[i][0],
                         axis=tuple(x for x in range(3) if x!= idx))
                  for i in ordering]
        binspecs = [histo_spec[i][1+idx] for i in ordering]
        baselines = gen_histo_stacked_baseline(histos)

        plotters = []
        for lbl, hist, bins, bot, clr in \
                zip(labels, histos, binspecs, baselines, colors):
            add_args = ax_add_args_histo(lbl, clr, baseline=bot)
            plotters.append(
                lambda fig, ax, b=bins, h=hist+bot, add=add_args: plot_histo(
                    b, h, add, figure=fig, axis=ax, show_legend=False))

        plot_top(plotters, f'{output_dir}/{prefix}_{v}.{suffix}',
                 xlabel=bin_names[idx], title=title_gen(prefix),
                 legend_add_args={
                     'numpoints': 1, 'loc': 'best', 'frameon': 'true'})


########
# Main #
########

if __name__ == '__main__':
    mplhep.style.use('LHCb2')
    args = parse_input()
    ntp = uproot.open(args.input)

    # group histograms
    for pref in args.prefix:
        histos = {name_cleanup(k): ntp[k].to_numpy() for k in ntp
                  if find_key(k, args.prefix) and
                  name_cleanup(k).endswith(args.suffix)}

        plot(histos, args.vars, args.labels, args.output)
