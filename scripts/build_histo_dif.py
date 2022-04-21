#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Thu Apr 21, 2022 at 05:19 PM -0400
#
# Description: build decay-in-flight histos

import numpy as np
import uproot

from argparse import ArgumentParser

from pyTuplingUtils.boolean.eval import BooleanEvaluator
from pyTuplingUtils.utils import gen_histo
from pyTuplingUtils.plot import plot_histo, ax_add_args_histo


################
# Configurable #
################

HISTO_NAME = 'dif.root'
PARTICLES = {
    'kSmearing': 'K',
    'piSmearing': 'pi',
}
TREE_NAME = 'mytuple/DecayTree'

TRUE_P_BRS = ['TRUEP_X', 'TRUEP_Y', 'TRUEP_Z']
RECO_P_BRS = ['PX', 'PY', 'PZ']

CUTS = {
    'kSmearing': 'abs(K_TRUEID) == 321',
    'piSmearing': 'abs(pi_TRUEID) == 211'
}


#######################
# Command line parser #
#######################

def parse_input():
    parser = ArgumentParser(description='build decay-in-flight histos.')

    parser.add_argument('-k', '--input-K', required=True,
                        help='specify Kaon input ntuple.')

    parser.add_argument('-p', '--input-pi', required=True,
                        help='specify pion input ntuple.')

    parser.add_argument('-o', '--output', required=True,
                        help='specify output dir.')

    parser.add_argument('--plot', action='store_true',
                        help='plot momentum ratios, for finding binnings.')

    return parser.parse_args()


###########
# Helpers #
###########

def plot(br, title, filename, nbins=40):
    histo, bins = gen_histo(br, bins=nbins)
    plot_histo(bins, histo,
               ax_add_args_histo(label='ratio', color='cornflowerblue'),
               output = filename, title=title)


########
# Main #
########

if __name__ == '__main__':
    args = parse_input()

    ntps = [args.input_K, args.input_pi]
    histos = list(PARTICLES.keys())
    for (idx, n) in enumerate(ntps):
        evaluator = BooleanEvaluator(n, TREE_NAME)
        histo = histos[idx]
        particle = PARTICLES[histo]

        cut = evaluator.eval(CUTS[histo])
        true_brs = [evaluator.eval(f'{particle}_{b}') for b in TRUE_P_BRS]
        reco_brs = [evaluator.eval(f'{particle}_{b}') for b in RECO_P_BRS]

        for i in range(len(RECO_P_BRS)):
            # reco / true !!!
            ratio = np.true_divide(
                reco_brs[i], true_brs[i], out=np.zeros_like(reco_brs[i]),
                where=true_brs[i] != 0)
            # apply cut
            ratio = ratio[cut]

            if args.plot:
                title = f'{particle}: {RECO_P_BRS[i]} / {TRUE_P_BRS[i]}'
                filename = f'{args.output}/{particle}_{RECO_P_BRS[i]}_{TRUE_P_BRS[i]}.png'
                plot(ratio, title, filename)
