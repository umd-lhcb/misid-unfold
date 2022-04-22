#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Thu Apr 21, 2022 at 08:09 PM -0400
#
# Description: build decay-in-flight histos

import numpy as np
import uproot
import mplhep

from argparse import ArgumentParser

from pyTuplingUtils.boolean.eval import BooleanEvaluator
from pyTuplingUtils.utils import gen_histo
from pyTuplingUtils.plot import plot_histo, ax_add_args_histo


################
# Configurable #
################

DATA_RANGE = {
    'k_smr': [
        ((0.85, 1.15), 100),
        ((0.85, 1.15), 100),
        ((0.95, 1.05), 100),
    ],
    'pi_smr': [
        ((0.55, 1.05), 100),
        ((0.55, 1.05), 100),
        ((0.57, 1.03), 100),
    ],
}

ADD_LABELS = ['x', 'y', 'z']

OUTPUT_NTP_NAME = 'dif.root'
PARTICLES = {
    'k_smr': 'K',
    'pi_smr': 'pi',
}
TREE_NAME = 'mytuple/DecayTree'

TRUE_P_BRS = ['TRUEP_X', 'TRUEP_Y', 'TRUEP_Z']
RECO_P_BRS = ['PX', 'PY', 'PZ']

CUTS = {
    'k_smr': 'abs(K_TRUEID) == 321',
    'pi_smr': 'abs(pi_TRUEID) == 211'
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

def plot(br, title, filename, nbins=40, plot_range=(0, 1)):
    histo, bins = gen_histo(br, bins=nbins, data_range=plot_range)
    plot_histo(bins, histo,
               ax_add_args_histo(label='ratio', color='cornflowerblue'),
               output = filename, title=title, show_legend=False)


########
# Main #
########

if __name__ == '__main__':
    mplhep.style.use('LHCb2')
    args = parse_input()
    output_ntp = uproot.recreate(f'{args.output}/{OUTPUT_NTP_NAME}')

    ntps = [args.input_K, args.input_pi]
    ptcls = list(PARTICLES.keys())
    for (idx, n) in enumerate(ntps):
        evaluator = BooleanEvaluator(n, TREE_NAME)
        ptcl = ptcls[idx]
        prefix = PARTICLES[ptcl]

        true_brs = [evaluator.eval(f'{prefix}_{b}') for b in TRUE_P_BRS]
        reco_brs = [evaluator.eval(f'{prefix}_{b}') for b in RECO_P_BRS]
        global_cut = evaluator.eval(CUTS[ptcl])
        add_cuts = [global_cut]
        output_brs = []
        output_tree = dict()

        for i in range(len(RECO_P_BRS)):
            # reco / true !!!
            ratio = np.true_divide(
                reco_brs[i], true_brs[i], out=np.zeros_like(reco_brs[i]),
                where=true_brs[i] != 0)
            add_cuts.append(ratio > DATA_RANGE[ptcl][i][0][0])
            add_cuts.append(ratio < DATA_RANGE[ptcl][i][0][1])
            output_brs.append(ratio)

            if args.plot:
                title = f'{prefix}: {RECO_P_BRS[i]} / {TRUE_P_BRS[i]}'
                filename = f'{args.output}/{prefix}_{RECO_P_BRS[i]}_{TRUE_P_BRS[i]}.png'
                plot_range = DATA_RANGE[ptcl][i][0]
                nbins = DATA_RANGE[ptcl][i][1]
                # apply cut
                plot(ratio[global_cut], title, filename, nbins, plot_range)


        for idx, br in enumerate(output_brs):
            br_name = f'{ptcl}_{ADD_LABELS[idx]}'
            output_tree[br_name] = br[np.logical_and.reduce(add_cuts)]

        # now write the tree
        output_ntp[ptcl] = output_tree
