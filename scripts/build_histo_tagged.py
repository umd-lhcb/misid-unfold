#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Tue Mar 29, 2022 at 10:14 PM -0400
#
# Description: tagged histogram builder (T)

import numpy as np
import uproot

from argparse import ArgumentParser
from pathlib import Path
from os import makedirs
from yaml import safe_load
from pyTuplingUtils.boolean.eval import BooleanEvaluator


################
# Configurable #
################

HISTO_NAME = 'tagged.root'


#######################
# Command line parser #
#######################

def parse_input():
    parser = ArgumentParser(description='tagged histogram builder (T).')

    parser.add_argument('-c', '--config', required=True,
                        help='specify YAML config.')

    parser.add_argument('-o', '--output', required=True,
                        help='specify output dir.')

    return parser.parse_args()


###########
# Helpers #
###########

def abs_dir(path):
    return str(Path(path).parent.absolute())


# Convert a C++ styled boolean expr to Python style
def convert_boolean_expr(expr):
    return expr.replace('&&', '&').replace('||', '|')


def ntp_tree(name, tree=None, dir_abs_path=''):
    ntps = []

    for n in name:
        path, tree_in_file = n.split(':')
        ntps.append(f'{dir_abs_path}/{path}')
        if not tree:
            tree = tree_in_file

    return ntps, tree


def histo_name_gen(name):
    return f'{name}Tag'


#######################
# Computation helpers #
#######################

COMPUTE_FUNC = {
    'ETA': lambda p, pz: 0.5*np.log((p+pz) / (p-pz))
}


########
# Main #
########

if __name__ == '__main__':
    args = parse_input()
    config_dir_path = abs_dir(args.config)
    makedirs(args.output, exist_ok=True)
    ntp = uproot.recreate(f'{args.output}/{HISTO_NAME}')

    with open(args.config, 'r') as f:
        config = safe_load(f)

    for particle, files in config['input_ntps'].items():
        print(f'Working on {particle}...')
        evaluator = BooleanEvaluator(
            *ntp_tree(files, dir_abs_path=config_dir_path))
        evaluator.transformer.known_func.update(COMPUTE_FUNC)

        # load branches needed to build histos
        histo_brs = []
        for br_name, var in config['binning_alias']['offline'].items():
            var = [var] if not isinstance(var, list) else var

            if br_name in COMPUTE_FUNC:
                histo_brs.append(evaluator.eval(f'{br_name}({",".join(var)})'))
            else:
                histo_brs.append(evaluator.eval(var[0]))

        global_cut = evaluator.eval(config["global_cuts"]["offline"])

        for sp, cut in config['tags'].items():
            cut_expr = convert_boolean_expr(cut)
            print(f'  specie {histo_name_gen(sp)} has the following cuts: {cut_expr}')
            cut = evaluator.eval(cut_expr)

            # Make sure the evaluator is aware of the new variable
            evaluator.transformer.known_symb[sp] = cut
            evaluator.transformer.cache[sp] = cut

            # Now build histograms
            ntp[f'{particle}__{histo_name_gen(sp)}'] = np.histogramdd(
                histo_brs, bins=list(config['binning'].values()),
                weights=(cut & global_cut)
            )
