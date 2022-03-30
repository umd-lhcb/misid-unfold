#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Tue Mar 29, 2022 at 09:40 PM -0400
#
# Description: histogram merger (M)

from argparse import ArgumentParser
from pathlib import Path
from os import makedirs
from os.path import basename
from yaml import safe_load
from glob import glob

import uproot
import numpy as np


################
# Configurable #
################

HISTO_NAME = 'merged.root'


#######################
# Command line parser #
#######################

def parse_input():
    parser = ArgumentParser(description='histogram merger (M).')

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


def merge_true_to_tag(output_ntp, path_prefix, path):
    for n in glob(f'{path_prefix}/{path}'):
        input_ntp = uproot.open(n)
        histo = list(input_ntp['eff'].to_numpy())

        if 'MuTag' in n:
            print('Special treatment for MuTag efficiency...')
            aux_ntp = uproot.open(n.replace('MuTag', 'GlobalTag'))
            aux_histo = aux_ntp['eff'].to_numpy()
            # We need to divide out the efficiency(!isMuon)
            histo[0] = histo[0] / aux_histo[0]

        histo[0][np.isnan(histo[0])] = 0  # no need to keep nan
        output_ntp[Path(n).stem] = histo


def merge_extra(output_ntp, path_prefix, spec):
    for path in spec:
        input_ntp = uproot.open(f'{path_prefix}/{path}')
        for src, tgt in spec[path].items():
            output_ntp[tgt] = input_ntp[src]


########
# Main #
########

KNOWN_MERGERS = {
    'true_to_tag': merge_true_to_tag,
    'extra': merge_extra,
}


if __name__ == '__main__':
    args = parse_input()
    path_prefix = abs_dir(args.config)
    makedirs(args.output, exist_ok=True)
    output_ntp = uproot.recreate(f'{args.output}/{HISTO_NAME}')

    with open(args.config, 'r') as f:
        config = safe_load(f)

    for mode, args in config['input_histos'].items():
        print(f'Merging {mode}...')
        KNOWN_MERGERS[mode](output_ntp, path_prefix, args)
