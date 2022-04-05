#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Tue Apr 05, 2022 at 12:15 PM -0400
#
# Description: histogram merger (M)

from argparse import ArgumentParser
from pathlib import Path
from os import makedirs
from yaml import safe_load

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

    parser.add_argument('-y', '--year', default='2016',
                        help='specify year.')

    return parser.parse_args()


###########
# Helpers #
###########

def abs_dir(path):
    return str(Path(path).parent.absolute())


def merge_true_to_tag(output_ntp, path_prefix, path, config):
    ptcl_true = list(config['pidcalib_config']['samples'])
    ptcl_tag = list(config['tags'])

    # Handle normal trueToTag first
    for p_true in ptcl_true:
        for p_tag in ptcl_tag:
            histo_name = f'{p_true}To{p_tag.capitalize()}'
            input_ntp = uproot.open(f'{path_prefix}/{path}/{histo_name}.root')
            histo = list(input_ntp['eff'].to_numpy())

            histo[0][np.isnan(histo[0])] = 0  # no need to keep nan
            output_ntp[histo_name] = histo

    # Handle additional ntuples
    for p_true in ptcl_true:
        for p_addon, suffixes \
                in config['pidcalib_config']['tags_addon'].items():
            histo_name_base = f'{p_true}To{p_addon.capitalize()}'
            aux_histos = []

            for suf in suffixes:
                histo_name = histo_name_base + '_' + suf
                input_ntp = uproot.open(f'{path_prefix}/{path}/{histo_name}.root')
                histo = list(input_ntp['eff'].to_numpy())

                aux_histos.append(histo)

            histo_ratio = aux_histos[0][0] / aux_histos[1][0]

            histo_ratio[np.isnan(histo_ratio)] = 0
            output_ntp[histo_name_base] = tuple([histo_ratio]+aux_histos[0][1:])


def merge_extra(output_ntp, path_prefix, spec, config):
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

    for mode, args in config['input_histos'][int(args.year)].items():
        print(f'Merging {mode}...')
        KNOWN_MERGERS[mode](output_ntp, path_prefix, args, config)
