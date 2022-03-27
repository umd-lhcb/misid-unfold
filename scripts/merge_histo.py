#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sun Mar 27, 2022 at 06:28 PM -0400
#
# Description: histogram merger (M)

from argparse import ArgumentParser
from pathlib import Path
from os import makedirs
from yaml import safe_load
from glob import glob

import uproot


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


def merge_true_to_tag(output_ntp, path_prefix, paths):
    for p in [f'{path_prefix}/{i}' for i in paths]:
        for n in glob(p):
            input_ntp = uproot.open(n)

            output_ntp[n.replace('.root', '')] = input_ntp['eff']


def merge_extra(output_ntp, path_prefix, config):
    for path in config:
        input_ntp = uproot.open(f'{path_prefix}/{path}')
        for src, tgt in config[path]:
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

    for mode, args in config['input_histos']:
        print(f'Merging {mode}...')
        KNOWN_MERGERS[mode](output_ntp, path_prefix, config)
