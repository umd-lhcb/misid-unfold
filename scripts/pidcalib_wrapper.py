#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sat Mar 26, 2022 at 04:57 AM -0400
#
# Description: pidcalib2 wrapper (P)

import json

from argparse import ArgumentParser
from yaml import safe_load


################
# Configurable #
################

JSON_BIN_FILENAME = 'binning.json'


#######################
# Command line parser #
#######################

def parse_input():
    parser = ArgumentParser(description='pidcalib2 wrapper (P).')

    parser.add_argument('-c', '--config', required=True,
                        help='specify YAML config.')

    parser.add_argument('-o', '--output', required=True,
                        help='specify output dir.')

    parser.add_argument('-d', '--dry-run', action='store_true',
                        help='printout command to execute w/o actually executing.')

    return parser.parse_args()


###########
# Helpers #
###########

def get_particles(tags, alias, addons=['mu']):
    raw = [i.replace('misid_', '') for i in tags] + addons
    return [alias[i] if i in alias else i for i in raw]


def dump_binning(yml_bins, bin_alias, particles, output):
    binning = {bin_alias[k]: v for k, v in yml_bins.items()}
    with open(output, 'w') as f:
        return json.dump({p: binning for p in particles}, f)


########
# Main #
########

if __name__ == '__main__':
    args = parse_input()

    with open(args.config, 'r') as f:
        config = safe_load(f)

    # Dump custom binning schema (a JSON file, consumed by pidcalib later)
    particles = get_particles(config['tags'], config['particle_alias'])
    dump_binning(config['binning'], config['binning_alias']['pidcalib'],
                 particles, f'{args.output}/{JSON_BIN_FILENAME}')
