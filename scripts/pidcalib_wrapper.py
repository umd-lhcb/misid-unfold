#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Mon Apr 04, 2022 at 07:08 PM -0400
#
# Description: pidcalib2 wrapper (P)

import json
import re

from argparse import ArgumentParser
from os import system, chdir, makedirs
from os.path import basename, dirname, abspath
from glob import glob
from yaml import safe_load


################
# Configurable #
################

CURR_DIR = dirname(abspath(__file__))
JSON_BIN_FILENAME = 'binning.json'
SAMPLE_ALIAS = lambda p: 'Electron' if p == 'e_B_Jpsi' else 'Turbo'
BINNING_ALIAS = {
    'P': lambda sample: 'Brunel_P',
    'ETA': lambda sample: 'Brunel_ETA',
    'nTracks': lambda sample: \
        'nTracks' if sample == 'e_B_Jpsi' else 'nTracks_Brunel',
}


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

    parser.add_argument('-D', '--debug', action='store_true',
                        help='enable debug mode by looping over a small subset of pidcalib samples.')

    parser.add_argument('-y', '--year', default='16',
                        help='specify year (only the last 2 digits).')

    parser.add_argument('-m', '--mode', default='true_to_tag',
                        help='specify working mode.')

    return parser.parse_args()


###########
# Helpers #
###########

def bin_alias(config):
    # FIXME: Special treatement for electron samples
    def inner(var, part):
        if part == 'e_B_Jpsi' and var == 'nTracks':
            return 'nTracks'
        return config['binning_alias']['pidcalib'][var]
    return inner


def dump_binning(yml_bins, config, samples, output):
    alias = bin_alias(config)

    with open(output, 'w') as f:
        json.dump(
            {s: {alias(k, s): v for k, v in yml_bins.items()}
             for s in samples}, f)


def run_cmd(cmd, dry_run=False):
    if dry_run:
        print(cmd)
    else:
        ret_val = system(cmd)
        if ret_val:
            raise ValueError('Command execution failed.')


def true_to_tag_gen(part_true, part_sample, part_tag_arr, global_cuts, pid_cuts,
                    year, output_folder,
                    dry_run=False, debug=False, polarity='down'):
    cuts = ''
    for gc, pc, nm in zip(global_cuts, pid_cuts, part_tag_arr):
        cuts += f' --cut "{gc}" --pid-cut "{pc}" --pkl-name {part_true}TrueTo{nm.capitalize()}Tag.pkl'

    # FIXME: nTrack name
    ntracks = 'nTracks_Brunel' if part_true != 'e' else 'nTracks'

    folder_name = f'{part_true}True-{year}'
    cmd = fr'''lb-conda pidcalib {CURR_DIR}/make_eff_histo_mod.py \
    --output-dir {output_folder}/{folder_name} \
    --sample {SAMPLE_ALIAS(part_true)}{year} --magnet {polarity} \
    --particle {part_sample} \
    --bin-var Brunel_P --bin-var Brunel_ETA --bin-var {ntracks} \
    --binning-file ./tmp/{JSON_BIN_FILENAME}'''

    cmd += cuts
    if debug:
        cmd += ' --max-files 3'  # debug only

    run_cmd(cmd, dry_run)

    if dry_run:
        return 1

    # Convert pkl -> root, rename and relocate
    for pkl in glob(f'{output_folder}/{folder_name}/*.pkl'):
        run_cmd(f'lb-conda pidcalib pidcalib2.pklhisto2root "{pkl}"')
    for ntp in glob(f'{output_folder}/{folder_name}/*.root'):
        run_cmd(f'cp "{ntp}" ./{basename(ntp)}')

    return 0


##############################
# Helpers for cut generation #
##############################

def cut_replacement(tagged_cuts):
    result = dict()

    for p, cut in tagged_cuts.items():
        for p_else in tagged_cuts:
            cut = re.sub(rf'\b{p_else}\b', f"({tagged_cuts[p_else]})", cut)
        result[p] = cut

    return result


def true_to_tag_directive_gen(tagged, tagged_addon, year, output_folder):
    result = []

    for p_true in tagged:  # yes, it's the true particle
        if p_true not in config['particle_alias']:
            continue  # e.g. we don't use pidcalib2 for ghost

        global_cuts = []
        pid_cuts = []

        for pid_cut_addon in tagged.values():
            global_cuts.append(config['global_cuts']['kinematic'])

            pid_cut = config['global_cuts']['tags'] + '&' + pid_cut_addon
            pid_cuts.append(pid_cut)

        # now handle 'tags_addon'
        for cut in tagged_addon.values():
            global_cuts.append(config['global_cuts']['kinematic'])
            pid_cuts.append(cut)

        result.append([
            p_true, config['particle_alias'][p_true],
            list(tagged) + list(tagged_addon),
            global_cuts, pid_cuts, year, output_folder
        ])

    return result  # we keep 'tags_addon' as-is


########
# Main #
########

if __name__ == '__main__':
    args = parse_input()

    with open(args.config, 'r') as f:
        config = safe_load(f)

    # Go to output dir
    makedirs(args.output, exist_ok=True)
    chdir(args.output)
    makedirs('raw_histos', exist_ok=True)
    makedirs('tmp', exist_ok=True)

    # Find particle names known to pidcalib
    samples = config['pidcalib_samples']

    # Dump custom binning schema (a JSON file, consumed by pidcalib later)
    dump_binning(config['binning'], config, ptcl, f'./tmp/{JSON_BIN_FILENAME}')

    # Generate efficiency histograms with pidcalib2
    if args.mode == 'true_to_tag':
        directives = true_to_tag_directive_gen(
            ptcl_tagged, config['tags_addon'], args.year, 'raw_histos')
        for d in directives:
            true_to_tag_gen(*d, dry_run=args.dry_run, debug=args.debug)
    else:
        print(f'unknown mode: {args.mode}')
