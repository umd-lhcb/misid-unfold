#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Tue Apr 05, 2022 at 01:48 AM -0400
#
# Description: pidcalib2 wrapper (P)

import json
import re

from argparse import ArgumentParser
from os import system, chdir, makedirs
from os.path import basename, dirname, abspath
from glob import glob
from collections import namedtuple
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

    parser.add_argument('-y', '--year', default='2016',
                        help='specify year (format: 20XX).')

    return parser.parse_args()


###########
# Helpers #
###########

PidDirective = namedtuple(
    'PidDirective',
    [
        'sample_name', 'sample_file', 'year', 'polarity', 'bin_vars',
        'cut_arr', 'pid_cut_arr',
        'folder_name', 'pkl_names'
    ])


def dump_binning(yml_bins, samples, output):
    with open(output, 'w') as f:
        json.dump(
            {s: {BINNING_ALIAS[bin_name](s): bin_range
                 for bin_name, bin_range in yml_bins.items()}
             for s in samples}, f)


def run_cmd(cmd, dry_run=False):
    if dry_run:
        print(cmd)
    else:
        ret_val = system(cmd)
        if ret_val:
            raise ValueError('Command execution failed.')


def true_to_tag_gen(directive, dry_run=False, debug=False, polarity='down'):
    bin_vars = ' ' + ' '.join(f'--bin-var {b}' for b in directive.bin_vars)
    cuts = ''
    for gc, pc, nm in \
            zip(directive.cut_arr, directive.pid_cut_arr, directive.pkl_names):
        cuts += f' --cut "{gc}" --pid-cut "{pc}" --pkl-name {nm}.pkl'

    cmd = fr'''lb-conda pidcalib {CURR_DIR}/make_eff_histo_mod.py \
    --output-dir {directive.folder_name} \
    --sample {directive.sample_file} --magnet {directive.polarity} \
    --particle {directive.sample_name} \
    --binning-file ./tmp/{JSON_BIN_FILENAME}'''

    cmd += cuts
    cmd += bin_vars
    if debug:
        cmd += ' --max-files 3'  # debug only

    run_cmd(cmd, dry_run)

    if dry_run:
        return 1

    # Convert pkl -> root, rename and relocate
    for pkl in glob(f'{directive.folder_name}/*.pkl'):
        run_cmd(f'lb-conda pidcalib pidcalib2.pklhisto2root "{pkl}"')
    for ntp in glob(f'{directive.folder_name}/*.root'):
        run_cmd(f'cp "{ntp}" ./{basename(ntp)}')

    return 0


##############################
# Helpers for cut generation #
##############################

def cut_replacement(tagged):
    result = dict()

    for p, cut in tagged.items():
        cut = re.sub(r'!(\w+)', r'\1 == 0.0', cut)  # NOTE: workaround for negate expression
        for p_else in tagged:
            cut = re.sub(rf'\b{p_else}\b', f"({tagged[p_else]})", cut)
        result[p] = cut

    return result


def true_to_tag_directive_gen(config, year, output_folder, polarity='down'):
    result = []

    for p_true, sample_name in config['pidcalib_config']['samples'].items():
        cut_arr = []
        pid_cut_arr = []
        pkl_names = []
        bin_vars = [BINNING_ALIAS[b](sample_name) for b in config['binning']]

        folder_name = f'{output_folder}/{p_true}To-{year}'
        sample_file = SAMPLE_ALIAS(sample_name) + year[2:]

        # handle nominal tags first
        for p_tag, pid_cut in config['tags'].items():
            cut_arr.append(config['pidcalib_config']['tags']['cut'])
            pid_cut_arr.append(pid_cut)
            pkl_names.append(f'{p_true}To{p_tag.capitalize()}')

        # now handle ad-hoc tags
        for p_tag, subconfig in config['pidcalib_config']['tags_addon'].items():
            for sub_tag in subconfig:
                cut_arr.append(subconfig[sub_tag]['cut'])
                pid_cut_arr.append(subconfig[sub_tag]['pid_cut'])
                pkl_names.append(f'{p_true}To{p_tag.capitalize()}_{sub_tag}')

        result.append(PidDirective(
            sample_name, sample_file, year, polarity, bin_vars,
            cut_arr, pid_cut_arr, folder_name, pkl_names))

    return result


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

    # Dump custom binning schema (a JSON file, consumed by pidcalib later)
    dump_binning(
        config['binning'], config['pidcalib_config']['samples'].values(),
        f'./tmp/{JSON_BIN_FILENAME}')

    # Generate efficiency histograms with pidcalib2
    config['tags'] = cut_replacement(config['tags'])
    directives = true_to_tag_directive_gen(config, args.year, args.output)

    # in case of a dry run
    if args.dry_run:
        for d in directives:
            print(d.sample_name)
            for f in d._fields:
                if f != 'sample_name':
                    val = getattr(d, f)
                    if isinstance(val, str):
                        print(f'  {f}: {val}')
                    else:
                        print(f'  {f}:')
                        for i in val:
                            print(f'    {i}')

    for d in directives:
        true_to_tag_gen(d, dry_run=args.dry_run, debug=args.debug)
