#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sat Mar 26, 2022 at 08:41 PM -0400
#
# Description: pidcalib2 wrapper (P)

import json

from argparse import ArgumentParser
from os import system, chdir, makedirs
from glob import glob
from yaml import safe_load


################
# Configurable #
################

JSON_BIN_FILENAME = 'binning.json'
SAMPLE_ALIAS = lambda p: 'Electron' if p == 'e' else 'Turbo'
REPLACEMENT_RULES = {
    '&&': '&',
    '||': '|',
    'mu_': '',
    'ProbNNpi': 'MC15TuneV1_ProbNNpi',
    'ProbNNghost': 'MC15TuneV1_ProbNNghost',
    'PIDK': 'DLLK',
    'PIDp': 'DLLp',
    'PIDe': 'DLLe'
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

    parser.add_argument('-y', '--year', default='16',
                        help='specify year (only the last 2 digits).')

    return parser.parse_args()


###########
# Helpers #
###########

def dump_binning(yml_bins, bin_alias, particles, output):
    binning = {bin_alias[k]: v for k, v in yml_bins.items()}
    with open(output, 'w') as f:
        return json.dump({p: binning for p in particles}, f)


def run_cmd(cmd, debug=False):
    if debug:
        print(cmd)
    else:
        system(cmd)


def pidcalib_gen(part_true, part_tag, part_sample, global_cut, pid_cut,
                 year, output_folder,
                 debug=False, polarity='down'):
    folder_name = f'{part_true}-to-{part_tag}-{year}'
    cmd = fr'''lb-conda pidcalib pidcalib2.make_eff_hists \
    --output-dir {output_folder}/{folder_name} \
    --sample {SAMPLE_ALIAS(part_true)}{year} --magnet {polarity} \
    --particle {part_sample} \
    --pid-cut "{pid_cut}" \
    --cut "{global_cut}" \
    --bin-var Brunel_P --bin-var Brunel_ETA --bin-var nTracks_Brunel \
    --binning-file ./tmp/{JSON_BIN_FILENAME}
    '''

    run_cmd(cmd, debug)

    # Convert pkl -> root, rename and relocate
    pkl = glob(f'{output_folder}/{folder_name}/*.pkl')
    pkl = 'fake.pkl' if not pkl else pkl[0]  # for debugging
    run_cmd(f'lb-conda pidcalib2.pklhisto2root {pkl}', debug=debug)

    ntp = glob(f'{output_folder}/{folder_name}/*.root')
    ntp = 'fake.root' if not ntp else ntp[0]  # for debugging
    run_cmd(f'cp {ntp} {folder_name}.root', debug=debug)


##############################
# Helpers for cut generation #
##############################

def cut_replacement(tagged_cuts):
    result = dict()

    for p, cut in tagged_cuts.items():
        for p_else in tagged_cuts:
            cut = cut.replace(p_else, f"({tagged_cuts[p_else]})")

        for src, tgt in REPLACEMENT_RULES.items():
            cut = cut.replace(src, tgt)

        result[p.replace('misid_', '')] = cut

    return result


def gen_pidcalib_sample_directive(config, year, output_folder, addon=['mu']):
    result = []
    tagged = cut_replacement(config['tags'])

    for p_true in list(tagged)+addon:
        if p_true not in config['particle_alias']:
            continue  # we don't use pidcalib2 for ghost

        for p_tag, pid_cut in tagged.items():
            global_cut = config['global_cuts']['mu'] if p_true == 'mu' else \
                config['global_cuts']['non_mu']
            result.append([
                p_true, p_tag, config['particle_alias'][p_true],
                global_cut, pid_cut, year, output_folder
            ])

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

    # Find particle names known to pidcalib
    particles_raw = [i.replace('misid_', '') for i in config['tags']] + ['mu']
    particles = [config['particle_alias'][i]
                 if i in config['particle_alias'] else i for i in particles_raw]

    # Dump custom binning schema (a JSON file, consumed by pidcalib later)
    dump_binning(config['binning'], config['binning_alias']['pidcalib'],
                 particles, f'./tmp/{JSON_BIN_FILENAME}')

    # Generate efficiency histograms with pidcalib2
    sample_directives = gen_pidcalib_sample_directive(
        config, args.year, 'raw_histos')
    for d in sample_directives:
        #  pidcalib_gen(*d, debug=True)
        pidcalib_gen(*d)
