#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sun Mar 27, 2022 at 12:21 AM -0400
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

    parser.add_argument('-m', '--mode', default='true_to_tag',
                        help='specify working mode.')

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


def true_to_tag_gen(part_true, part_sample, part_tag_arr, global_cuts, pid_cuts,
                 year, output_folder,
                 debug=False, polarity='down'):
    cuts = ''
    for gc, pc in zip(global_cuts, pid_cuts):
        cuts += f' --cut "{gc}" --pid-cut "{pc}"'

    folder_name = f'{part_true}True-{year}'
    cmd = fr'''lb-conda pidcalib pidcalib2.make_eff_hists \
    --output-dir {output_folder}/{folder_name} \
    --sample {SAMPLE_ALIAS(part_true)}{year} --magnet {polarity} \
    --particle {part_sample} \
    --bin-var Brunel_P --bin-var Brunel_ETA --bin-var nTracks_Brunel \
    --binning-file ./tmp/{JSON_BIN_FILENAME}'''

    cmd += cuts
    cmd += ' --max-files 3'  # debug only

    run_cmd(cmd, debug)

    if debug:
        return 1

    # Convert pkl -> root, rename and relocate
    rename_rules = {cut.replace(' ', ''): part
                    for cut, part in zip(pid_cuts, part_tag_arr)}

    for pkl in glob(f'{output_folder}/{folder_name}/*.pkl'):
        run_cmd(f'lb-conda pidcalib2.pklhisto2root "{pkl}"')
        for ntp in glob(f'{output_folder}/{folder_name}/*.root'):
            for rule, p in rename_rules:
                if rule in ntp:
                    run_cmd(f'cp "{ntp}" {part_true}To{p.capitalize()}-{year}.root')

    return 0


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


def true_to_tag_directive_gen(config, year, output_folder, addon=['mu']):
    result = []
    tagged = cut_replacement(config['tags'])

    for p_true in list(tagged)+addon:
        if p_true not in config['particle_alias']:
            continue  # we don't use pidcalib2 for ghost

        global_cuts = []
        pid_cuts = []

        for pid_cut_addon in tagged.values():
            global_cuts.append(config['global_cuts']['kinematic'])

            pid_cut = config['global_cuts']['mu'] if p_true == 'mu' else \
                config['global_cuts']['non_mu']
            pid_cut += pid_cut_addon
            pid_cuts.append(pid_cut)

        result.append([
            p_true, config['particle_alias'][p_true], tagged,
            global_cuts, pid_cuts, year, output_folder
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
    if args.mode == 'true_to_tag':
        directives = true_to_tag_directive_gen(config, args.year, 'raw_histos')
        for d in directives:
            true_to_tag_gen(*d, debug=args.dry_run)
    else:
        print(f'unknown mode: {args.mode}')
