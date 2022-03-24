#!/usr/bin/env python
#
# Description: tagged histogram builder (T)

from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
from pyTuplingUtils.boolean.eval import BooleanEvaluator


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
    return f'{name.replace("misid_", "")}Tag'


########
# Main #
########

if __name__ == '__main__':
    args = parse_input()
    config_dir_path = abs_dir(args.config)

    with open(args.config, 'r') as f:
        config = safe_load(f)

    for particle, files in config['input_ntps'].items():
        evaluator = BooleanEvaluator(
            *ntp_tree(files, dir_abs_path=config_dir_path))

        specie_cuts = dict()
        for sp, cut in config['tags'].items():
            cut_expr = f'{" & ".join(config["add_global_cuts"])} & {convert_boolean_expr(cut)}'
            print('specie {histo_name_gen(sp)} has the following cuts: {cut_expr}')
            cut = evaluator.eval(cut_expr)

            # Make sure the evaluator is aware of the new variable
            evaluator.transformer.known_symb[sp] = cut
            evaluator.transformer.cache[sp] = cut

            specie_cuts[histo_name_gen(sp)] = cut
