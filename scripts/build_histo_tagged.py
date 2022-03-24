#!/usr/bin/env python
#
# Description: tagged histogram builder (T)

from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load


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


if __name__ == '__main__':
    args = parse_input()

    config_dir_path = abs_dir(args.config)

    print(config_dir_path)
