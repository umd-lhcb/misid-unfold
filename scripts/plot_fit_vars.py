#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Tue Apr 19, 2022 at 03:55 PM -0400
#
# Description: plot fit variables w/ w/o decay-in-flight smearing

import uproot

from argparse import ArgumentParser


################
# Configurable #
################

PLOT_VARS = ['q2', 'mmiss2', 'el']
MISID_WTS = ['wmis_norm']
MISID_TAGS = {
    r'$K$ tag': 'is_misid_kTagToMuTag',
    r'$\pi$ tag': 'is_misid_piTagToMuTag',
    r'$p$ tag': 'is_misid_pTagToMuTag',
    r'$e$ tag': 'is_misid_eTagToMuTag',
    r'ghost tag': 'is_misid_gTagToMuTag',
}


#######################
# Command line helper #
#######################

def parse_input():
    parser = ArgumentParser(
        description='plot fit variables w/ w/o decay-in-flight smearing.')

    parser.add_argument('-i', '--input', help='specify main input ntuple.')

    parser.add_argument('-a', '--aux', help='specify auxilliary ntuple containing misID weights.')

    parser.add_argument('-o', '--output', help='specify output folder.')

    parser.add_argument('-t', '--tree', help='specify tree name.',
                        default='tree')

    return parser.parse_args()


###########
# Helpers #
###########

def load_vars(ntp, tree, variables):
    variables_to_load = [v for v in variables if v in ntp[tree]]
    print(variables_to_load)


########
# Main #
########

if __name__ == '__main__':
    args = parse_input()
    ntp_main = uproot.open(args.input)

    all_vars = PLOT_VARS + MISID_WTS + list(MISID_TAGS.values())

    print(all_vars)
    print(PLOT_VARS)
    load_vars(ntp_main, args.tree, all)
