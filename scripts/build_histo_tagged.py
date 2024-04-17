#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Tue Apr 05, 2022 at 02:15 AM -0400
#
# Description: tagged histogram builder (T)

import numpy as np
import uproot

from argparse import ArgumentParser
from pathlib import Path
from os import makedirs
from yaml import safe_load
from pyTuplingUtils.boolean.eval import BooleanEvaluator


################
# Configurable #
################

HISTO_NAME = "tagged.root"


#######################
# Command line parser #
#######################


def parse_input():
    parser = ArgumentParser(description="tagged histogram builder (T).")

    parser.add_argument("-c", "--config", required=True, help="specify YAML config.")

    parser.add_argument("-o", "--output", required=True, help="specify output dir.")

    parser.add_argument("-y", "--year", default="2016", help="specify year.")

    parser.add_argument(
        "--ctrl-sample", action="store_true", help="Use control sample uBDT cut."
    )

    return parser.parse_args()


###########
# Helpers #
###########


def abs_dir(path):
    return str(Path(path).parent.absolute())


def ntp_tree(name, dir_abs_path=""):
    ntps = []
    path, tree = name.split(":")
    ntps.append(f"{dir_abs_path}/{path}")
    return ntps, tree


def histo_name_gen(name):
    return f"{name}Tag"


########
# Main #
########

if __name__ == "__main__":
    args = parse_input()
    config_dir_path = abs_dir(args.config)
    if (args.ctrl_sample):
        # Assuming .root extension
        HISTO_NAME = HISTO_NAME[0:-5] + "_misid_ctrl" + HISTO_NAME[-5:]
    makedirs(args.output, exist_ok=True)
    ntp = uproot.recreate(f"{args.output}/{HISTO_NAME}")

    with open(args.config, "r") as f:
        config = safe_load(f)

    for particle, subconfig in config["input_ntps"][int(args.year)].items():
        print(f"Working on {particle}...")
        evaluator = BooleanEvaluator(
            *ntp_tree(subconfig["files"], dir_abs_path=config_dir_path)
        )

        # load branches needed to build histos
        histo_brs = []
        for br_name in config["binning"]:
            histo_brs.append(evaluator.eval(br_name))

        global_cut = evaluator.eval(subconfig["cuts"])

        for sp, cut_expr in config["tags"].items():
            print(f"  specie {histo_name_gen(sp)} has the following cuts: {cut_expr}")
            cut = evaluator.eval(cut_expr)

            # Make sure the evaluator is aware of the new variable
            evaluator.transformer.known_symb[sp] = cut
            evaluator.transformer.cache[sp] = cut

            # Now build histograms
            ntp[f"{particle}__{histo_name_gen(sp)}"] = np.histogramdd(
                histo_brs,
                bins=list(config["binning"].values()),
                weights=(cut & global_cut),
            )
