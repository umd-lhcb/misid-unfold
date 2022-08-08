#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Mon Aug 08, 2022 at 01:45 AM -0400
#
# Description: efficiency histogram builder (E)

import numpy as np

from argparse import ArgumentParser
from pathlib import Path
from os import makedirs
from itertools import product
from yaml import safe_load
from statsmodels.stats.proportion import proportion_confint

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True  # Don't hijack argparse!
ROOT.PyConfig.DisableRootLogon = True  # Don't read .rootlogon.py

from ROOT import TChain, TFile, RDataFrame
from ROOT.RDF import TH3DModel, TH2DModel, TH1DModel
from ROOT.std import vector


################
# Configurable #
################

HISTO_NAME = "eff.root"


#######################
# Command line parser #
#######################


def parse_input():
    parser = ArgumentParser(description="tagged histogram builder (T).")

    parser.add_argument("-c", "--config", required=True, help="specify YAML config.")

    parser.add_argument("-o", "--output", required=True, help="specify output dir.")

    parser.add_argument("-y", "--year", default="2016", help="specify year.")

    return parser.parse_args()


#################
# Histo helpers #
#################


def div_with_confint(num, denom):
    if denom < 1:
        return 0, 0

    ratio = num / denom
    intv = proportion_confint(num, denom, method="beta", alpha=0.32)  # Clopper-Pearson
    # Use the larger error bar and pretend its a Gaussian
    err_bar = max([abs(x - ratio) for x in intv])
    return ratio, err_bar


def histo_builder(binning_scheme, ntps, cuts=None, tree="tree", name="eff"):
    bin_vars = []
    bin_nums = []
    bin_edges = []

    for var_name, spec in binning_scheme.items():
        bin_vars.append(var_name)
        nums = 0
        edges = vector("double")()
        for s in spec:
            nums = nums + 1
            edges.push_back(s)
        bin_nums.append(nums)
        bin_edges.append(edges)

    chain = TChain(tree)
    for n in ntps:
        chain.Add(n)

    df = RDataFrame(chain)
    if cuts:
        df = df.Filter(cuts)

    h3_model = TH3DModel(
        name,
        name,
        bin_nums[0],
        bin_edges[0].data(),
        bin_nums[1],
        bin_edges[1].data(),
        bin_nums[2],
        bin_edges[2].data(),
    )
    return df.Histo3D(h3_model, *bin_vars)  # this is a smart pointer


def compute_efficiency(histo_all, histo_passed):
    xbins = histo_all.GetNbinsX()
    ybins = histo_all.GetNbinsY()
    zbins = histo_all.GetNbinsZ()

    histo_eff = histo_passed.Clone()
    for i, j, k in product(
        range(1, xbins + 1), range(1, ybins + 1), range(1, zbins + 1)
    ):
        idx = histo_all.GetBin(i, j, k)
        n_all = histo_all.GetBinContent(idx)
        n_passed = histo_passed.GetBinContent(idx)

        eff, err = div_with_confint(n_passed, n_all)
        histo_eff.SetBinContent(idx, eff)
        histo_eff.SetBinError(idx, err)

    return histo_eff


###############
# Path helper #
###############


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
