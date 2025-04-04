#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Thu Sep 29, 2022 at 12:33 AM -0400
#
# Description: efficiency histogram builder (E)

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
from ROOT.RDF import TH3DModel
from ROOT.std import vector


#######################
# Command line parser #
#######################


def parse_input():
    parser = ArgumentParser(description="efficiency histogram builder (E).")

    parser.add_argument("-c", "--config", required=True, help="specify YAML config.")

    parser.add_argument("-o", "--output", required=True, help="specify output dir.")

    parser.add_argument("-y", "--year", type=int, default=2016, help="specify year.")

    parser.add_argument(
        "--ctrl-sample", action="store_true", help="Use control sample uBDT cut."
    )

    return parser.parse_args()


#################
# Histo helpers #
#################
ONE_SIGMA = 0.682689492137

def div_with_confint(num, denom):
    if denom < 1:
        return 0, 0

    ratio = num / denom
    intv = proportion_confint(num, denom, method="beta", alpha=1-ONE_SIGMA)  # Clopper-Pearson
    # Use the larger error bar and pretend its a Gaussian
    err_bar = max([abs(x - ratio) for x in intv])
    return ratio, err_bar


def histo_builder(binning_scheme, df, cuts=None, name="eff"):
    bin_vars = []
    bin_nums = []
    bin_edges = []

    for var_name, spec in binning_scheme.items():
        nums = -1
        edges = vector("double")()
        for s in spec:
            nums = nums + 1
            edges.push_back(s)
        bin_vars.append(var_name)
        bin_nums.append(nums)
        bin_edges.append(edges)

    if cuts is not None:
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

    histo = df.Histo3D(h3_model, *bin_vars)  # this returns a pointer and is lazy
    histo.GetTitle()  # make it execute immediately
    return histo


def compute_efficiency(histo_all, histo_passed, proj_p=False):
    if proj_p:
        # Project 3D histogram onto P axis
        # Meant for e conditional efficiency, where statistics is too low
        histo_passed_p = histo_passed.ProjectionX()
        histo_all_p = histo_all.ProjectionX()

        xbins = histo_all_p.GetNbinsX()

        histo_eff = histo_passed_p.Clone()
        for i in range(0, xbins + 2):
            n_all = histo_all_p.GetBinContent(i)
            n_passed = histo_passed_p.GetBinContent(i)

            eff, err = div_with_confint(n_passed, n_all)
            histo_eff.SetBinContent(i, eff)
            histo_eff.SetBinError(i, err)
    else:
        xbins = histo_all.GetNbinsX()
        ybins = histo_all.GetNbinsY()
        zbins = histo_all.GetNbinsZ()

        histo_eff = histo_passed.Clone()
        # take over/underflow bin into account
        for i, j, k in product(
            range(0, xbins + 2), range(0, ybins + 2), range(0, zbins + 2)
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


########
# Main #
########

if __name__ == "__main__":
    args = parse_input()
    config_dir_path = abs_dir(args.config)
    makedirs(args.output, exist_ok=True)

    with open(args.config, "r") as f:
        config = safe_load(f)
    binning_scheme = config["binning"]
    pid_config = config["local_pid_config"][args.year]
    tagged_config = config["tags"]

    for p_true, subconfig in pid_config.items():
        input_ntps, tree = ntp_tree(subconfig["samples"], config_dir_path)
        chain = TChain(tree)
        for n in input_ntps:
            chain.Add(n)
        df = RDataFrame(chain)

        # handle 'tags'
        if "tags" in subconfig:
            print(f"Handling tagged for {p_true}")
            global_cut = subconfig["tags"]["cut"].replace("&", "&&")
            histo_all = histo_builder(binning_scheme, df, global_cut, "all")
            print(f"  # of event passing 'cut': {histo_all.GetEntries()}")
            print(f"  cut: {global_cut}")
            for p_tag, pid_cut in tagged_config.items():
                print(f"    Handling {p_tag}")
                cuts = global_cut + " && " + pid_cut.replace("&", "&&")
                df = df.Define(p_tag, cuts)  # this is to make expr '!pi' work
                histo_passed = histo_builder(binning_scheme, df, p_tag, "passed")
                print(f"    pid_cut: {cuts}")
                print(f"    # of event passing 'pid_cut': {histo_passed.GetEntries()}")

                ntp_name = f"{p_true}TrueTo{p_tag.capitalize()}Tag.root"
                ntp_out = TFile(f"{args.output}/{ntp_name}", "recreate")

                eff = compute_efficiency(histo_all.GetPtr(), histo_passed.GetPtr())
                eff.SetName("eff")
                eff.Write()
                histo_all.Write()
                histo_passed.Write()

        # handle 'tags_addon'
        if "tags_addon" in subconfig:
            for p_tag, suffixes in subconfig["tags_addon"].items():
                print(f"Handling ad-hoc tagged for {p_true}")
                for suf in suffixes:
                    subsubconfig = subconfig["tags_addon"][p_tag][suf]
                    ntp_name = f"{p_true}TrueTo{p_tag.capitalize()}Tag_{suf}.root"
                    print(f"  Generating {ntp_name}")

                    global_cut = subsubconfig["cut"].replace("&", "&&")
                    histo_all = histo_builder(binning_scheme, df, global_cut, "all")
                    print(f"  cut: {global_cut}")
                    print(f"  # of event passing 'cut': {histo_all.GetEntries()}")

                    if args.ctrl_sample:
                        cuts = global_cut + " && " + subsubconfig["pid_cut"]["misid_ctrl"].replace("&", "&&")
                    else:
                        cuts = global_cut + " && " + subsubconfig["pid_cut"]["default"].replace("&", "&&")
                    histo_passed = histo_builder(binning_scheme, df, cuts, "passed")
                    print(f"  pid_cut: {cuts}")
                    print(
                        f"  # of event passing 'pid_cut': {histo_passed.GetEntries()}"
                    )

                    ntp_out = TFile(f"{args.output}/{ntp_name}", "recreate")

                    projection_p = subsubconfig["projection_p"] if "projection_p" in subsubconfig else False
                    eff = compute_efficiency(histo_all.GetPtr(), histo_passed.GetPtr(), projection_p)
                    eff.SetName("eff")
                    eff.Write()
                    histo_all.Write()
                    histo_passed.Write()
