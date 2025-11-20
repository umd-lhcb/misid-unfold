#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Mon Sep 12, 2022 at 05:02 AM -0400
#
# Description: histogram merger (M)

import itertools
import numpy as np

from argparse import ArgumentParser
from pathlib import Path
from os import makedirs
from os.path import basename
from glob import glob
from yaml import safe_load
from scipy.special import erf, erfinv

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True  # Don't hijack argparse!
ROOT.PyConfig.DisableRootLogon = True  # Don't read .rootlogon.py


################
# Configurable #
################

HISTO_NAME = "merged.root"
BAD_ERROR_THRESH = 0.2


#######################
# Command line parser #
#######################


def parse_input():
    parser = ArgumentParser(description="histogram merger (M).")

    parser.add_argument("-c", "--config", required=True, help="specify YAML config.")

    parser.add_argument("-o", "--output", required=True, help="specify output dir.")

    return parser.parse_args()


#################
# Histo helpers #
#################


def get_bin_info(histo, axis="x"):
    nbins = getattr(histo, f"GetNbins{axis.upper()}")()
    axis = getattr(histo, f"Get{axis.upper()}axis")()
    bin_edges = ROOT.std.vector("double")()

    for i in range(1, nbins + 2):
        bin_edges.push_back(getattr(axis, "GetBinLowEdge")(i))

    return nbins, bin_edges


def get_axis_title(histo, axis="x"):
    axis = getattr(histo, f"Get{axis.upper()}axis")()
    return getattr(axis, "GetTitle")()


def prep_root_histo(name, histo_orig):
    histo = None
    histo_axis_nbins = None

    if "TH1" in str(type(histo_orig)):
        histo = ROOT.TH1D(name, name, 3, 0, 1)
        nbins_x, bin_edges_x = get_bin_info(histo_orig)
        histo_axis_nbins = nbins_x

        histo.SetBins(nbins_x, bin_edges_x.data())
        histo.GetXaxis().SetTitle(get_axis_title(histo_orig))

    elif "TH2" in str(type(histo_orig)):
        histo = ROOT.TH2D(name, name, 3, 0, 1, 3, 0, 1)
        nbins_x, bin_edges_x = get_bin_info(histo_orig, "x")
        nbins_y, bin_edges_y = get_bin_info(histo_orig, "y")
        histo_axis_nbins = (nbins_x, nbins_y)

        histo.SetBins(nbins_x, bin_edges_x.data(), nbins_y, bin_edges_y.data())
        histo.GetXaxis().SetTitle(get_axis_title(histo_orig, "x"))
        histo.GetYaxis().SetTitle(get_axis_title(histo_orig, "y"))

    elif "TH3" in str(type(histo_orig)):
        histo = ROOT.TH3D(name, name, 3, 0, 1, 3, 0, 1, 3, 0, 1)
        nbins_x, bin_edges_x = get_bin_info(histo_orig, "x")
        nbins_y, bin_edges_y = get_bin_info(histo_orig, "y")
        nbins_z, bin_edges_z = get_bin_info(histo_orig, "z")
        histo_axis_nbins = (nbins_x, nbins_y, nbins_z)

        histo.SetBins(
            nbins_x,
            bin_edges_x.data(),
            nbins_y,
            bin_edges_y.data(),
            nbins_z,
            bin_edges_z.data(),
        )
        histo.GetXaxis().SetTitle(get_axis_title(histo_orig, "x"))
        histo.GetYaxis().SetTitle(get_axis_title(histo_orig, "y"))
        histo.GetZaxis().SetTitle(get_axis_title(histo_orig, "z"))

    else:
        raise Exception("histogram not supported by ROOT")

    return histo, histo_axis_nbins


def recenter_dist(mean, std):
    half = 0.5 * (
        erf((1 - mean) / (std * np.sqrt(2))) + erf((0 - mean) / (std * np.sqrt(2)))
    )
    shifted = erfinv(half) * std * np.sqrt(2) + mean

    if abs(std) > BAD_ERROR_THRESH:
        print("    WARNING: Very large std!")
    if mean < 0 or mean > 1:
        print("    INFO: Raw mean not in [0, 1]!")
    if shifted < 0:
        print("    URGENT: Shifted mean < 0!")
    if max(abs(shifted / mean), abs(mean / shifted)) > 5:
        print("    WARNING: Raw and shifted means are significantly different!")

    print(f"    Raw mean ± std: {mean:.7f} ± {std:.7f}. Shifted mean: {shifted:.7f}")
    return shifted


def rebuild_root_histo(name, histo_orig, recenter=True):
    histo, histo_axis_nbins = prep_root_histo(name, histo_orig)

    indices_ranges = [list(range(1, n + 1)) for n in histo_axis_nbins]
    for idx in itertools.product(*indices_ranges):
        print(f"  Working on index: {idx}")
        value = histo_orig.GetBinContent(*idx)
        error = histo_orig.GetBinError(*idx)

        if recenter and not np.isnan(value) and not np.isnan(error):
            if error > 0.0:
                value = recenter_dist(value, error)
            else:
                print(f"    Zero or negative efficiency uncertainty found ({error}). Manually set efficiency to 0.0.")
                value = 0.0
        else:
            if np.isnan(value):
                value = 0.0
            if np.isnan(error):
                error = 0.0
            print(
                "    nan or 0.0 efficiency encountered. Manually set efficiency to 0.0."
            )

        histo.SetBinContent(histo.GetBin(*idx), value)
        histo.SetBinError(histo.GetBin(*idx), error)  # use the raw error

    return histo


def divide_histo(name, histo_nom, histo_denom):
    histo, histo_axis_nbins = prep_root_histo(name, histo_nom)

    indices_ranges = [list(range(1, n + 1)) for n in histo_axis_nbins]
    for idx in itertools.product(*indices_ranges):
        print(f"  Working on index: {idx}")
        nom = histo_nom.GetBinContent(*idx)
        nom_err = histo_nom.GetBinError(*idx)

        denom = histo_denom.GetBinContent(*idx)
        denom_err = histo_denom.GetBinError(*idx)

        nom_ok = not (np.isnan(nom) or np.isnan(nom_err))
        denom_ok = not (np.isnan(denom) or np.isnan(denom_err))

        if not (denom_ok and nom_ok):
            histo.SetBinContent(histo.GetBin(*idx), 0.0)
            histo.SetBinError(histo.GetBin(*idx), 0.0)
        else:
            if denom > 0.0 and nom_err > 0.0 and denom_err > 0.0:
                nom = recenter_dist(nom, nom_err)
                denom = recenter_dist(denom, denom_err)
                value = nom / denom
                error = value * np.sqrt( (nom_err / nom)**2 + (denom_err / denom)**2 )
            else:
                value = 0.0
                error = 0.0
            histo.SetBinContent(histo.GetBin(*idx), value)
            histo.SetBinError(histo.GetBin(*idx), error)

    return histo


def multiply_histo_in_place(histo1, histo2):
    print(f"  Multiplying histo {histo1.GetName()} and {histo2.GetName()}")
    _, histo_axis_nbins = prep_root_histo("tmp", histo1)

    indices_ranges = [list(range(1, n + 1)) for n in histo_axis_nbins]
    for idx in itertools.product(*indices_ranges):
        print(f"  Working on index: {idx}")
        # This only applies for e efficiency, where the PIDCalib efficiency
        # which lacks mu_BDT is multiplied by the conditional e mu_BDT efficiency
        # Due to low statistics in the MC sample, it is only computed as function of P,
        # and saved in a 1D histogram
        idx_p = idx[0]
        fac1 = histo1.GetBinContent(*idx)
        fac2 = histo2.GetBinContent(idx_p)
        print(f"    original: {fac1}")
        print(f"    new fac:  {fac2}")
        if np.isnan(fac1):
            fac1 = 0
        if np.isnan(fac2):
            fac2 = 0

        err1 = histo1.GetBinContent(*idx)
        err2 = histo2.GetBinContent(idx_p)
        if np.isnan(err1):
            err1 = 0
        if np.isnan(err2):
            err2 = 0

        if fac1 == 0 or fac2 == 0:
            err = 0
        else:
            err = fac1 * fac2 * np.sqrt((err1 / fac1) ** 2 + (err2 / fac2) ** 2)

        histo1.SetBinContent(histo1.GetBin(*idx), fac1 * fac2)
        histo1.SetBinError(histo1.GetBin(*idx), err)


########################
# Histo helpers: fixup #
########################


def is_valid_neighbor(idx_list1, idx_list2):
    diff = np.abs(np.subtract(idx_list1, idx_list2))
    if diff[diff == 1].size == 1:
        return True
    return False


def get_nearby_idx(idx, nbins):
    result = []

    for i, bin_limit in zip(idx, nbins):
        valid_idx = [i]
        if i > 1:
            valid_idx.append(i - 1)
        if i < bin_limit:
            valid_idx.append(i + 1)
        result.append(valid_idx)

    # this is a 3x3 cube (in 3D case, but it works for N-D)
    neighbors_raw = itertools.product(*result)
    # remove edges and center of the cube (center of the cube is the input idx
    # itself) (so remove 8x2+4+1 = 21 points)
    neighbors = [i for i in neighbors_raw if is_valid_neighbor(i, idx)]
    return neighbors


def compute_weighted_average(means, stds):
    wt = [1 / x**2 for x in stds]
    weighted_mean = np.sum(np.multiply(wt, means)) / np.sum(wt)
    error = np.sqrt(1 / np.sum(wt))
    return weighted_mean, error


def fix_bad_bins_in_histo(histo, bad_err_thresh=0.2):
    _, histo_axis_nbins = prep_root_histo("_tmp", histo)

    indices_ranges = [list(range(1, n + 1)) for n in histo_axis_nbins]
    for idx in itertools.product(*indices_ranges):
        mean = histo.GetBinContent(*idx)
        error = histo.GetBinError(*idx)
        if abs(error) > BAD_ERROR_THRESH:
            print(
                f"  FIX: bin {idx} = {mean:.7f} ± {error:.7f} has a large error, replace w/ weighted average of nearby bins"
            )

            means = []
            stds = []
            nearby_idx = get_nearby_idx(idx, histo_axis_nbins)
            for ni in nearby_idx:
                print(f"    Looking at nearby bin with index {ni}")
                bin_std = abs(histo.GetBinError(*ni))
                bin_mean = histo.GetBinContent(*ni)

                if bin_std > BAD_ERROR_THRESH:
                    print(f"    std = {bin_std:.7f} > {BAD_ERROR_THRESH}, skipping...")
                    continue
                if bin_mean < 0:
                    print(f"    mean = {bin_mean:.7f} < 0, skipping...")
                    continue
                if bin_std == 0.0:
                    print(
                        "    std = 0.0, perhaps the sample does not cover this bin, skipping..."
                    )
                    continue

                means.append(bin_mean)
                stds.append(bin_std)
                print(f"    Use {bin_mean:.7f} ± {bin_std:.7f}")

            new_mean, new_error = compute_weighted_average(means, stds)
            print(f"    Use weighted average: {new_mean:.7f} ± {new_error:.7f}")

            histo.SetBinContent(histo.GetBin(*idx), new_mean)
            histo.SetBinError(histo.GetBin(*idx), new_error)


###########
# Helpers #
###########


def abs_dir(path):
    return str(Path(path).parent.absolute())


def merge_true_to_tag(output_ntp, path_prefix, paths, config, year):
    # Build a dictionary of 'output histoname' -> related histos
    histo_map = dict()
    for p in paths:
        for h in glob(f"{path_prefix}/{p}/*.root"):
            histo_name = basename(h).split(".")[0]
            main_name = histo_name.split("_")[0] + "_" + year[2:]
            # Create a separate nom/denom pair for vmu effs
            if "vmu" in histo_name:
                main_name += "_vmu"
            if main_name not in histo_map:
                histo_map[main_name] = [h]
            else:
                histo_map[main_name].append(h)
            # denom is the same for vmu and iso, so in the denom iteration
            # add it to both histo lists
            if "denom" in histo_name:
                main_name += "_vmu"
                if main_name not in histo_map:
                    histo_map[main_name] = [h]
                else:
                    histo_map[main_name].append(h)

    for main_name, inputs in histo_map.items():
        # copy histo (mostly) verbatim when there's only 1 histo
        if len(inputs) == 1:
            print(f"Copy {main_name} and shift means...")
            input_ntp = ROOT.TFile(inputs[0])
            histo_orig = input_ntp.Get("eff")
            histo_out = rebuild_root_histo(main_name, histo_orig)
            fix_bad_bins_in_histo(histo_out)

            output_ntp.cd()
            histo_out.Write()

        else:
            print(f"Combining {main_name}...")
            path_nom = [p for p in inputs if "_nom" in p][0]
            path_denom = [p for p in inputs if "_denom" in p][0]

            ntp_nom = ROOT.TFile(path_nom)
            ntp_denom = ROOT.TFile(path_denom)

            histo_nom = ntp_nom.Get("eff")
            histo_denom = ntp_denom.Get("eff")
            histo_ratio = divide_histo(main_name, histo_nom, histo_denom)
            fix_bad_bins_in_histo(histo_ratio)

            for p in [i for i in inputs if "_nom" not in i and "_denom" not in i]:
                ntp_addon = ROOT.TFile(p)
                histo_addon = ntp_addon.Get("eff")
                multiply_histo_in_place(histo_ratio, histo_addon)

            output_ntp.cd()
            histo_ratio.Write()


def merge_extra(output_ntp, path_prefix, spec, config, year):
    for path in spec:
        input_ntp = ROOT.TFile(f"{path_prefix}/{path}")
        for src, tgt in spec[path].items():
            print(f"Copy {tgt} from {src} and shift means...")
            histo_src = input_ntp.Get(src)
            histo_tgt = rebuild_root_histo(tgt, histo_src)
            fix_bad_bins_in_histo(histo_tgt)

            output_ntp.cd()
            histo_tgt.Write()


########
# Main #
########

KNOWN_MERGERS = {
    "true_to_tag": merge_true_to_tag,
    "extra": merge_extra,
}


if __name__ == "__main__":
    args = parse_input()
    path_prefix = abs_dir(args.config)
    makedirs(args.output, exist_ok=True)

    with open(args.config, "r") as f:
        config = safe_load(f)

    output_ntp = ROOT.TFile(f"{args.output}/{HISTO_NAME}", "RECREATE")
    for year in config["input_histos"]:
        for mode, params in config["input_histos"][year].items():
            if mode not in KNOWN_MERGERS:
                print(f'WARNING: Unknown mode: "{mode}". Skipping...')
                continue

            print(f"Merging {mode} using {params}...")
            KNOWN_MERGERS[mode](output_ntp, path_prefix, params, config, str(year))

    output_ntp.Close()
