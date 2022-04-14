#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Wed Apr 13, 2022 at 11:16 PM -0400
#
# Description: histogram merger (M)

import itertools
import numpy as np

from argparse import ArgumentParser
from pathlib import Path
from os import makedirs
from yaml import safe_load
from scipy.special import erf, erfinv

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True  # Don't hijack argparse!
ROOT.PyConfig.DisableRootLogon = True  # Don't read .rootlogon.py


################
# Configurable #
################

HISTO_NAME = 'merged.root'
BAD_ERROR_THRESH = 0.2


#######################
# Command line parser #
#######################

def parse_input():
    parser = ArgumentParser(description='histogram merger (M).')

    parser.add_argument('-c', '--config', required=True,
                        help='specify YAML config.')

    parser.add_argument('-o', '--output', required=True,
                        help='specify output dir.')

    parser.add_argument('-y', '--year', default='2016',
                        help='specify year.')

    return parser.parse_args()


#################
# Histo helpers #
#################

def get_bin_info(histo, axis='x'):
    nbins = getattr(histo, f'GetNbins{axis.upper()}')()
    axis = getattr(histo, f'Get{axis.upper()}axis')()
    bin_edges = ROOT.std.vector('double')()

    for i in range(1, nbins+2):
        bin_edges.push_back(getattr(axis, 'GetBinLowEdge')(i))

    return nbins, bin_edges


def get_axis_title(histo, axis='x'):
    axis = getattr(histo, f'Get{axis.upper()}axis')()
    return getattr(axis, 'GetTitle')()


def prep_root_histo(name, histo_orig):
    histo = None
    histo_axis_nbins = None

    if 'TH1' in str(type(histo_orig)):
        histo = ROOT.TH1D(name, name, 3, 0, 1)
        nbins_x, bin_edges_x = get_bin_info(histo_orig)
        histo_axis_nbins = (nbins_x)

        histo.SetBins(nbins_x, bin_edges_x.data())
        histo.GetXaxis().SetTitle(get_axis_title(histo_orig))

    elif 'TH2' in str(type(histo_orig)):
        histo = ROOT.TH2D(name, name, 3, 0, 1, 3, 0, 1)
        nbins_x, bin_edges_x = get_bin_info(histo_orig, 'x')
        nbins_y, bin_edges_y = get_bin_info(histo_orig, 'y')
        histo_axis_nbins = (nbins_x, nbins_y)

        histo.SetBins(nbins_x, bin_edges_x.data(),
                      nbins_y, bin_edges_y.data())
        histo.GetXaxis().SetTitle(get_axis_title(histo_orig, 'x'))
        histo.GetYaxis().SetTitle(get_axis_title(histo_orig, 'y'))

    elif 'TH3' in str(type(histo_orig)):
        histo = ROOT.TH3D(name, name, 3, 0, 1, 3, 0, 1, 3, 0, 1)
        nbins_x, bin_edges_x = get_bin_info(histo_orig, 'x')
        nbins_y, bin_edges_y = get_bin_info(histo_orig, 'y')
        nbins_z, bin_edges_z = get_bin_info(histo_orig, 'z')
        histo_axis_nbins = (nbins_x, nbins_y, nbins_z)

        histo.SetBins(nbins_x, bin_edges_x.data(),
                      nbins_y, bin_edges_y.data(),
                      nbins_z, bin_edges_z.data())
        histo.GetXaxis().SetTitle(get_axis_title(histo_orig, 'x'))
        histo.GetYaxis().SetTitle(get_axis_title(histo_orig, 'y'))
        histo.GetZaxis().SetTitle(get_axis_title(histo_orig, 'z'))

    else:
        raise Exception("histogram not supported by ROOT")

    return histo, histo_axis_nbins


def recenter_dist(mean, std):
    half = 0.5*(erf((1-mean)/(std*np.sqrt(2))) + erf((0-mean)/(std*np.sqrt(2))))
    shifted = erfinv(half)*std*np.sqrt(2) + mean

    if abs(std) > BAD_ERROR_THRESH:
        print('    WARNING: Very large std!')
    if mean < 0 or mean > 1:
        print('    INFO: Raw mean not in [0, 1]!')
    if shifted < 0:
        print('    URGENT: Shifted mean < 0!')
    if max(abs(shifted / mean), abs(mean / shifted)) > 5:
        print('    WARNING: Raw and shifted means are significantly different!')

    print(f'    Raw mean ± std: {mean:.7f} ± {std:.7f}. Shifted mean: {shifted:.7f}')
    return shifted


def rebuild_root_histo(name, histo_orig, recenter=True):
    histo, histo_axis_nbins = prep_root_histo(name, histo_orig)

    indices_ranges = [list(range(1, n+1)) for n in histo_axis_nbins]
    for idx in itertools.product(*indices_ranges):
        print(f'  Working on index: {idx}')
        value = histo_orig.GetBinContent(*idx)
        value = 0.0 if np.isnan(value) else value
        error = histo_orig.GetBinError(*idx)
        error = 0.0 if np.isnan(error) else error

        if recenter and value*error != 0:
            value = recenter_dist(value, error)
        else:
            print(f'    nan or 0.0 efficiency encountered. Manually set efficiency to 0.0.')

        histo.SetBinContent(histo.GetBin(*idx), value)
        histo.SetBinError(histo.GetBin(*idx), error)  # use the raw error

    return histo


def divide_histo(name, histo_nom, histo_denom):
    histo, histo_axis_nbins = prep_root_histo(name, histo_nom)

    indices_ranges = [list(range(1, n+1)) for n in histo_axis_nbins]
    for idx in itertools.product(*indices_ranges):
        print(f'  Working on index: {idx}')
        nom = histo_nom.GetBinContent(*idx)
        nom = 0.0 if np.isnan(nom) else nom
        nom_err = histo_nom.GetBinError(*idx)
        nom_err = 0.0 if np.isnan(nom_err) else nom_err

        denom = histo_denom.GetBinContent(*idx)
        denom = 0.0 if np.isnan(denom) else denom
        denom_err = histo_denom.GetBinError(*idx)
        denom_err = 0.0 if np.isnan(denom_err) else denom_err

        if denom == 0.0 or nom == 0.0:
            histo.SetBinContent(histo.GetBin(*idx), 0.0)
            histo.SetBinError(histo.GetBin(*idx), 0.0)

        else:
            nom = recenter_dist(nom, nom_err)
            denom = recenter_dist(denom, denom_err)
            value = nom / denom
            histo.SetBinContent(histo.GetBin(*idx), value)
            histo.SetBinError(histo.GetBin(*idx), nom_err)

    return histo


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
            valid_idx.append(i-1)
        if i < bin_limit:
            valid_idx.append(i+1)
        result.append(valid_idx)

    # this is a 3x3 cube (in 3D case, but it works for N-D)
    neighbors_raw = itertools.product(*result)
    # remove edges and center of the cube (center of the cube is the input idx
    # itself) (so remove 8x2+4+1 = 21 points)
    neighbors = [i for i in neighbors_raw if is_valid_neighbor(i, idx)]
    return  neighbors


def compute_weighted_average(means, stds):
    wt = [1/x**2 for x in stds]
    weighted_mean = np.sum(np.multiply(wt, means)) / np.sum(wt)
    error = np.sqrt(1 / np.sum(wt))
    return weighted_mean, error


def fix_bad_bins_in_histo(histo, bad_err_thresh=0.2):
    _, histo_axis_nbins = prep_root_histo('_tmp', histo)

    indices_ranges = [list(range(1, n+1)) for n in histo_axis_nbins]
    for idx in itertools.product(*indices_ranges):
        mean = histo.GetBinContent(*idx)
        error = histo.GetBinError(*idx)
        if abs(error) > BAD_ERROR_THRESH:
            print(f'  FIX: bin {idx} = {mean:.7f} ± {error:.7f} has a large error, replace w/ weighted average of nearby bins')

            means = []
            stds = []
            nearby_idx = get_nearby_idx(idx, histo_axis_nbins)
            for ni in nearby_idx:
                print(f'    Looking at nearby bin with index {ni}')
                bin_std = abs(histo.GetBinError(*ni))
                bin_mean = histo.GetBinContent(*ni)

                if bin_std > BAD_ERROR_THRESH:
                    print(f'    std = {bin_std:.7f} > {BAD_ERROR_THRESH}, skipping...')
                    continue
                if bin_mean < 0:
                    print(f'    mean = {bin_mean:.7f} < 0, skipping...')
                    continue
                if bin_std == 0.0:
                    print('    std = 0.0, perhaps the sample does not cover this bin, skipping...')
                    continue

                means.append(bin_mean)
                stds.append(bin_std)
                print(f'    Use {bin_mean:.7f} ± {bin_std:.7f}')

            new_mean, new_error = compute_weighted_average(means, stds)
            print(f'    Use weighted average: {new_mean:.7f} ± {new_error:.7f}')

            histo.SetBinContent(histo.GetBin(*idx), new_mean)
            histo.SetBinError(histo.GetBin(*idx), new_error)


###########
# Helpers #
###########

def abs_dir(path):
    return str(Path(path).parent.absolute())


def merge_true_to_tag(output_ntp, path_prefix, path, config):
    ptcl_true = list(config['pidcalib_config']['samples'])
    ptcl_tag = list(config['tags'])

    # Handle normal trueToTag first
    for p_true in ptcl_true:
        for p_tag in ptcl_tag:
            histo_name = f'{p_true}TrueTo{p_tag.capitalize()}Tag'
            print(f'Copy {histo_name} and shift means...')

            input_ntp = ROOT.TFile(f'{path_prefix}/{path}/{histo_name}.root')
            histo_orig = input_ntp.Get('eff')
            histo_out = rebuild_root_histo(histo_name, histo_orig)
            fix_bad_bins_in_histo(histo_out)

            output_ntp.cd()
            histo_out.Write()

    # Handle additional ntuples
    for p_true in ptcl_true:
        for p_addon in config['pidcalib_config']['tags_addon']:
            histo_name = f'{p_true}TrueTo{p_addon.capitalize()}Tag'
            print(f'Copy {histo_name} and shift means...')

            ntp_nom = ROOT.TFile(f'{path_prefix}/{path}/{histo_name}_nom.root')
            ntp_denom = ROOT.TFile(
                f'{path_prefix}/{path}/{histo_name}_denom.root')

            histo_nom = ntp_nom.Get('eff')
            histo_denom = ntp_denom.Get('eff')
            histo_ratio = divide_histo(histo_name, histo_nom, histo_denom)
            fix_bad_bins_in_histo(histo_ratio)

            output_ntp.cd()
            histo_ratio.Write()


def merge_extra(output_ntp, path_prefix, spec, config):
    for path in spec:
        input_ntp = ROOT.TFile(f'{path_prefix}/{path}')
        for src, tgt in spec[path].items():
            print(f'Copy {tgt} from {src} and shfit means...')
            histo_src = input_ntp.Get(src)
            histo_tgt = rebuild_root_histo(tgt, histo_src)
            fix_bad_bins_in_histo(histo_tgt)

            output_ntp.cd()
            histo_tgt.Write()


########
# Main #
########

KNOWN_MERGERS = {
    'true_to_tag': merge_true_to_tag,
    'extra': merge_extra,
}


if __name__ == '__main__':
    args = parse_input()
    path_prefix = abs_dir(args.config)
    makedirs(args.output, exist_ok=True)

    with open(args.config, 'r') as f:
        config = safe_load(f)

    output_ntp = ROOT.TFile(f'{args.output}/{HISTO_NAME}', 'RECREATE')
    for mode, params in config['input_histos'][int(args.year)].items():
        print(f'Merging {mode}...')
        KNOWN_MERGERS[mode](output_ntp, path_prefix, params, config)
    output_ntp.Close()
