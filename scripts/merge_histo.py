#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Tue Apr 12, 2022 at 01:22 PM -0400
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
    return erfinv(half)*std*np.sqrt(2) + mean


def rebuild_root_histo(name, histo_orig, recenter=True):
    histo, histo_axis_nbins = prep_root_histo(name, histo_orig)

    indices_ranges = [list(range(1, n+1)) for n in histo_axis_nbins]
    for idx in itertools.product(*indices_ranges):
        value = histo_orig.GetBinContent(*idx)
        value = 0.0 if np.isnan(value) else value
        error = histo_orig.GetBinError(*idx)
        error = 0.0 if np.isnan(error) else error

        if recenter:
            value = recenter_dist(value, error)

        histo.SetBinContent(histo.GetBin(*idx), value)

    return histo


def divide_histo(name, histo_nom, histo_denom):
    histo, histo_axis_nbins = prep_root_histo(name, histo_nom)

    indices_ranges = [list(range(1, n+1)) for n in histo_axis_nbins]
    for idx in itertools.product(*indices_ranges):
        print('=====')
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

        else:
            print(nom, nom_err)
            print(denom, denom_err)
            nom = recenter_dist(nom, nom_err)
            denom = recenter_dist(denom, denom_err)

            print(nom, nom_err)
            print(denom, denom_err)
            value = nom / denom
            print(value)
            histo.SetBinContent(histo.GetBin(*idx), value)

    return histo


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
            input_ntp = ROOT.TFile(f'{path_prefix}/{path}/{histo_name}.root')
            histo_orig = input_ntp.Get('eff')
            histo_out = rebuild_root_histo(histo_name, histo_orig)

            output_ntp.cd()
            histo_out.Write()

    # Handle additional ntuples
    for p_true in ptcl_true:
        for p_addon in config['pidcalib_config']['tags_addon']:
            histo_name = f'{p_true}TrueTo{p_addon.capitalize()}Tag'

            ntp_nom = ROOT.TFile(f'{path_prefix}/{path}/{histo_name}_nom.root')
            ntp_denom = ROOT.TFile(
                f'{path_prefix}/{path}/{histo_name}_denom.root')

            histo_nom = ntp_nom.Get('eff')
            histo_denom = ntp_denom.Get('eff')
            histo_ratio = divide_histo(histo_name, histo_nom, histo_denom)

            output_ntp.cd()
            histo_ratio.Write()


def merge_extra(output_ntp, path_prefix, spec, config):
    for path in spec:
        input_ntp = ROOT.TFile(f'{path_prefix}/{path}')
        for src, tgt in spec[path].items():
            histo_src = input_ntp.Get(src)
            histo_tgt = rebuild_root_histo(tgt, histo_src)

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
