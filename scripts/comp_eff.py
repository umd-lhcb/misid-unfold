#!/usr/bin/env python

from argparse import ArgumentParser
from itertools import product
from math import sqrt

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True  # Don't hijack argparse!
ROOT.PyConfig.DisableRootLogon = True  # Don't read .rootlogon.py

from ROOT import TFile


################
# Configurable #
################

SPECIES = [f"gTrueTo{s}Tag" for s in ["Pi", "K", "P", "E", "G"]]
REF_PREFIX = "ghost"
COMP_PREFIX = "ghost_Jpsi"


#######################
# Command line parser #
#######################


def parse_input():
    parser = ArgumentParser(description="tagged histogram builder (T).")
    parser.add_argument("ntp", help="specify input ntuple.")
    parser.add_argument(
        "-t",
        "--threshold",
        default=0.05,
        help="specify a threshold abs. diff. above which is considered bad.",
    )
    return parser.parse_args()


########
# Main #
########

if __name__ == "__main__":
    args = parse_input()
    ntp = TFile.Open(args.ntp)
    thresh = args.threshold

    for s in SPECIES:
        name_ref = f"{REF_PREFIX}__{s}"
        name_comp = f"{COMP_PREFIX}__{s}"

        histo_ref = ntp.Get(name_ref)
        histo_comp = ntp.Get(name_comp)

        xbins = histo_ref.GetNbinsX()
        ybins = histo_ref.GetNbinsY()
        zbins = histo_ref.GetNbinsZ()
        ndof = xbins * ybins * zbins

        bad_bins = []
        chi2 = 0

        for i, j, k in product(
            range(1, xbins + 1), range(1, ybins + 1), range(1, zbins + 1)
        ):
            idx = histo_ref.GetBin(i, j, k)

            val_ref = histo_ref.GetBinContent(idx)
            err_ref = histo_ref.GetBinError(idx)

            val_comp = histo_comp.GetBinContent(idx)
            err_comp = histo_comp.GetBinError(idx)

            diff = abs(val_ref - val_comp)
            err = sqrt(err_ref**2 + err_comp**2)

            if err < 1e-6:
                continue

            chi2 = chi2 + diff**2 / err**2
            if diff > thresh:
                bad_bins.append(((i, j, k), diff, val_ref, val_comp))

        print(f"Comparing {s}...")
        print(f"  chi2 = {chi2}, ndof = {ndof}, chi2ndof = {chi2/ndof}")
        print(
            f"  These are the worst {len(bad_bins)} bins: (bin idx, abs. diff, RDX. val, RJpsi val)"
        )
        for idx, abs_diff, val_ref, val_comp in sorted(
            bad_bins, key=lambda x: x[1], reverse=True
        ):
            print(f"    {idx}\t{abs_diff}\t{val_ref}\t{val_comp}")