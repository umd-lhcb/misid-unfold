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
        type=float,
        default=0.05,
        help="specify a threshold abs. diff. above which is considered bad.",
    )
    parser.add_argument(
        "-i", "--idx", nargs="+", default=None, help="specify idx to print."
    )
    return parser.parse_args()


########
# Main #
########

if __name__ == "__main__":
    args = parse_input()
    ntp = TFile.Open(args.ntp)
    thresh = args.threshold
    idx_to_print = (
        [tuple(int(i) for i in s.split(",")) for s in args.idx] if args.idx else []
    )

    for s in SPECIES:
        name_ref = f"{REF_PREFIX}__{s}"
        name_comp = f"{COMP_PREFIX}__{s}"

        histo_ref = ntp.Get(name_ref)
        histo_comp = ntp.Get(name_comp)

        xbins = histo_ref.GetNbinsX()
        ybins = histo_ref.GetNbinsY()
        zbins = histo_ref.GetNbinsZ()
        ndof = xbins * ybins * zbins

        forced_print = []
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

            if (i, j, k) in idx_to_print:
                forced_print.append(
                    ((i, j, k), diff, val_ref, err_ref, val_comp, err_comp)
                )

            if err < 1e-6:
                continue

            chi2 = chi2 + diff**2 / err**2
            # if diff > thresh:
            if diff > err_ref + err_comp:
                bad_bins.append(((i, j, k), diff, val_ref, err_ref, val_comp, err_comp))

        print(f"Comparing {s}...")
        print(f"  chi2 = {chi2:.3g}, ndof = {ndof}, chi2ndof = {chi2/ndof:.3g}")

        if forced_print:
            print(
                f"  These bins are requested: (bin idx, abs. diff, RDX. val, RJpsi val)"
            )
            for idx, abs_diff, val_ref, err_ref, val_comp, err_comp in forced_print:
                print(
                    f"    {idx}\t{abs_diff:.3f}\t{val_ref:.3f}±{err_ref:.3f}\t{val_comp:.3f}±{err_comp:.3f}"
                )

        print(
            f"  These are the worst {len(bad_bins)} bins: (bin idx, abs. diff, RDX. val, RJpsi val)"
        )
        for idx, abs_diff, val_ref, err_ref, val_comp, err_comp in sorted(
            bad_bins, key=lambda x: x[1], reverse=True
        ):
            print(
                f"    {idx}\t{abs_diff:.3f}\t{val_ref:.3f}±{err_ref:.3f}\t{val_comp:.3f}±{err_comp:.3f}"
            )
