#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Fri Sep 23, 2022 at 04:33 AM -0400
#
# Description: build decay-in-flight histos

import numpy as np
import uproot
import mplhep

from argparse import ArgumentParser

from pyTuplingUtils.boolean.eval import BooleanEvaluator
from pyTuplingUtils.utils import gen_histo
from pyTuplingUtils.plot import plot_histo, ax_add_args_histo


################
# Configurable #
################

OUTPUT_NTP_NAME = "dif.root"
TREE_NAME = "mytuple/DecayTree"

ADD_LABELS = ["x", "y", "z"]
TRUE_P_BRS = ["TRUEP_X", "TRUEP_Y", "TRUEP_Z"]
RECO_P_BRS = ["PX", "PY", "PZ"]

####
PARTICLES = {
    "k_smr": "K",
    "pi_smr": "pi",
    "k_smr_ubdt_veto": "K",
    "pi_smr_ubdt_veto": "pi",
}

CUTS = {
    "k_smr": "K_TRACK_GhostProb < 0.5 & K_isMuon & K_PIDmu > 2.0 & BDTmuCut > 0.25 & abs(K_TRUEID) == 321",
    "pi_smr": "pi_TRACK_GhostProb < 0.5 & pi_isMuon & pi_PIDmu > 2.0 & BDTmuCut > 0.25 & abs(pi_TRUEID) == 211",
    "k_smr_ubdt_veto": "K_TRACK_GhostProb < 0.5 & K_isMuon & K_PIDmu > 2.0 & BDTmuCut < 0.25 & abs(K_TRUEID) == 321",
    "pi_smr_ubdt_veto": "pi_TRACK_GhostProb < 0.5 & pi_isMuon & pi_PIDmu > 2.0 & BDTmuCut < 0.25 & abs(pi_TRUEID) == 211",
}
####
PLOT_NBINS = 100

PLOT_PARTICLE_ALIASES = {
    "k_smr": r"$K$ (passing UBDT cut)",
    "pi_smr": r"$\pi$ (passing UBDT cut)",
    "k_smr_ubdt_veto": r"$K$ (failing UBDT cut)",
    "pi_smr_ubdt_veto": r"$\pi$ (failing UBDT cut)",
}

PLOT_VAR_ALIASES = {
    "PX": r"$p_x^{reco} / p_x^{true}$",
    "PY": r"$p_y^{reco} / p_y^{true}$",
    "PZ": r"$p_z^{reco} / p_z^{true}$",
}


#######################
# Command line parser #
#######################


def parse_input():
    parser = ArgumentParser(description="build decay-in-flight histos.")

    parser.add_argument("ntps", nargs="+", help="specify input ntuples.")
    parser.add_argument("-o", "--output", required=True, help="specify output dir.")

    parser.add_argument(
        "--plot",
        action="store_true",
        help="plot momentum ratios, for finding binnings.",
    )

    return parser.parse_args()


###########
# Helpers #
###########


def plot(br, title, filename, nbins=40):
    histo, bins = gen_histo(br, bins=nbins)
    plot_histo(
        bins,
        histo,
        ax_add_args_histo(label="ratio", color="cornflowerblue"),
        output=filename,
        title=title,
        show_legend=False,
    )


########
# Main #
########

if __name__ == "__main__":
    mplhep.style.use("LHCb2")
    args = parse_input()
    output_ntp = uproot.recreate(f"{args.output}/{OUTPUT_NTP_NAME}")

    ptcls = list(PARTICLES.keys())
    for (idx, n) in enumerate(args.ntps):
        evaluator = BooleanEvaluator(n, TREE_NAME)
        ptcl = ptcls[idx]
        prefix = PARTICLES[ptcl]
        print(f"Working on {ptcl}, with a branch prefix {prefix}")
        print(f"  Input ntuple: {n}")
        print(f"  Output tree name: {ptcl}")

        true_brs = [evaluator.eval(f"{prefix}_{b}") for b in TRUE_P_BRS]
        reco_brs = [evaluator.eval(f"{prefix}_{b}") for b in RECO_P_BRS]

        print(f"  Global cuts: {CUTS[ptcl]}")
        global_cut = evaluator.eval(CUTS[ptcl])
        output_brs = []
        output_tree = dict()

        for i in range(len(RECO_P_BRS)):
            # reco / true !!!
            ratio = np.true_divide(
                reco_brs[i],
                true_brs[i],
                out=np.zeros_like(reco_brs[i]),
                where=true_brs[i] != 0,
            )
            output_brs.append(ratio)

            if args.plot:
                title = (
                    f"{PLOT_PARTICLE_ALIASES[ptcl]}: {PLOT_VAR_ALIASES[RECO_P_BRS[i]]}"
                )
                filename = f"{args.output}/{ptcl}_{RECO_P_BRS[i]}_{TRUE_P_BRS[i]}.pdf"
                plot(ratio[global_cut], title, filename, PLOT_NBINS)

        for idx, br in enumerate(output_brs):
            br_name = f"{ptcl}_{ADD_LABELS[idx]}"
            output_tree[br_name] = br[global_cut]

        # now write the tree
        output_ntp[ptcl] = output_tree
