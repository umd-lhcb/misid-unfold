#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Fri Sep 23, 2022 at 04:33 AM -0400
#
# Description: build decay-in-flight histos

import numpy as np
import uproot
import mplhep
import math

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

####
def get_cuts(p):
    # Included HLT2, Stripping and Offline cuts, except:
    # - uBDT cut, which is added explicitely below.
    # - IPCHI2 cut, which shouldn't affect the smearing AND rejects most
    #   signal events on our dominantly prompt charm MC.
    cuts =    f" ETA( {p}_P, {p}_PZ ) > 1.7"
    cuts += f" & ETA( {p}_P, {p}_PZ ) < 5"
    cuts += f" & {p}_TRACK_GhostProb < 0.5"
    cuts += f" & {p}_isMuon"
    cuts += f" & {p}_PIDmu > 2.0"
    cuts += f" & {p}_PIDe < 1"
    cuts += f" & {p}_P > 3 * GeV"
    cuts += f" & {p}_P < 100 * GeV"
    cuts +=  " & nPVs > 0"
    cuts +=  " & nSPDHits < 450"
    cuts +=  " & LOG10pp(K_PX, K_PY, K_PZ, pi_PX, pi_PY, pi_PZ) > -6.5"
    # Truth-matching cuts
    cuts +=  " & abs(D_TRUEID) == 421"
    cuts +=  " & pi_TRUEORIGINVERTEX_TYPE == 2"
    cuts +=  " & K_TRUEORIGINVERTEX_TYPE == 2"
    #
    # cuts += f" & abs({p}_TRUEP_X) > 50"
    # cuts += f" & abs({p}_TRUEP_Y) > 50"
    return cuts
CUTS = {
    "k_smr":            get_cuts("K")  + " & abs(pi_TRUEID) == 211 & BDTmuCut > 0.25",
    "pi_smr":           get_cuts("pi") + " & abs(K_TRUEID)  == 321 & BDTmuCut > 0.25",
    "k_smr_ubdt_veto":  get_cuts("K")  + " & abs(pi_TRUEID) == 211 & BDTmuCut < 0.25",
    "pi_smr_ubdt_veto": get_cuts("pi") + " & abs(K_TRUEID)  == 321 & BDTmuCut < 0.25",
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
        print(f"\nWorking on {ptcl}, with a branch prefix {prefix}")
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

        true_px = evaluator.eval(f"{prefix}_TRUEP_X")
        true_py = evaluator.eval(f"{prefix}_TRUEP_Y")
        true_pz = evaluator.eval(f"{prefix}_TRUEP_Z")
        true_pt = evaluator.eval(f"{prefix}_TRUEPT")
        true_p = np.sqrt( np.add( np.power(true_pt, 2), np.power(true_pz, 2) ) )
        reco_px = evaluator.eval(f"{prefix}_PX")
        reco_py = evaluator.eval(f"{prefix}_PY")
        reco_pz = evaluator.eval(f"{prefix}_PZ")
        reco_pt = evaluator.eval(f"{prefix}_PT")
        reco_p = evaluator.eval(f"{prefix}_P")

        # Get delta_theta
        true_theta = np.arccos(np.divide(true_pz,true_p))
        reco_theta = np.arccos(np.divide(reco_pz,reco_p))
        delta_theta = np.subtract(reco_theta, true_theta)

        output_brs.append(delta_theta)
        output_tree[f"{ptcl}_dTheta"] = delta_theta[global_cut]

        # Get delta_phi
        true_abs_phi = np.arccos(np.divide(true_px,true_pt))
        reco_abs_phi = np.arccos(np.divide(reco_px,reco_pt))
        # For y < 0, phi is actually -1 * acos(px/pt)
        true_py_sign = np.sign(true_py)
        true_phi = np.multiply(true_py_sign, true_abs_phi)
        reco_py_sign = np.sign(reco_py)
        reco_phi = np.multiply(reco_py_sign, reco_abs_phi)
        delta_phi = np.subtract(reco_phi, true_phi)
        # Set delta_phi -> 2pi - delta_phi for delta_phi > pi
        twopis = np.full_like(delta_phi, 2. * math.pi)
        delta_phi = np.subtract(twopis, delta_phi, out=delta_phi, where=delta_phi>math.pi)
        # Set delta_phi -> delta_phi + 2pi for delta_phi <= -pi
        delta_phi = np.add(delta_phi, twopis, out=delta_phi, where=delta_phi<=-math.pi)

        output_brs.append(delta_phi)
        output_tree[f"{ptcl}_dPhi"] = delta_phi[global_cut]

        # Get P ratio
        ratio_p = np.divide(reco_p, true_p, out=np.zeros_like(reco_p), where=true_p>0)

        output_brs.append(ratio_p)
        output_tree[f"{ptcl}_rP"] = ratio_p[global_cut]

        # now write the tree
        output_ntp[ptcl] = output_tree
