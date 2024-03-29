#!/usr/bin/env python
# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Fri Jun 24, 2022 at 03:06 AM -0400
#
# Description: histogram plotter (for this project)

import re
import sys
import uproot
import mplhep
import numpy as np
import matplotlib.pyplot as plt

from argparse import ArgumentParser
from pyTuplingUtils.utils import gen_histo_stacked_baseline
from pyTuplingUtils.plot import plot_prepare, plot_histo, ax_add_args_histo


################
# Configurable #
################

# generated w/ https://gka.github.io/palettes of the input colors:
#   #00429d, #96ffea, #8bd189
#   #8bd189, #ff005e, #93003a
DEFAULT_COLORS = ["#00429d", "#5585b7", "#8bd189", "#d15d5f", "#93003a"]

KNOWN_PARTICLES = {
    "d0": r"$D^0$",
    "dst": r"$D^*$",
    "e": r"$e$",
    "pi": r"$\pi$",
    "p": r"$p$",
    "k": r"$K$",
    "g": "ghost",
    "mu": r"$\mu$",
    "d0_usb": r"$D^0$ (USB)",
    "dst_usb": r"$D^*$ (USB)",
    "dst_comb_d": r"$D^*$ ($D^*$ comb)",
    "dst_comb_dsb": r"$D^*$ ($D^*$ comb SB)",
}

TAG_ORDERING = ["pi", "k", "p", "e", "g"]

GET_PARTICLE = (
    lambda x: KNOWN_PARTICLES[x.lower()] if x.lower() in KNOWN_PARTICLES else x
)


#######################
# Command line helper #
#######################


def parse_input():
    parser = ArgumentParser(description="histogram plotter (for this project).")

    parser.add_argument("-i", "--input", help="specify input ntuple.")

    parser.add_argument("-o", "--output", help="specify output folder.")

    parser.add_argument(
        "-p",
        "--prefix",
        nargs="+",
        default=["D0"],
        help="specify histo prefixes to plot.",
    )

    parser.add_argument(
        "-s", "--suffix", default="True", help="specify histo suffix to plot."
    )

    parser.add_argument(
        "-v",
        "--vars",
        nargs="+",
        default=["p", "eta", "ntracks"],
        help="specify histo binning variables.",
    )

    parser.add_argument(
        "-l",
        "--labels",
        nargs="+",
        default=[r"$p$ [MeV]", r"$\eta$", r"nTracks"],
        help="specify histo binning variables in displayed format.",
    )

    parser.add_argument(
        "--show-title",
        nargs="+",
        type=int,
        default=[1, 1, 1],
        help="control display of title for each individual plot.",
    )

    parser.add_argument(
        "--show-legend",
        nargs="+",
        type=int,
        default=[1, 1, 1],
        help="control display of legend for each individual plot.",
    )

    parser.add_argument(
        "--extension", default="png", help="specify output file extension."
    )

    return parser.parse_args()


###########
# Helpers #
###########


def find_key(name, pref, suffix, sep="__"):
    for p in TAG_ORDERING:
        if pref + sep + p + suffix in name:
            return True
    return False


def name_cleanup(name):
    return name.split(";")[0]


def prefix_gen(histo_spec):
    raw_name = list(histo_spec)[0]
    ptcl, descr = raw_name.split("__")
    mod = ""
    if "true" in descr.lower():
        mod = "true"
    elif "tag" in descr.lower():
        mod = "tag"
    return f"{ptcl}-{mod}"


def label_gen(name, known=KNOWN_PARTICLES):
    name = name_cleanup(name).split("__")[1]
    regex = re.search(r"(\w+)(True|Tag)", name)
    return f"{GET_PARTICLE(regex.group(1))} {regex.group(2).lower()}"


def title_gen(name, known=KNOWN_PARTICLES):
    ptcl, descr = name.split("-")
    return f"{GET_PARTICLE(ptcl)}, {descr.lower()}"


########
# Plot #
########


def plot(
    histo_spec,
    bin_vars,
    bin_names,
    output_dir,
    ordering,
    show_title_lst,
    show_legend_lst,
    colors=DEFAULT_COLORS,
    ext="pdf",
):
    prefix = prefix_gen(histo_spec)

    for idx, v in enumerate(bin_vars):
        labels = [label_gen(i) for i in ordering]
        histos = [
            np.sum(histo_spec[i][0], axis=tuple(x for x in range(3) if x != idx))
            for i in ordering
        ]
        binspecs = [histo_spec[i][1 + idx] for i in ordering]
        baselines = gen_histo_stacked_baseline(histos)
        show_title = show_title_lst[idx]
        show_legend = show_legend_lst[idx]

        plotters = []
        for lbl, hist, bins, bot, clr in zip(
            labels, histos, binspecs, baselines, colors
        ):
            add_args = ax_add_args_histo(lbl, clr, baseline=bot)
            plotters.append(
                lambda fig, ax, b=bins, h=hist + bot, add=add_args: plot_histo(
                    b, h, add, figure=fig, axis=ax, show_legend=False
                )
            )

            title = title_gen(prefix) if show_title else "   "
            fig, ax, _ = plot_prepare(
                xlabel=bin_names[idx], title=title, show_legend=False
            )

            for p in plotters:
                p(fig, ax)

            if show_legend:
                handles, leg_lbls = ax.get_legend_handles_labels()
                ax.legend(
                    handles[::-1],
                    leg_lbls[::-1],
                    numpoints=1,
                    loc="best",
                    frameon="true",
                )

            fig.savefig(f"{output_dir}/{prefix}_{v}.{ext}")
            plt.close("all")


########
# Main #
########

if __name__ == "__main__":
    mplhep.style.use("LHCb2")
    args = parse_input()
    ntp = uproot.open(args.input)

    # group histograms
    for pref in args.prefix:
        histos = {
            name: ntp[k].to_numpy()
            for k in ntp
            if find_key(name := name_cleanup(k), pref, args.suffix)
        }

        if len(histos) == 0:
            print(f"No matching histo for prefix {pref}!")
            sys.exit(1)

        print("Histos to plot:")
        for h in histos:
            print(f"  {h}")

        ordering = [f"{pref}__{p}{args.suffix}" for p in TAG_ORDERING]
        plot(
            histos,
            args.vars,
            args.labels,
            args.output,
            ordering,
            args.show_title,
            args.show_legend,
            ext=args.extension,
        )
