#!/usr/bin/env python
###############################################################################
# (c) Copyright 2021 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
# NOTE: This is modified by me to add ability to specify output filename

"""Module to make LHCb PID efficiency histograms.

This module creates histograms that can be used to estimate the PID
efficiency of a user's sample.

Examples:
    Create a single efficiency histogram for a single PID cut::

        $ python -m src.pidcalib2.make_eff_hists --sample=Turbo18 --magnet=up \
            --particle=Pi --pid-cut="DLLK > 4" --bin-var=P --bin-var=ETA \
            --bin-var=nSPDHits --output-dir=pidcalib_output

    Create multiple histograms in one run (most of the time is spent reading
    in data, so specifying multiple cuts is much faster than running
    make_eff_hists sequentially)::

        $ python -m src.pidcalib2.make_eff_hists --sample=Turbo16 --magnet=up \
            --particle=Pi --pid-cut="DLLK > 0" --pid-cut="DLLK > 4" \
            --pid-cut="DLLK > 6" --bin-var=P --bin-var=ETA \
            --bin-var=nSPDHits --output-dir=pidcalib_output
"""


import argparse
import logging
import pathlib
import pickle
import re
import sys

import logzero
from logzero import logger as log

from pidcalib2 import argparse_actions, binning, utils

try:
    from pidcalib2.version import version  # type: ignore
except ImportError:
    version = "N/A"


def decode_arguments(args):
    """Decode CLI arguments."""
    parser = argparse.ArgumentParser(
        allow_abbrev=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-s",
        "--sample",
        help="calibration sample; see pidcalib2.make_eff_hists --list configs",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-m",
        "--magnet",
        help="magnet polarity",
        required=True,
        choices=["up", "down"],
    )
    parser.add_argument(
        "-p",
        "--particle",
        help="particle type; see pidcalib2.make_eff_hists --list configs",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--pid-cut",
        help=(
            "PID cut string, e.g., 'DLLK < 4.0' (-i can be used multiple times for "
            "multiple cuts)."
        ),
        action="append",
        dest="pid_cuts",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--cut",
        help=(
            "arbitrary cut string, e.g., 'Dst_IPCHI2 < 10.0' (-c can be used multiple "
            "times for multiple cuts)."
        ),
        action="append",
        dest="cuts",
    )
    parser.add_argument(
        "-b",
        "--bin-var",
        help=(
            "binning variable (-b can be used multiple times for a multi-dimensional "
            "binning)"
        ),
        action="append",
        dest="bin_vars",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--binning-file",
        help="file where new/alternative binnings are defines",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="pidcalib_output",
        help="directory where to save output files",
    )
    parser.add_argument(
        "-l",
        "--list",
        action=argparse_actions.ListValidAction,
        help="list all [configs, aliases]",
    )
    parser.add_argument(
        "-d",
        "--local-dataframe",
        help="(debug) read a calibration DataFrame from file",
    )
    parser.add_argument(
        "-f",
        "--file-list",
        help="(debug) read calibration file paths from a text file",
    )
    parser.add_argument(
        "-a",
        "--samples-file",
        help="(debug) read calibration sample lists from a custom file",
    )
    parser.add_argument(
        "-n",
        "--max-files",
        type=int,
        help="(debug) a max number of files to read",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="(debug) increase verbosity",
    )

    parser.add_argument(
        "-z",
        "--pkl-name",
        action="append",
        help="specify output pkl filename"
    )
    parser.add_argument("-V", "--version", action="version", version=version)
    return parser.parse_args(args)


def make_eff_hists(config: dict) -> None:
    """Create sWeighted PID calibration histograms and save them to disk.

    Calibration samples from EOS are read and relevant branches extracted to
    a DataFrame. Each PID cut is applied to the DataFrame in turn and the
    results are histogrammed (each event with its associated sWeight).
    Particle type and binning variables are used to select an appropriate
    predefined binning. The efficiency histograms are saved to a requested
    output directory.

    Args:
        config: A configuration dictionary. See decode_arguments(args) for
            details.
    """
    if config["verbose"]:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)

    pattern = re.compile(r"\s+")
    config["pid_cuts"] = [
        re.sub(pattern, "", pid_cut) for pid_cut in config["pid_cuts"]
    ]

    config["version"] = version
    log.info("Running PIDCalib2 make_eff_hists with the following config:")
    utils.log_config(config)

    if not config["cuts"]:
        log.warning(
            "No --cut specified. Cuts on PID samples should match your sample cuts."
        )

    output_dir = pathlib.Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    binning.check_and_load_binnings(
        config["particle"],
        config["bin_vars"],
        config["binning_file"] if "binning_file" in config else None,
    )

    if config["local_dataframe"]:
        hists = utils.create_histograms_from_local_dataframe(config)
    else:
        hists_list = utils.create_histograms(config)
        hists = utils.add_hists(list(hists_list.values()))

    eff_hists = utils.create_eff_histograms(hists)
    eff_hists_filtered = {k: v for k, v in eff_hists.items() if k.startswith("eff_")}

    for name, output_filename in zip(eff_hists_filtered, config["pkl_name"]):
        if name.startswith("eff_"):
            cut = name.replace("eff_", "")
            eff_hist_path = output_dir / output_filename
            with open(eff_hist_path, "wb") as f:
                pickle.dump(eff_hists[f"eff_{cut}"], f)
                pickle.dump(eff_hists[f"passing_{cut}"], f)
                pickle.dump(eff_hists["total"], f)

            log.info(f"Efficiency histogram of cut {cut} saved to '{eff_hist_path}'")


def main():
    config = vars(decode_arguments(sys.argv[1:]))
    make_eff_hists(config)


if __name__ == "__main__":
    main()
