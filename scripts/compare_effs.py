#!/usr/bin/env python3
#
# Author: Lucas Meyer Garcia
# Compare PID efficiencies.
# Largely based on run_fit.py (see rdx-run2-analysis)

import sys
import datetime
from os import listdir
from os.path import exists
from subprocess import Popen, PIPE, STDOUT
from argparse import ArgumentParser

#################################
# Command line arguments parser #
#################################


def parseInput():
    parser = ArgumentParser(
        description='Run fit with configuration from yaml file.')

    parser.add_argument(
        'pathRef',
        help='Specify path to the directory containing the new efficiencies.')
    parser.add_argument(
        'pathNew',
        help='Specify path to the directory containing the old efficiencies.')

    return parser.parse_args()


#######################
# Auxiliary functions #
#######################


def cTerm(msg, color):
    num = 30
    if color == 'red': num = 91
    if color == 'green': num = 92
    if color == 'yellow': num = 93
    if color == 'blue': num = 94
    if color == 'magenta': num = 95
    if color == 'cyan': num = 96
    return f'\033[{num};1m{msg}\033[0m'


## Run shell command and print output to terminal and logfile
def runCmdLog(cmd, outFile, mode='w'):
    print('\n' + cTerm(' '.join(cmd) + ' | tee ' + outFile, 'green') + '\n')
    f = open(outFile, mode)
    with Popen(cmd,
               stdout=PIPE,
               stderr=STDOUT,
               bufsize=1,
               universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
            f.write(line)
    return p.returncode


## Run shell command and print output only to terminal
def runCmd(cmd):
    print('\n' + cTerm(' '.join(cmd), 'green') + '\n')
    with Popen(cmd,
               stdout=PIPE,
               stderr=STDOUT,
               bufsize=1,
               universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    return p.returncode


########
# Main #
########

if __name__ == '__main__':
    args = parseInput()

    # Executables
    exe = './bin/compareEffs'

    # Configuration variables
    timeStamp = datetime.datetime.now().strftime('%y_%m_%d_%H_%M')
    pathRef = args.pathRef
    pathNew = args.pathNew
    outDir = f'gen/compare_effs_{timeStamp}'
    cmdBase = [exe]

    # Get efficiencies to be compared
    # Only run for file that exist in both ref and new folders
    files = [
        fName for fName in listdir(pathNew)
        if ('.root' in fName and exists(f'{pathRef}/{fName}'))
    ]

    # Compiling
    exitCode = runCmd(['make'])
    if exitCode != 0:
        print(cTerm('\nCompilation failed, exiting\n', 'red'))
        sys.exit(1)

    for fName in files:
        fileRef = f'{pathRef}/{fName}'
        fileNew = f'{pathNew}/{fName}'
        outDirEff = outDir + '/' + fName[0:-5]  # Remove .root extension

        # Prepare output directory
        runCmd(['mkdir', '-p', outDirEff])

        # Produce comparison
        cmd = cmdBase + ['-o', outDirEff, '-r', fileRef, '-n', fileNew]
        runCmd(cmd)
