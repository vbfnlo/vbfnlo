# This file can be used to modify the input files vbfnlo.dat in all regress/runs subdirs to either a fixed value, or to the value given in the result.out files of each subdir

import argparse
import glob
import logging
import os.path
import re

# custom modules
import parsedat

from os.path import exists


def set_config(args):
    alldirs = [os.path.dirname(entry) for name in args.testdirs for entry in glob.glob(name + '/vbfnlo.dat', recursive=True)]
    logging.info("Found dirs: " + ", ".join(alldirs))

    for checkdir in alldirs:
        inputfile = os.path.join(checkdir, 'result.out')
        outputfile = os.path.join(checkdir, 'vbfnlo.dat')

        EcmSet = '13600d0'
        MhSet = '125.0d0'
        mtopSet = '172.5d0'

        if not exists(outputfile):
            print("Missing file: ", outputfile)
            continue
        if not exists(inputfile):
            print("Missing file: ", inputfile)
            continue

        if args.stat == "low":
            lo_points  = 10
            nlo_points = 10
            lo_its  = 1
            nlo_its = 1
        elif args.stat == "mid":
            lo_points  = 20
            nlo_points = 20
            lo_its  = 6
            nlo_its = 6
        elif args.stat == "high":
            lo_points  = extract_var(" LO_POINTS",inputfile)
            nlo_points = extract_var(" NLO_POINTS",inputfile)
            lo_its  = extract_var(" LO_ITERATIONS",inputfile)
            nlo_its = extract_var(" NLO_ITERATIONS",inputfile)

        inputvars = {
            'LO_POINTS': str(lo_points) + '\t\t! number of points for LO calculation (= 2^...',
            'NLO_POINTS': str(nlo_points) + '\t\t! number of points for real-emissions calc. (= 2^...) ',
            'LO_ITERATIONS': str(lo_its) + '\t\t! number of iterations for LO calculation',
            'NLO_ITERATIONS': str(nlo_its) + '\t\t! number of iterations for real-emissions calc.',
            'ECM': EcmSet + 5*'\t' + '! collider center-of-mass energy',
            'HMASS': MhSet + 4*'\t' + '! Higgs mass',
            'TOPMASS': mtopSet + 4*'\t' + '! Top mass'
        }

        if args.stat == "high":
            logging.info(" New data written to from " + inputfile + " to "  + outputfile)
        else:
            logging.info(" New data written to "+ outputfile)

        for key, val in inputvars.items():
            parsedat.replacemachine(outputfile, key, val, keep=False)

def extract_var(string, filename):
    file = open(filename, "r")
    for line in file:
        if string in line:
            return [int(i) for i in line.split() if i.isdigit()][0]
    raise Exception("Expected variables not found in " + str(filename) + ".")


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--verbose',
                        '-v',
                        action='count',
                        help='Generate more output.',
                        default=0)
    parser.add_argument('--stat',
                        default="high",
                        help='Choose whether to run with low, mid or high statistics.'
                             ' Default: high')
    parser.add_argument('--testdirs',
                        nargs="*",
                        default=["runs/**"],
                        help='List of test directories, can be shell glob.'
                             ' Default: runs/**')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    if args.verbose == 1:
        logging.getLogger().setLevel(logging.INFO)
    elif args.verbose >= 2:
        logging.getLogger().setLevel(logging.DEBUG)
    set_config(args)
