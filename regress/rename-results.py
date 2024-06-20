# This file can be used to rename vbfnlo output files to keep the as a backup

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
        inputfile = os.path.join(checkdir, args.inname)
        outputfile = os.path.join(checkdir, args.outname)
        os.rename(inputfile, outputfile)
        logging.info(" New data written to " + outputfile)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--verbose',
                        '-v',
                        action='count',
                        help='Generate more output.',
                        default=0)
    parser.add_argument('--inname',
                        default="result.out",
                        help='Name of input file.'
                             ' Default: result.out')
    parser.add_argument('--outname',
                        default="result_backup.out",
                        help='Name of output/renamed file.'
                             ' Default: result_backup.out')
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
