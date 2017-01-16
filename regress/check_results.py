#!/usr/bin/env python3

# standard modules
import glob
import logging
import argparse
import configparser
import shutil

# custom modules
import parsedat
from test_runs import get_cross_section, print_diff, print_short_diff, get_diff

import make_runs

logging.getLogger().setLevel(logging.WARNING)


def check_results(args):
    alldirs = [entry for name in args.testdirs for entry in glob.glob(name)]
    logging.debug("Run dirs: " + ", ".join(alldirs))
    diffs = {}
    lastdiffs = {}
    for checkdir in alldirs:
        diffs[checkdir] = {}
        logging.debug('Checking %s.' % checkdir)
        outfiles = glob.glob(checkdir + '/' + args.jobprefix + '*.o*')
        if not outfiles:
            logging.info('No run result files found for %s.' % checkdir)
            continue
        config = configparser.ConfigParser()
        config.read(checkdir + '/test.ini')
        confdir = checkdir
        if config.getboolean('copy', 'defaults', fallback=False):
            confdir = '../src/'
        try:
            conf = parsedat.readconf(confdir + '/vbfnlo.dat')
            if args.reference is None:
                ref = get_cross_section(checkdir + '/result.out', conf)
            else:
                reffiles = glob.glob(checkdir + '/' + args.reference + '*.o*')
                reffile = sorted(reffiles)[-1]
                ref = get_cross_section(reffile, conf)
        except FileNotFoundError:
            logging.warning("No input or reference run found in %s." % checkdir)
            continue
        except AssertionError:
            logging.warning("No cross section found in reference run in %s." % checkdir)
            continue
        except IndexError:
            logging.warning("No matching reference run found in %s." % checkdir)
            continue
        for outfile in outfiles:
            try:
                run = get_cross_section(outfile, conf)
            except AssertionError as e:
                logging.info('In %s: ' % outfile + str(e))
                continue
            logging.debug("Comparing with %s." % outfile)
            diffs[checkdir][outfile] = get_diff(run, ref)
        if len(diffs[checkdir]) < 1:
            logging.warning("No working output file found out of %s. Skipping %s."
                            % (', '.join(outfiles), checkdir))
            continue
        # get most recent run
        lastdiff = sorted(diffs[checkdir].items())[-1]
        if args.replace:
            shutil.copy(lastdiff[0], checkdir + '/result.out')
        lastdiffs[checkdir] = lastdiff[1]

    runs_by_agreement = sorted(
        lastdiffs.items(),
        key=lambda x: abs(x[1].sigdiff),
    )
    # print(runs_by_agreement)
    for checkdir, rundiff in runs_by_agreement[::-1]:
        for run, diff in diffs[checkdir].items():
            if args.verbose >= 1:
                print(run)
                print_diff(diff)
            else:
                print_short_diff(diff, " " + run)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--verbose',
                        '-v',
                        action='count',
                        help='Generate more output.',
                        default=0)
    parser.add_argument('--output',
                        default="result_10_10.out",
                        help='Filename for run output.')
    parser.add_argument('--jobprefix',
                        default=make_runs.QSUB_PREFIX,
                        help="Prefix of run output.")
    parser.add_argument('--reference',
                        help="Prefix of reference output. Defaults to 'result.out'")
    parser.add_argument('--replace',
                        action='store_true',
                        help="Replace current refrence runs with the new ones.")
    parser.add_argument('testdirs',
                        nargs="*",
                        default="[0-9]*",
                        help='List of test directories, can be shell glob.')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()

    if args.verbose == 1:
        logging.getLogger().setLevel(logging.INFO)
    elif args.verbose >= 2:
        logging.getLogger().setLevel(logging.DEBUG)
    check_results(args)
