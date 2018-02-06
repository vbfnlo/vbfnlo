#!/usr/bin/env python3

# standard modules
import glob
import tempfile
import logging
from distutils.util import strtobool
import argparse
import configparser
import subprocess

# custom modules
import parsedat
from test_runs import setup_dir, get_run_command, get_cross_section

logging.getLogger().setLevel(logging.WARNING)

QSUB_PREFIX = 'vbfnlo_regress_'


def make_runs(args):
    alldirs = [entry for name in args.testdirs for entry in glob.glob(name)]
    logging.info("Found dirs: " + ", ".join(alldirs))
    for checkdir in alldirs:
        logging.info('Running %s.' % checkdir)
        config = configparser.ConfigParser()
        config.read(checkdir + '/test.ini')
        if args.stat == 'low':
            with tempfile.TemporaryDirectory() as rundir:
                command = get_run_command(rundir, bin=args.bin)
                if command is None:
                    continue
                conf = setup_start_run(config, checkdir, rundir, args)
                logging.debug('Starting %s in %s.', ' '.join(command), rundir)
                try:
                    out = subprocess.check_output(command, cwd=rundir).decode()
                    write_result(conf, checkdir, out, args)
                except subprocess.CalledProcessError as exce:
                    print("Running %s in %s (%s) failed:" % (command, rundir, checkdir))
                    print('returncode: ' + str(exce.returncode))
                    print('stdout: ' + str(exce.output))
        else:
            rundir = checkdir
            command = get_run_command(rundir, bin=args.bin)
            if command is None:
                continue
            conf = setup_start_run(config, checkdir, rundir, args)
            logging.debug('Starting run.')
            send_qsub(command, rundir, args)


def setup_start_run(config, checkdir, rundir, args):
    try:
        setup_dir(config, checkdir, rundir, stat=args.stat)
        conf = parsedat.readconf(rundir + '/vbfnlo.dat')
    except (FileNotFoundError, AssertionError):
        logging.warning("No input or reference run found in %s.", checkdir)
        return None, None
    if strtobool(conf['EWCOR_SWITCH']):
        logging.warning('Not using electroweak corrections at the moment.'
                        ' Skipping run.')
        return None, None
    return conf


def write_result(conf, checkdir, out, args):
    try:
        run = get_cross_section(out, conf)
        assert abs(run.xsec) > 1e-12
    except:
        logging.error(out)
        return
    with open(checkdir + '/' + args.output, 'w') as f:
        f.write(out)
        logging.debug('Finished run.')


def send_qsub(command, rundir, args):
    qsubcmd = ['qsub'] 
    if args.qsubargs:
        qsubcmd += args.qsubargs.split(' ')
    qsubcmd += ['-cwd', '-N', args.jobprefix + rundir.replace('/', '_')] + command
    subprocess.call(qsubcmd, cwd=rundir)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--verbose',
                        '-v',
                        action='count',
                        help='Generate more output.',
                        default=0)
    parser.add_argument('--output',
                        default="result_10_10.out",
                        help='Filename for run output.'
                             ' Default: result_10_10.out')
    parser.add_argument('--bin',
                        help='Binary to use for runs.'
                             'If empty, env VBFNLOPATH is used')
    parser.add_argument('--stat',
                        default="low",
                        help='Choose whether to run with low or high statistics.'
                             ' Default: low')
    parser.add_argument('--jobprefix',
                        default=QSUB_PREFIX,
                        help="Prefix of run output.")
    parser.add_argument('--qsubargs',
                        default="",
                        help="additional arguments to give to qsub.")
    parser.add_argument('testdirs',
                        nargs="*",
                        default=["runs/[0-9]*"],
                        help='List of test directories, can be shell glob.'
                             ' Default: runs/[0-9]*')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    if args.verbose == 1:
        logging.getLogger().setLevel(logging.INFO)
    elif args.verbose >= 2:
        logging.getLogger().setLevel(logging.DEBUG)
    make_runs(args)
