#!/usr/bin/env python3

# standard modules
import glob
import os
import tempfile
import re
import collections
import shutil
import subprocess
import logging
from distutils.util import strtobool
import configparser
import math
import pytest

# custom modules
import parsedat

alldirs = glob.glob(os.path.dirname(os.path.realpath(__file__))+'/[0-9]*')
envtestdirs = os.getenv('TESTDIRS', None)
if envtestdirs:
    alldirs = envtestdirs.split(',')
Run = collections.namedtuple('Run', ['xsec', 'err'])
Diff = collections.namedtuple('Diff', ['ref', 'run', 'diff', 'err', 'reldiff', 'sigdiff'])


@pytest.mark.parametrize('checkdir', alldirs)
def test_run(checkdir):
    "run a specific directory"

    if not os.path.exists(checkdir + '/result_10_10.out'):
        pytest.skip("No low-statistic results for %s." % checkdir)

    with tempfile.TemporaryDirectory() as rundir:
        config = configparser.ConfigParser()
        config.read(checkdir + '/test.ini')
        setup_dir(config, checkdir, rundir)

        conf = parsedat.readconf(rundir + '/vbfnlo.dat')
        command = get_run_command(rundir, int(conf['PROCESS']))
        try:
            out = subprocess.check_output(command, cwd=rundir).decode()
        except subprocess.CalledProcessError as exce:
            if exce.returncode == 2:
                print("Process not compiled in.")
                if os.getenv('RUNALLPROCESSES', None) == 'yes':
                    assert False
                else:
                    pytest.skip("Process not compiled. Compile with --enable-processes=all to run the complete testsuite.")
            print("Running %s in %s (%s) failed:" % (command, rundir, checkdir))
            print('returncode: ' + str(exce.returncode))
            print('stdout: ' + str(exce.output))
            assert False

        run = get_cross_section(out, conf)
        resultfiles = glob.glob(checkdir + '/result_10_10*.out')
        referenceruns = [get_cross_section(fname, conf) for fname in resultfiles]
        diffs = [get_diff(run, ref) for ref in referenceruns]

        print(out)  # get's displayed on failure
        for d in diffs:
            print_diff(d)

        if config.getboolean('check', 'output', fallback=False):
            referenceout = []
            for fname in resultfiles:
                with open(fname) as f:
                    referenceout.append(f.read())
            compare_output(out, referenceout)


        # pass if any of the reference runs agrees with the result
        assert any([abs(d.reldiff) < 1e-4 for d in diffs])

        # TODO
        # for runs with more statstics use
        # assert abs(sigdiff) < 3.0 and abs(reldiff) < 1e-3


def setup_dir(config, checkdir, rundir, stat='low'):
    # check for existance of vbfnlo.dat except for [require] vbfnlo.dat = False in config

    if config.getboolean('copy', 'defaults', fallback=False):
        runfiles = glob.glob('../src/*.dat')
        logging.debug('Copying default %s.' % ', '.join(runfiles))
        for fname in runfiles:
            shutil.copy(fname, rundir)
    else:
        assert os.path.exists(checkdir + '/vbfnlo.dat'),\
            'No vbfnlo.dat found in %s' % checkdir

    if checkdir != rundir:
        runfiles = glob.glob(checkdir + '/*.dat')
        if runfiles:
            logging.info('Copying %s.' % ', '.join(runfiles))
        for fname in runfiles:
            shutil.copy(fname, rundir)
    if stat == 'low':
        # TODO: loop with increading statistics
        vars = {
            'LO_POINTS': '10',
            'NLO_POINTS': '10',
            'LO_ITERATIONS': '1',
            'NLO_ITERATIONS': '1',
        }
        for key, val in vars.items():
            parsedat.replacemachine(rundir + '/vbfnlo.dat',
                                    key, val, keep=False)


def get_run_command(rundir, procid):
    cmd = 'vbfnlo'

    vbfnlobin = os.getenv('VBFNLOPATH', None)
    if vbfnlobin:
        if os.path.isabs(vbfnlobin):
            cmd = [os.path.join(vbfnlobin, cmd)]
        else:
            srcdir = os.path.abspath(os.path.curdir)
            cmd = [os.path.join(srcdir, vbfnlobin, cmd)]

    cmdprefix = os.getenv('CMDPREFIX', None)
    if cmdprefix:
        cmd = cmdprefix.split(' ') + cmd
    return cmd


# regular expressions to extract cross section from stdout output
nlo_result_re = re.compile('final result at NLO\s*\n\s*'
                           'sigma = \s*(?P<xsec>.*?)\s*'
                           '\+-\s*(?P<err>.*?)\s*fb\s*(?P<percent>.*?)\s*%')
lo_result_re = re.compile('result \(LO\):\s*(?P<xsec>.*?)\s*'
                          '\+-\s*(?P<err>.*?)\s*fb\s*(?P<percent>.*?)\s*%')


def get_cross_section(source, conf):
    if len(source.splitlines()) == 1:
        # assume source is a filename
        with open(source) as f:
            out = f.read()
    else:
        # assume source contains the output lines
        out = source
    nlo = strtobool(conf['NLO_SWITCH'])
    if nlo:
        match = nlo_result_re.search(out)
    else:
        match = lo_result_re.search(out)
    assert match, "Could not find cross section in output"
    xsec = float(match.groupdict()['xsec'])
    err = float(match.groupdict()['err'])
    return Run(xsec, err)


def get_diff(run, reference):
    diff = run.xsec - reference.xsec
    err = math.sqrt(run.err**2 + reference.err**2)
    reldiff = (run.xsec / reference.xsec - 1)
    sigdiff = diff / err
    return Diff(reference, run, diff, err, reldiff, sigdiff)


def print_diff(d):
    print("Reference: %.2e +- %.2e fb" % (d.ref.xsec, d.ref.err))
    print("New Run:   %.2e +- %.2e fb" % (d.run.xsec, d.run.err))
    print("Diff:      %.2e +- %.2e fb (%+.3f %%  %.2f sigma)" %
          (d.diff, d.err, d.reldiff * 100, d.sigdiff))
    # print(d.run.err)
    # print(d.run.err.__class__)
    # print(d.ref.err)
    # print(d.ref.err.__class__)
    print("Errordiff: %+.2f" % (d.run.err / d.ref.err - 1))


def print_short_diff(d, note=""):
    print("%.2f sigma = %+.3f %% (%.2e +- %.2e fb)" %
          (abs(d.sigdiff), d.reldiff * 100, d.diff, d.err) + note)


def compare_output(out, referenceout):
    """compare output text with reference run

    remove lines containing VBFNLO, which also contain a git hash"""

    outclean = [line for line in out.splitlines() if 'VBFNLO' not in line]
    foundmatch = False
    for r in referenceout:
        refclean = [line for line in r.splitlines() if 'VBFNLO' not in line]
        if '\n'.join(outclean) == '\n'.join(refclean):
            foundmatch = True
            return
    if not foundmatch:
        assert '\n'.join(outclean) == '\n'.join(refclean)

if __name__ == "__main__":
    print("This script should not be run directly, but by calling py.test in the regress folder.")
