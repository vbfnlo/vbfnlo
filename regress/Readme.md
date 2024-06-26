Regress tests
=============

This folder collects example runs in the subdirectory `runs`.
Those are good starting points to use for new analysis runs and to compare the output of modified VBFNLO versions.

New versions of VBFNLO are tested to reproduce the output of the existing runs.

The tests can be run in 2 modes:
* low statistics (2**10 phase space points, ~5 minutes)
* full/high statistics (>2**20 points usually, few days)

The low statistics runs can be started via `make check`.


For developers
==============

`make check` uses the following helpers, that should not be called directly:
* `makecheck.sh`: For argument and dependency checks
* `test_runs.py`: Called by py.test to run the tests and check the output

'make check' only compares the low statistic runs 'result_10_10.out' and checks for perfect numerical agreement between these files and the results generated by the respective vbfnlo installation on the fly.
If such is not found, a 'FAIL' is returned. A test is skipped if no low stat result file is found.

To update the result files, run `make_runs.py`. See `make_runs.py --help` for possible arguments.

To add new runs, add a directory with input files and output in `result.out`. To generate the output `make_runs.py` can be used.
To have the tests run automatically there must be a low-stat file available as `results_10_10.out`.

To run the tests with full statistics, you probably need a computing cluster.
`make_runs.py --stat=high` starts the runs using sbatch. Adjust the code to use other cluster systems if needed.
The results can be collected and compared to the existing runs using 'check_results.py'.
Here, use --stat to define which results should be compared the high stat result.out or other references via --reference. 
Small deviations are expected here, e.g. depending on the use of mpirun or the version of FeynHiggs used.

Additionally, 'set-vbfnlo-config.py' can be used to set values for LO/NLO points and iterations in the vbfnlo.dat files for all runs subdirectories.
