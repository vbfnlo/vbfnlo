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

To update the result files, run `make_runs.py`. See `make_runs.py --help` for possible arguments.

To add new runs, add a directory with input files and output in `result.out`. To generate the output `make_runs.py` can be used.
To have the tests run automatically there must be a low-stat file available as `results_10_10.out`.

To run the tests with full statistics, you probably need a computing cluster.
`make_runs.py --stat=high` starts the runs using qsub. Adjust the code to use other cluster systems if needed.
The results can be collected and compared to the existing runs using `check_results.py`.
