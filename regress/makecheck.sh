#!/bin/bash

if ! which py.test > /dev/null; then
    echo "py.test is needed to run the VBFNLO tests. Please install py.test via a package manager or using pip."
    exit 1
fi

if [ ! -x $VBFNLOPATH/vbfnlo ]; then
    echo "Can't find vbfnlo in \"$(realpath $VBFNLOPATH)\"."
    echo "Either the path is set incorrectly or the compilation did not finish."
    exit 1
fi


if [ ! -e $TESTPATH/100_Hjj ]; then
    echo "Can't find tests in \"$(realpath $TESTPATH)\"."
    echo "Either the path is set incorrectly or the compilation did not finish."
    exit 1
fi

py.test -vv -rsx $TESTPATH/test_runs.py
