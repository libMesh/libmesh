#!/bin/sh
# Tests for ncgen4 using list of test cdl files from the cdl4
# directory, and comparing output to expected results in the expected4
# directory. Same as tst_ncgen4.sh but for classic format.
# Dennie Heimbigner

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

set -e
echo ""
verbose=0
export verbose

echo "*** Performing diff/cycle tests for classic format: k=1"
KFLAG=1 ; export KFLAG
bash ${srcdir}/tst_ncgen4_diff.sh
bash ${srcdir}/tst_ncgen4_cycle.sh
echo "SUCCESS!!"
exit 0


