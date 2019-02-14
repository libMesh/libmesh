#!/bin/sh
# Tests for ncgen4 using list of test cdl files from the cdl4
# directory, and comparing output to expected results in the expected4
# directory. Note that these tests are run for classic files in
# tst_ncgen4_classic.sh
# Dennis Heimbigner

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e

# To add a new test,
# 1. put the .cdl file in the 'cdl4' directory
# 2. put the result of running ncgen then ncdump
#    into the directory 'expected4' as .dmp
# 3. Modify the file tst_ncgen_shared.sh to add
#    the test to the end of the TESTS4 variable
# 4. Add the new files into cdl4/Makefile.am
#    and expected4/Makefile.am

verbose=1
export verbose

KFLAG=3 ; export KFLAG
echo "*** Performing diff tests: k=3"
bash  ${srcdir}/tst_ncgen4_diff.sh
echo "*** Performing cycle tests: k=3"
bash  ${srcdir}/tst_ncgen4_cycle.sh
KFLAG=4 ; export KFLAG
echo "*** Performing diff tests: k=4"
bash  ${srcdir}/tst_ncgen4_diff.sh
echo "*** Performing cycle tests: k=4"
bash  ${srcdir}/tst_ncgen4_cycle.sh
rm -rf ${RESULTSDIR}
echo "SUCCESS!!"
exit 0
