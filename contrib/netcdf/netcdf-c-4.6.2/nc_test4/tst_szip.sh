#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e

rm -f testnc.h5 testszip.nc szip_dump.cdl

echo "*** Test read of known szip file"
${NCDUMP} ${srcdir}/ref_szip.h5 >szip_dump.cdl
diff -w ${srcdir}/ref_szip.cdl ./szip_dump.cdl

echo "*** Testing tst_szip "
${execdir}/test_szip
echo "***Passed"

echo "*** Testing h5testszip "
${execdir}/h5testszip
echo "***Passed"

echo "*** Testing h5testszip on testszip.nc"
${execdir}/h5testszip ./testszip.nc
echo "***Passed"

rm -f testnc.h5 testszip.nc
rm -f szip_dump.cdl

exit
