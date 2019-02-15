#!/bin/sh
# This shell script runs the ncdump tests.
# $Id: run_tests.sh,v 1.18 2010/05/19 13:43:39 ed Exp $

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e
echo ""
echo "*** Testing ncgen and ncdump using some test CDL files."
echo "*** creating tst_small.nc from ref_tst_small.cdl..."
${NCGEN} -b -o tst_small.nc $srcdir/ref_tst_small.cdl
echo "*** creating tst_small.cdl from tst_small.nc..."
${NCDUMP} tst_small.nc > tst_small.cdl
diff -b -w tst_small.cdl $srcdir/ref_tst_small.cdl

echo "*** creating test0_ncdump.nc from test0.cdl..."
${NCGEN} -o test0_ncdump.nc -b $srcdir/test0.cdl
echo "*** creating test1_ncdump.cdl from test0_ncdump.nc..."
${NCDUMP} -n test1_ncdump test0_ncdump.nc > test1_ncdump.cdl
echo "*** creating test1_ncdump.nc from test1_ncdump.cdl..."
${NCGEN} -b test1_ncdump.cdl
echo "*** creating test2_ncdump.cdl from test1_ncdump.nc..."
${NCDUMP} test1_ncdump.nc > test2_ncdump.cdl
echo "*** checking that test1_ncdump.cdl and test2_ncdump.cdl are the same..."
diff -b -w test1_ncdump.cdl test2_ncdump.cdl

echo "*** All tests of ncgen and ncdump using test0.cdl passed!"
exit 0
