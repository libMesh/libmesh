#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# This shell tests non-sequential reads for bm_file.c

# $Id: run_bm_test2.sh,v 1.3 2007/12/27 03:43:13 ed Exp $

set -e
echo ""

echo "*** Testing the benchmarking program bm_file using start/count/inc..."
${execdir}/bm_file -h -f 3 -o tst_simple_3.nc -c 0:-1:0:10  -t 0:0 -u 0:1 -r 0:1 tst_simple.nc
${NCDUMP} tst_simple.nc > tst_simple.cdl
${NCDUMP} -n tst_simple tst_simple_3.nc > tst_simple_3.cdl
diff tst_simple.cdl tst_simple_3.cdl
echo '*** SUCCESS!!!'

echo "*** Testing the benchmarking program bm_file using start/count/inc with negative increment..."
${execdir}/bm_file -h -f 3 -o tst_simple_3.nc -c 0:-1:0:10  -t 0:9 -u 0:1 -r 0:-1 tst_simple.nc
${NCDUMP} tst_simple.nc > tst_simple.cdl
${NCDUMP} -n tst_simple tst_simple_3.nc > tst_simple_3.cdl
diff tst_simple.cdl tst_simple_3.cdl
echo '*** SUCCESS!!!'

echo "*** Testing the benchmarking program bm_file using start/count/inc with count > 1..."
${execdir}/bm_file -h -f 3 -o tst_simple_3.nc -c 0:-1:0:10  -t 0:0 -u 0:2 -r 0:2 tst_simple.nc
${NCDUMP} tst_simple.nc > tst_simple.cdl
${NCDUMP} -n tst_simple tst_simple_3.nc > tst_simple_3.cdl
diff tst_simple.cdl tst_simple_3.cdl
echo '*** SUCCESS!!!'

echo "*** Testing the benchmarking program bm_file using start/count/inc with uneven count..."
${execdir}/bm_file -h -f 3 -o tst_simple_3.nc -c 0:-1:0:10  -t 0:0 -u 0:4 -r 0:4 tst_simple.nc
${NCDUMP} tst_simple.nc > tst_simple.cdl
${NCDUMP} -n tst_simple tst_simple_3.nc > tst_simple_3.cdl
diff tst_simple.cdl tst_simple_3.cdl
echo '*** SUCCESS!!!'

exit 0
