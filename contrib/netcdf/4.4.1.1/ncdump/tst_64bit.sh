#!/bin/sh
if test "x$SETX" = x1 ; then echo "file=$0"; set -x ; fi
# This shell script runs the ncdump tests.
# $Id: tst_64bit.sh,v 1.9 2006/03/04 18:50:15 ed Exp $

echo ""
echo "*** Testing ncgen and ncdump with 64-bit offset format."
set -e
echo "*** creating test0.nc from test0.cdl..."
../ncgen/ncgen -b -k2 $srcdir/test0.cdl
echo "*** creating test1.cdl from test0.nc..."
./ncdump -n test1 test0.nc > test1.cdl
echo "*** creating test1.nc from test1.cdl..."
../ncgen/ncgen -b -k2 test1.cdl
echo "*** creating test2.cdl from test1.nc..."
./ncdump test1.nc > test2.cdl
cmp test1.cdl test2.cdl
echo "*** All ncgen and ncdump with 64-bit offset format tests passed!"

echo ""
echo "*** Testing ncgen and ncdump with CDF5 format."
set -e
echo "*** creating test0.nc from test0.cdl..."
../ncgen/ncgen -b -k5 $srcdir/test0.cdl
echo "*** creating test1.cdl from test0.nc..."
./ncdump -n test1 test0.nc > test1.cdl
echo "*** creating test1.nc from test1.cdl..."
../ncgen/ncgen -b -k5 test1.cdl
echo "*** creating test2.cdl from test1.nc..."
./ncdump test1.nc > test2.cdl
cmp test1.cdl test2.cdl
echo "*** All ncgen and ncdump with CDF5 format tests passed!"

exit 0
