#!/bin/sh
if test "x$SETX" = x1 ; then echo "file=$0"; set -x ; fi
# This shell script tests the output several previous tests.
# $Id: tst_output.sh,v 1.17 2010/05/14 16:21:15 ed Exp $


if test "x$srcdir" = x ; then
srcdir="."
fi

ECODE=0

echo ""
echo "*** Testing extended file format output."
set -e
echo "Test extended format output for a netcdf-3 file"
rm -f tmp
../ncgen/ncgen -k nc3 -b -o ./test.nc $srcdir/ref_tst_small.cdl
./ncdump -K test.nc >tmp
if ! grep 'classic mode=00000000' <tmp ; then
echo "*** Fail: extended format for a classic file"
ECODE=1
fi

echo "Test extended format output for a 64-bit offset netcdf-3 file"
rm -f tmp
../ncgen/ncgen -k nc6 -b -o ./test.nc $srcdir/ref_tst_small.cdl
./ncdump -K test.nc >tmp
if ! grep '64-bit offset mode=00000200' <tmp ; then
echo "*** Fail: extended format for a 64-bit classic file"
ECODE=1
fi

echo "Test extended format output for a 64-bit CDF-5 classic file"
rm -f tmp
../ncgen/ncgen -k5 -b -o ./test.nc $srcdir/ref_tst_small.cdl
./ncdump -K test.nc >tmp
if ! grep -F '64-bit data mode=00000020' <tmp ; then
echo "*** Fail: extended format for a 64-bit CDF-5 classic file"
ECODE=1
fi

rm -f tmp test.nc

exit $ECODE


