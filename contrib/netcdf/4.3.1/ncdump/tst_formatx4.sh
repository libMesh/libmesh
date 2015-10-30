#!/bin/sh
# This shell script tests the output several previous tests.
# $Id: tst_output.sh,v 1.17 2010/05/14 16:21:15 ed Exp $

if test "x$srcdir" = x ; then
srcdir="."
fi

ECODE=0
echo ""
echo "*** Testing extended file format output."
set -e
echo "Test extended format output for a netcdf-4 file"
rm -f tmp
../ncgen/ncgen -k3 -b -o ./test.nc $srcdir/ref_tst_small.cdl
./ncdump -K test.nc >tmp
if ! fgrep 'HDF5 mode=00001000' <tmp ; then
echo "*** Fail: extended format for a netcdf-4 file"
ECODE=1
fi

echo "Test extended format output for a classic netcdf-4 file"
rm -f tmp
../ncgen/ncgen -k4 -b -o ./test.nc $srcdir/ref_tst_small.cdl
./ncdump -K test.nc >tmp
if ! fgrep 'HDF5 mode=00001000' <tmp ; then
echo "*** Fail: extended format for a classic netcdf-4 file"
ECODE=1
fi

rm -f tmp test.nc

exit $ECODE

