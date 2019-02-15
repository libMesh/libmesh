#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

# This shell script tests the output several previous tests.

ECODE=0
echo ""
echo "*** Testing extended file format output."
set -e
echo "Test extended format output for a netcdf-4 file"
rm -f tmp_tst_formatx4
${NCGEN} -k nc4 -b -o ./tst_formatx4.nc $srcdir/ref_tst_small.cdl
${NCDUMP} -K tst_formatx4.nc >tmp_tst_formatx4
if ! grep 'HDF5 mode=00001000' <tmp_tst_formatx4 ; then
echo "*** Fail: extended format for a netcdf-4 file"
ECODE=1
fi

echo "Test extended format output for a classic netcdf-4 file"
rm -f tmp_tst_formatx4
${NCGEN} -k nc7 -b -o ./tst_formatx4.nc $srcdir/ref_tst_small.cdl
${NCDUMP} -K tst_formatx4.nc >tmp_tst_formatx4
if ! grep 'HDF5 mode=00001000' <tmp_tst_formatx4 ; then
echo "*** Fail: extended format for a classic netcdf-4 file"
ECODE=1
fi

rm -f tmp_tst_formatx4 tst_formatx4.nc

exit $ECODE

