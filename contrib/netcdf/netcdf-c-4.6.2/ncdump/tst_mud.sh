#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

# This shell script tests ncdump and ncgen on netCDF-4 variables with multiple 
# unlimited dimensions.

set -e

echo ""
echo "*** Testing ncdump output for multiple unlimited dimensions"
echo "*** creating netcdf file tst_mud4.nc from ref_tst_mud4.cdl ..."
${NCGEN} -4 -b -o tst_mud4.nc $srcdir/ref_tst_mud4.cdl
echo "*** creating tst_mud4.cdl from tst_mud4.nc ..."
${NCDUMP} tst_mud4.nc > tst_mud4.cdl
# echo "*** comparing tst_mud4.cdl with ref_tst_mud4.cdl..."
diff -b tst_mud4.cdl $srcdir/ref_tst_mud4.cdl
# echo "*** comparing annotation from ncdump -bc tst_mud4.nc with expected output..."
${NCDUMP} -bc tst_mud4.nc > tst_mud4-bc.cdl
diff -b tst_mud4-bc.cdl $srcdir/ref_tst_mud4-bc.cdl
# Now test with char arrays instead of ints
echo "*** creating netcdf file tst_mud4_chars.nc from ref_tst_mud4_chars.cdl ..."
${NCGEN} -4 -b -o tst_mud4_chars.nc $srcdir/ref_tst_mud4_chars.cdl
echo "*** creating tst_mud4_chars.cdl from tst_mud4_chars.nc ..."
${NCDUMP} tst_mud4_chars.nc > tst_mud4_chars.cdl
# echo "*** comparing tst_mud4_chars.cdl with ref_tst_mud4_chars.cdl..."
diff -b tst_mud4_chars.cdl $srcdir/ref_tst_mud4_chars.cdl
exit 0
# unused
# echo "*** comparing annotation from ncdump -bc tst_mud4_chars.nc with expected output..."
${NCDUMP} -bc tst_mud4_chars.nc > tst_mud4_chars-bc.cdl
# diff -b tst_mud4_chars-bc.cdl $srcdir/ref_tst_mud4_chars-bc.cdl
echo "*** All ncdump test output for multiple unlimited dimensions passed!"
exit 0
