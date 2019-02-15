#!/bin/sh
# This script runs UTF-8 tests for netCDF-4.
# Ward Fisher, Dennis Heimbigner, Ed Hartnett

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh
set -e

echo ""
echo "*** Testing netcdf-4 file with utf8 characters..."

rm -f tst_utf8_nc4.nc tst_utf8_nc4.cdl
${NCGEN} -4 -b -o tst_utf8_nc4.nc ${srcdir}/ref_tst_utf8_4.cdl
${NCDUMP} -n 'utf8' tst_utf8_nc4.nc > tst_utf8_nc4.cdl
diff -b -w tst_utf8_nc4.cdl ${srcdir}/ref_tst_utf8_4.cdl

echo "*** NetCDF-4 UTF8 testing passed!"
exit 0
