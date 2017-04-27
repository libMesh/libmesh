#!/bin/sh
if test "x$SETX" = x1 ; then echo "file=$0"; set -x ; fi
#
# Moving some netcdf-4 only tests here, out of tst_nccopy and run_utf8_tests.
# Without this, the tests fail when netcdf-4 is disabled.

set -e
if test "x$srcdir" = "x"; then
    srcdir=`dirname $0`;
fi
export srcdir
echo ""

rm -f utf8.nc utf8.cdl
echo "*** creating enhanced file with utf8 characters..."
../ncgen/ncgen -4 -b -o utf8.nc ${srcdir}/ref_tst_utf8_4.cdl
echo "*** dump and compare utf8 output..."
./ncdump utf8.nc > utf8.cdl
diff -b -w utf8.cdl ${srcdir}/ref_tst_utf8_4.cdl
