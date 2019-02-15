#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e

echo "*** Test netcdf-4 integer constant suffixes"

echo "*** creating inttags4.nc from inttags4.cdl..."
${NCGEN} -lb -k nc4 -o inttags4.nc $srcdir/inttags4.cdl
echo "*** creating tst_inttags4.cdl from inttags4.nc..."
${NCDUMP} inttags4.nc > tst_inttags4.cdl
echo "*** comparing tst_inttags4.cdl to ref_inttags4.nc..."
diff -b -w tst_inttags4.cdl $srcdir/ref_inttags4.cdl

rm inttags4.nc tst_inttags4.cdl

exit 0
