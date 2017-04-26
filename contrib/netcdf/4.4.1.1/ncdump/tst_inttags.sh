#!/bin/sh
if test "x$SETX" = x1 ; then echo "file=$0"; set -x ; fi

set -e

if test "x$srcdir" = x ; then
srcdir=`pwd`
fi

echo "*** Test integer constant suffixes"
echo "*** creating inttags.nc from inttags.cdl..."
../ncgen/ncgen -lb -o inttags.nc $srcdir/inttags.cdl
echo "*** creating tst_inttags.cdl from inttags.nc..."
./ncdump inttags.nc > tst_inttags.cdl
echo "*** comparing tst_inttags.cdl to ref_inttags.nc..."
diff -b -w tst_inttags.cdl $srcdir/ref_inttags.cdl

rm inttags.nc tst_inttags.cdl

exit 0
