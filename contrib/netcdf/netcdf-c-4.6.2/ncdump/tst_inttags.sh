#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

set -e


echo "*** Test integer constant suffixes"
echo "*** creating inttags.nc from inttags.cdl..."
${NCGEN} -lb -o inttags.nc $srcdir/inttags.cdl
echo "*** creating tst_inttags.cdl from inttags.nc..."
${NCDUMP} inttags.nc > tst_inttags.cdl
echo "*** comparing tst_inttags.cdl to ref_inttags.nc..."
diff -b -w tst_inttags.cdl $srcdir/ref_inttags.cdl

rm inttags.nc tst_inttags.cdl

exit 0
