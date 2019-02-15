#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

# This shell script runs an ncgen buf in handling character _Fillvalue.
set -e
echo ""
echo "*** Testing char _Fillvalue bug fix."

# echo "*** Create tst_charfill.nc from tst_charfill.cdl..."
rm -f tst_charfill.nc tmp_tst_charfill.cdl
${NCGEN} -b -o tst_charfill.nc $srcdir/tst_charfill.cdl
# echo "*** dumping tst_charfill.nc to tmp_tst_charfill.cdl..."
${NCDUMP} tst_charfill.nc > tmp_tst_charfill.cdl
# echo "*** comparing tmp_tst_charfill.cdl with ref_tst_charfill.cdl..."
diff -b tmp_tst_charfill.cdl $srcdir/ref_tst_charfill.cdl

echo "*** All char _Fillvalue  bug test for netCDF-4 format passed!"
exit 0
