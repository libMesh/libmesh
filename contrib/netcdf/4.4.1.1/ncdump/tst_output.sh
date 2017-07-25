#!/bin/sh
if test "x$SETX" = x1 ; then echo "file=$0"; set -x ; fi
# This shell script tests the output from several previous tests.

set -e

if test "x$srcdir" = x ; then
srcdir=`pwd`
fi

echo ""
echo "*** Testing ncgen and ncdump test output for classic format."
echo "*** creating ctest1.cdl from ctest0.nc..."
./ncdump -n c1 ctest0.nc | sed 's/e+0/e+/g' > ctest1.cdl
echo "*** creating c0.nc from c0.cdl..."
../ncgen/ncgen -b -o c0.nc $srcdir/../ncgen/c0.cdl
echo "*** creating c1.cdl from c0.nc..."
./ncdump -n c1 c0.nc | sed 's/e+0/e+/g' > c1.cdl
echo "*** comparing ncdump of C program output (ctest1.cdl) with c1.cdl..."
diff -b c1.cdl ctest1.cdl
echo "*** test output for ncdump -k"
test `./ncdump -k c0.nc` = "classic";
../ncgen/ncgen -k `./ncdump -k c0.nc` -b -o c0tmp.nc $srcdir/../ncgen/c0.cdl
cmp c0tmp.nc c0.nc

echo "*** test output for ncdump -x"
echo "*** creating tst_ncml.nc from tst_ncml.cdl"
../ncgen/ncgen -b -o tst_ncml.nc $srcdir/tst_ncml.cdl
echo "*** creating c1.ncml from tst_ncml.nc"
./ncdump -x tst_ncml.nc | sed 's/e-00/e-0/g' > c1.ncml
echo "*** comparing ncdump -x of generated file with ref1.ncml ..."
diff -b c1.ncml $srcdir/ref1.ncml

echo "*** test output for ncdump -s"
echo "*** creating tst_mslp.nc from tst_mslp.cdl"
../ncgen/ncgen -b -o tst_mslp.nc $srcdir/tst_mslp.cdl
echo "*** creating tst_format_att.cdl from tst_mslp.nc"
./ncdump -s tst_mslp.nc > tst_format_att.cdl
echo "*** comparing ncdump -s of generated file with ref_tst_format_att.cdl ..."
diff -b tst_format_att.cdl $srcdir/ref_tst_format_att.cdl

echo "*** All ncgen and ncdump test output for classic format passed!"

echo "*** Testing ncgen and ncdump test output for 64-bit offset format."
echo "*** creating ctest1.cdl from test0_64.nc..."
./ncdump -n c1 ctest0_64.nc | sed 's/e+0/e+/g' > ctest1_64.cdl
echo "*** creating c0.nc from c0.cdl..."
../ncgen/ncgen -k nc6 -b -o c0.nc $srcdir/../ncgen/c0.cdl
echo "*** creating c1.cdl from c0.nc..."
./ncdump -n c1 c0.nc | sed 's/e+0/e+/g' > c1.cdl
echo "*** comparing ncdump of C program output (ctest1_64.cdl) with c1.cdl..."
diff -b c1.cdl ctest1_64.cdl
echo "*** test output for ncdump -k"
test "`./ncdump -k c0.nc`" = "64-bit offset";
../ncgen/ncgen -k nc6 -b -o c0tmp.nc $srcdir/../ncgen/c0.cdl
cmp c0tmp.nc c0.nc

echo "*** test output for ncdump -s"
echo "*** creating tst_mslp_64.nc from tst_mslp.cdl"
../ncgen/ncgen -k nc6 -b -o tst_mslp_64.nc $srcdir/tst_mslp.cdl
echo "*** creating tst_format_att_64.cdl from tst_mslp_64.nc"
./ncdump -s tst_mslp_64.nc | sed 's/e+0/e+/g' > tst_format_att_64.cdl
echo "*** comparing ncdump -s of generated file with ref_tst_format_att_64.cdl ..."
diff -b tst_format_att_64.cdl $srcdir/ref_tst_format_att_64.cdl

echo "*** All ncgen and ncdump test output for 64-bit offset format passed!"
exit 0
