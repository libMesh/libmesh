#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# This shell script tests the output from several previous tests.
set -e

echo ""
echo "*** Testing ncgen and ncdump test output for classic format."

echo "*** Testing that ncgen produces correct C code from c0.cdl."
${execdir}/ref_ctest
${NCGEN} -lc -o ctest0.nc $srcdir/../ncgen/c0.cdl > tst_output_ctest.c
diff -b tst_output_ctest.c $srcdir/ref_ctest.c

echo "*** creating ctest1.cdl from tst_output_ctest0.nc..."
${NCDUMP} -n c1 ${builddir}/ctest0.nc | sed 's/e+0/e+/g' > tst_output_ctest1.cdl
echo "*** creating tst_output_c0.nc from c0.cdl..."
${NCGEN} -b -o tst_output_c0.nc ${ncgenc0}
echo "*** creating tst_output_c1.cdl from tst_output_c0.nc..."
${NCDUMP} -n c1 ${builddir}/tst_output_c0.nc | sed 's/e+0/e+/g' > tst_output_c1.cdl
echo "*** comparing tst_output_c1.cdl with ref_ctest1_nc4c.cdl..."
diff -b tst_output_c1.cdl $srcdir/ref_ctest1_nc4c.cdl
echo "*** comparing ncdump of C program output (ctest1.cdl) with c1.cdl..."
diff -b tst_output_c1.cdl tst_output_ctest1.cdl
echo "*** test output for ncdump -k"
KIND=`${NCDUMP} -k tst_output_c0.nc`
test "$KIND" = "classic";
${NCGEN} -k $KIND -b -o tst_output_c0tmp.nc ${ncgenc0}
cmp tst_output_c0tmp.nc tst_output_c0.nc

echo "*** test output for ncdump -x"
echo "*** creating tst_ncml.nc from tst_ncml.cdl"
${NCGEN} -b -o tst_ncml.nc $srcdir/tst_ncml.cdl
echo "*** creating c1.ncml from tst_ncml.nc"
${NCDUMP} -x tst_ncml.nc | sed 's/e-00/e-0/g' > c1.ncml
echo "*** comparing ncdump -x of generated file with ref1.ncml ..."
diff -b c1.ncml $srcdir/ref1.ncml

echo "*** test output for ncdump -s"
echo "*** creating tst_mslp.nc from tst_mslp.cdl"
${NCGEN} -b -o tst_mslp.nc $srcdir/tst_mslp.cdl
echo "*** creating tst_format_att.cdl from tst_mslp.nc"
${NCDUMP} -s tst_mslp.nc > tst_format_att.cdl
echo "*** comparing ncdump -s of generated file with ref_tst_format_att.cdl ..."
diff -b tst_format_att.cdl $srcdir/ref_tst_format_att.cdl

echo "*** All ncgen and ncdump test output for classic format passed!"

echo "*** Testing that ncgen with c0.cdl for 64-bit offset format."
${execdir}/ref_ctest64
${NCGEN}  -k2 -lc -o ctest0_64.nc $srcdir/../ncgen/c0.cdl > tst_output_ctest64.c
diff -b tst_output_ctest64.c $srcdir/ref_ctest64.c

echo "*** Testing ncgen and ncdump test output for 64-bit offset format."
echo "*** creating ctest1_64.cdl from test0_64.nc..."
${NCDUMP} -n c1 ctest0_64.nc | sed 's/e+0/e+/g' > tst_output_ctest1_64.cdl
echo "*** creating tst_output_c0_64.nc from c0.cdl..."
${NCGEN} -k nc6 -b -o tst_output_c0_64.nc ${ncgenc0}
echo "*** creating tst_output_c1_64.cdl from tst_output_c0_64.nc..."
${NCDUMP} -n c1 tst_output_c0_64.nc | sed 's/e+0/e+/g' > tst_output_c1_64.cdl
echo "*** comparing ncdump of C program output (ctest1_64.cdl) with tst_output_c1_64.cdl..."
diff -b tst_output_c1_64.cdl tst_output_ctest1_64.cdl
echo "*** test output for ncdump -k"
test "`${NCDUMP} -k tst_output_c0_64.nc`" = "64-bit offset";
${NCGEN} -k nc6 -b -o tst_output_c0_64_tmp.nc ${ncgenc0}
cmp tst_output_c0_64_tmp.nc tst_output_c0_64.nc

echo "*** test output for ncdump -s"
echo "*** creating tst_mslp_64.nc from tst_mslp.cdl"
${NCGEN} -k nc6 -b -o tst_mslp_64.nc ${srcdir}/tst_mslp.cdl
echo "*** creating tst_format_att_64.cdl from tst_mslp_64.nc"
${NCDUMP} -s tst_mslp_64.nc | sed 's/e+0/e+/g' > tst_format_att_64.cdl
echo "*** comparing ncdump -s of generated file with ref_tst_format_att_64.cdl ..."
diff -b tst_format_att_64.cdl $srcdir/ref_tst_format_att_64.cdl

echo "*** All ncgen and ncdump test output for 64-bit offset format passed!"
exit 0
