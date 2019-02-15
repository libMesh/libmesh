#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

# This shell script tests ncdump -g option for specifying groups for
# which data is to be output.

set -e
echo ""
echo "*** Testing ncdump -g output for specifying group subsets"
echo "*** creating netcdf file tst_grp_spec.nc from ref_tst_grp_spec0.cdl..."
${NCGEN} -knc4 -o tmp_all.nc -b $srcdir/ref_tst_grp_spec0.cdl
echo "*** creating tmp_subset.cdl from tmp_all.nc with ncdump -g ..."
${NCDUMP} -g g1,g4 -n tmp_subset tmp_all.nc > tmp_subset.cdl
echo "*** comparing tmp_subset.cdl with ref_tst_grp_spec.cdl..."
diff -b tmp_subset.cdl $srcdir/ref_tst_grp_spec.cdl

echo "*** All ncdump test output for -g option passed!"
exit 0
