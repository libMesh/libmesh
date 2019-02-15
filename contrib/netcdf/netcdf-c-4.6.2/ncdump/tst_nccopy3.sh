#!/bin/sh
# For a netCDF-3 build, test nccopy on netCDF files in this
# directory. This test depends on a bunch of other ncdump tests
# running first, to produce the data files that are used to test
# nccopy.
# Dennis Heimbigner

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e
echo ""

# get some config.h parameters
if test -f ${top_builddir}/config.h ; then
  if fgrep -e '#define ENABLE_CDF5 1' ${top_builddir}/config.h >/dev/null ; then
    HAVE_CDF5=1
  else
    HAVE_CDF5=0
  fi
else
  echo "Cannot locate config.h"
  exit 1
fi

TESTFILES='tst_output_c0 tst_output_c0tmp ctest0 ctest0_64 test0_offset test1_offset
 tst_calendars tst_mslp tst_mslp_64 tst_ncml tst_small tst_utf8'

if test "x$HAVE_CDF5" = x1 ; then
    TESTFILES="$TESTFILES small small2"
fi


echo "*** Testing netCDF-3 features of nccopy on ncdump/*.nc files"
for i in $TESTFILES ; do
    echo "*** Testing nccopy $i.nc nccopy3_copy_of_$i.nc ..."
ls -l $i.nc
    ${NCCOPY} $i.nc nccopy3_copy_of_$i.nc
    ${NCDUMP} -n nccopy3_copy_of_$i $i.nc > tmp_tst_nccopy3.cdl
    ${NCDUMP} nccopy3_copy_of_$i.nc > nccopy3_copy_of_$i.cdl
    diff nccopy3_copy_of_$i.cdl tmp_tst_nccopy3.cdl
    rm nccopy3_copy_of_$i.nc nccopy3_copy_of_$i.cdl tmp_tst_nccopy3.cdl
done
echo "*** Testing nccopy -u"
${NCGEN} -b $srcdir/tst_brecs.cdl
# convert record dimension to fixed-size dimension
$NCCOPY -u tst_brecs.nc nccopy3_copy_of_tst_brecs.nc
${NCDUMP} -n nccopy3_copy_of_tst_brecs tst_brecs.nc | sed '/ = UNLIMITED ;/s/\(.*\) = UNLIMITED ; \/\/ (\(.*\) currently)/\1 = \2 ;/' > tmp_tst_nccopy3.cdl
${NCDUMP} nccopy3_copy_of_tst_brecs.nc >  nccopy3_copy_of_tst_brecs.cdl
diff -b nccopy3_copy_of_tst_brecs.cdl tmp_tst_nccopy3.cdl
rm nccopy3_copy_of_tst_brecs.cdl tmp_tst_nccopy3.cdl tst_brecs.nc nccopy3_copy_of_tst_brecs.nc

echo "*** All nccopy tests passed!"
exit 0
