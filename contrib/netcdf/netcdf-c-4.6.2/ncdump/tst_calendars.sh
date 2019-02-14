#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

# This shell script tests ncdump -t option for CF calendar attributes

set -e
echo ""
echo "*** Testing ncdump -t output for times with CF calendar attribute"
echo "*** creating netcdf file tst_calendars.nc from tst_calendars.cdl..."
${NCGEN} -b -o tst_calendars.nc $srcdir/tst_calendars.cdl
echo "*** creating tst_times.cdl from tst_calendars.nc with ncdump -t ..."
${NCDUMP} -n tst_times -t tst_calendars.nc > tst_times.cdl
echo "*** comparing tst_times.cdl with ref_times.cdl..."
diff -b tst_times.cdl $srcdir/ref_times.cdl
echo ""
echo "*** Testing ncdump -t output"
echo "*** creating tst_mud4.cdl from tst_mud4.nc ..."

# Test 360, 365 and 366 day calendars specifically
TSTS="test_360_day_1900 test_365_day_1900 test_366_day_1900"
for t in $TSTS ; do
  rm -f ./${t}.cdl
  echo "create: ${t}.cdl from ${t}.nc"
  ${NCDUMP} -n ${t} -t $srcdir/ref_${t}.nc > ${t}.cdl
  echo "compare: ${t}.cdl ref_${t}.cdl"
  diff -b ${t}.cdl $srcdir/ref_${t}.cdl
  rm -f ${t}.cdl
done

echo "*** All ncdump test output for -t option with CF calendar atts passed!"

exit 0
