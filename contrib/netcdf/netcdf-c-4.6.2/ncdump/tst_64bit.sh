#!/bin/sh

set -e

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# This shell script runs the ncdump tests.
# get some config.h parameters
if test -f ${top_builddir}/config.h ; then
  if fgrep -e '#define ENABLE_CDF5 1' ${top_builddir}/config.h >/dev/null ; then
    CDF5=1
  else
    CDF5=0
  fi
else
  echo "Cannot locate config.h"
  exit 1
fi

echo ""
echo "*** Testing ncgen and ncdump with 64-bit offset format."
set -e
echo "*** creating test0_offset.nc from test0.cdl..."
${NCGEN} -b -k2 -o test0_offset.nc $srcdir/test0.cdl
echo "*** creating test1_offset.cdl from test0_offset.nc..."
${NCDUMP} -n test1_offset test0_offset.nc > test1_offset.cdl
echo "*** creating test1_offset.nc from test1_offset.cdl..."
${NCGEN} -b -k2 -o test1_offset.nc test1_offset.cdl
echo "*** creating test2_offset.cdl from test1.nc..."
${NCDUMP} test1_offset.nc > test2_offset.cdl
cmp test1_offset.cdl test2_offset.cdl
echo "*** All ncgen and ncdump with 64-bit offset format tests passed!"


if test "x$CDF5" = x1 ; then
echo ""
echo "*** Testing ncgen and ncdump with CDF5 format."
set -e
echo "*** creating test0_cdf5.nc from test0.cdl..."
${NCGEN} -b -k5 -o test0_cdf5.nc $srcdir/test0.cdl
echo "*** creating test1_cdf5.cdl from test0_cdf5.nc..."
${NCDUMP} -n test1_cdf5 test0_cdf5.nc > test1_cdf5.cdl
echo "*** creating test1_cdf5.nc from test1_cdf5.cdl..."
${NCGEN} -b -k5 -o test1_cdf5.nc test1_cdf5.cdl
echo "*** creating test2_cdf5.cdl from test1_cdf5.nc..."
${NCDUMP} test1_cdf5.nc > test2_cdf5.cdl
cmp test1_cdf5.cdl test2_cdf5.cdl
echo "*** All ncgen and ncdump with CDF5 format tests passed!"
fi

exit 0
