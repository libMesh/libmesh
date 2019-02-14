#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# This shell script runs the ncdump tests.
# get some config.h parameters
if test -f ${top_builddir}/config.h ; then
  if fgrep -e '#define ENABLE_CDF5 1' ${top_builddir}/config.h >/dev/null ; then
    ENABLE_CDF5=1
  else
    ENABLE_CDF5=0
  fi
else
  echo "Cannot locate config.h"
  exit 1
fi

# This shell script tests the output several previous tests.

ECODE=0

echo ""
echo "*** Testing extended file format output."
set -e
echo "Test extended format output for a netcdf-3 file"
rm -f tmp_tst_formatx3
${NCGEN} -k nc3 -b -o ./tst_formatx3.nc $srcdir/ref_tst_small.cdl
${NCDUMP} -K tst_formatx3.nc >tmp_tst_formatx3
if ! grep 'classic mode=00000000' <tmp_tst_formatx3 ; then
echo "*** Fail: extended format for a classic file"
ECODE=1
fi

echo "Test extended format output for a 64-bit offset netcdf-3 file"
rm -f tmp_tst_formatx3
${NCGEN} -k nc6 -b -o ./tst_formatx3.nc $srcdir/ref_tst_small.cdl
${NCDUMP} -K tst_formatx3.nc >tmp_tst_formatx3
if ! grep '64-bit offset mode=00000200' <tmp_tst_formatx3 ; then
echo "*** Fail: extended format for a 64-bit classic file"
ECODE=1
fi


# Only do following test if ENABLE_CDF5 is true.

if test "x$ENABLE_CDF5" = x1 ; then
    echo "Test extended format output for a 64-bit CDF-5 classic file"
    rm -f tmp_tst_formatx3
    ${NCGEN} -k5 -b -o ./tst_formatx3.nc $srcdir/ref_tst_small.cdl
    ${NCDUMP} -K tst_formatx3.nc >tmp_tst_formatx3
    if ! grep -F '64-bit data mode=00000020' <tmp_tst_formatx3 ; then
        echo "*** Fail: extended format for a 64-bit CDF-5 classic file"
        ECODE=1
    fi
fi

rm -f tmp_tst_formatx3 tst_formatx3.nc

exit $ECODE
