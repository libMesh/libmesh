#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

verbose=1
set -e

# Setup
PASS=1

# Define the .cdl files to test
CLASSIC="small ref_tst_nans ref_tst_utf8"
EXTENDED="ref_nc_test_netcdf4 ref_tst_comp ref_tst_opaque_data"

rm -fr ./results_tst_inmemory_nc3
mkdir ./results_tst_inmemory_nc3

# Dump classic files two ways and compare
dotest() {
K=$1
for f in $2 ; do
  echo "Testing ${f}"
  ${NCGEN} -$K -o ./results_tst_inmemory_nc3/${f}.nc ${srcdir}/${f}.cdl
  ${NCDUMP} ./results_tst_inmemory_nc3/${f}.nc > ./results_tst_inmemory_nc3/${f}.cdl
  ${NCDUMP} -Xm ./results_tst_inmemory_nc3/${f}.nc > ./results_tst_inmemory_nc3/${f}.cdx
  diff -w ./results_tst_inmemory_nc3/${f}.cdl ./results_tst_inmemory_nc3/${f}.cdx &> ./results_tst_inmemory_nc3/${f}.diff
  if test -s ./results_tst_inmemory_nc3/${f}.diff ; then
    echo "***FAIL: $f"
    PASS=0
  fi
done
}

dotest "3" "$CLASSIC"

# Cleanup
rm -fr results_tst_inmemory_nc3

if test "x$PASS" = x1 ; then
  echo "*** PASS all tests"
  CODE=0
else
  CODE=1
fi
exit $CODE
