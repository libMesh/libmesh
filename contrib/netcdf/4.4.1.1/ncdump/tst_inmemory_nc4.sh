#!/bin/sh
if test "x$SETX" = x1 ; then echo "file=$0"; set -x ; fi
verbose=1
set -e
if test "x$builddir" = "x"; then builddir=`pwd`; fi
if test "x$srcdir" = "x"; then srcdir=`dirname $0`; fi

# Make buildir absolute
cd $builddir
builddir=`pwd`

# Make srcdir be absolute
cd $srcdir
srcdir=`pwd`
cd $builddir

# Setup
PASS=1

# Define the .cdl files to test
CLASSIC="small ref_tst_nans ref_tst_utf8"
EXTENDED="ref_nc_test_netcdf4 ref_tst_comp ref_tst_opaque_data"

rm -fr ./results
mkdir ./results

# Dump classic files two ways and compare
dotest() {
K=$1
for f in $2 ; do
  echo "Testing ${f}"
  ${builddir}/../ncgen/ncgen -$K -o ./results/${f}.nc ${srcdir}/${f}.cdl
  ./ncdump ./results/${f}.nc > ./results/${f}.cdl
  ./ncdump -Xm ./results/${f}.nc > ./results/${f}.cdx
  diff -w ./results/${f}.cdl ./results/${f}.cdx &> ./results/${f}.diff
  if test -s ./results/${f}.diff ; then
    echo "***FAIL: $f"
    PASS=0
  fi
done
}

dotest "3" "$CLASSIC"
dotest "5" "$EXTENDED5"

if test -f ${builddir}/../config.h ; then
  if fgrep -e '#define USE_NETCDF4 1' ${builddir}/../config.h >/dev/null ; then
    dotest "4" "$EXTENDED4"
  fi
fi

# Cleanup
rm -fr results

if test "x$PASS" = x1 ; then
  echo "*** PASS all tests"
  CODE=0
else
  CODE=1
fi
exit $CODE
