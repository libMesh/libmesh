#!/bin/sh

echo "*** Test Maximum dimension sizes X mode"

set -x

if test "x$SETX" = x1 ; then echo "file=$0"; set -x ; fi
# This shell script tests max dimension sizes X mode

RETURN=0

if test "x$srcdir" = "x"; then
    srcdir=`dirname $0`; 
fi
# add hack for sunos
export srcdir;

echo ""

rm -f tst_dimsize_classic.nc tst_dimsize_64offset.nc tst_dimsize_64data.nc

echo "*** Generate: tst_dimsize_classic.nc tst_dimsize_64offset.nc tst_dimsize_64data.nc"
./tst_dimsizes

echo "*** Verify that ncdump can read dimsizes"

rm -fr ./tmp
if ../ncdump/ncdump -h tst_dimsize_classic.nc > ./tmp ; then
echo "*** PASS: ncdump tst_dimsize_classic.nc"
else
echo "*** FAIL: ncdump tst_dimsize_classic.nc"
RETURN=1
fi

rm -fr ./tmp
if ../ncdump/ncdump -h tst_dimsize_64offset.nc > ./tmp ; then
echo "*** PASS: ncdump tst_dimsize_64offset.nc"
else
echo "*** FAIL: ncdump tst_dimsize_64offset.nc"
RETURN=1
fi

if test -f tst_dimsize_64data.nc ; then
  rm -fr ./tmp
  if ../ncdump/ncdump -h tst_dimsize_64data.nc > ./tmp ; then
    echo "*** PASS: ncdump tst_dimsize_64data.nc"
  else
    echo "*** FAIL: ncdump tst_dimsize_64data.nc"
  RETURN=1
  fi
fi

# Cleanup
rm -f tmp tst_dimsize_classic.nc tst_dimsize_64offset.nc tst_dimsize_64data.nc

exit $RETURN

