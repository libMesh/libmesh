#!/bin/sh

#set -e

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

echo "*** Test Maximum dimension sizes X mode"

# This shell script tests max dimension sizes X mode

RETURN=0

echo ""

rm -f tst_dimsize_classic.nc tst_dimsize_64offset.nc tst_dimsize_64data.nc

echo "*** Generate: tst_dimsize_classic.nc tst_dimsize_64offset.nc tst_dimsize_64data.nc"
${execdir}/tst_dimsizes

echo "*** Verify that ncdump can read dimsizes"

rm -fr ./tmp_tst_dimsizes
if ${NCDUMP} -h tst_dimsize_classic.nc > ./tmp_tst_dimsizes ; then
echo "*** PASS: ncdump tst_dimsize_classic.nc"
else
echo "*** FAIL: ncdump tst_dimsize_classic.nc"
RETURN=1
fi

rm -fr ./tmp_tst_dimsizes
if ${NCDUMP} -h tst_dimsize_64offset.nc > ./tmp_tst_dimsizes ; then
echo "*** PASS: ncdump tst_dimsize_64offset.nc"
else
echo "*** FAIL: ncdump tst_dimsize_64offset.nc"
RETURN=1
fi

if test -f tst_dimsize_64data.nc ; then
  rm -fr ./tmp_tst_dimsizes
  if ${NCDUMP} -h tst_dimsize_64data.nc > ./tmp_tst_dimsizes ; then
    echo "*** PASS: ncdump tst_dimsize_64data.nc"
  else
    echo "*** FAIL: ncdump tst_dimsize_64data.nc"
  RETURN=1
  fi
fi

# Cleanup
rm -f tmp_tst_dimsizes tst_dimsize_classic.nc tst_dimsize_64offset.nc tst_dimsize_64data.nc

exit $RETURN
