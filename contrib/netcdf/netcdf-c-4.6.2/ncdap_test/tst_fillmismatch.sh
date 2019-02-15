#!/bin/sh

#export NCPATHDEBUG=1

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

F="fillmismatch.nc"
EXPECTED="${srcdir}/expected3"

URL='file://'
URL="${URL}${srcdir}/testdata3/$F"

# First check that without [fillmismatch], we get a failure
rm -f ./tmp_tst_mismatch
if ${NCDUMP} "${URL}" > ./tmp_tst_mismatch 2>&1 ; then
echo "*** Fail: ${NCDUMP} ${URL} passed"
exit 1
else
echo "*** XFail: ${NCDUMP} ${URL} failed"
fi

# Now check that with [fillmismatch], we get sucess
URL="[fillmismatch]${URL}"
rm -f ./tmp_tst_mismatch
if ${NCDUMP} "${URL}" > ./tmp_tst_mismatch ; then
echo "*** Pass: ${NCDUMP} ${URL} passed"
else
echo "*** Fail: ${NCDUMP} ${URL} failed"
exit 1
fi

# Verify result
diff -w ${EXPECTED}/$F.dmp ./tmp_tst_mismatch
#cleanup
rm -f ./tmp_tst_mismatch
exit
