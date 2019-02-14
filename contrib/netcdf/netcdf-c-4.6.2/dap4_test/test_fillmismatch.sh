#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e
. ${srcdir}/d4test_common.sh

echo "test_fillmismatch.sh:"

F="test_fillmismatch.nc"

URL='[dap4]file://'
URL="${URL}${srcdir}/misctestfiles/$F"

# First check that without [fillmismatch], we get a failure
rm -f ./tmp_dap4_mismatch
if ${NCDUMP} -h "${URL}" > ./tmp_dap4_mismatch 2>&1 ; then
echo "*** Fail: ${NCDUMP} ${URL} passed"
exit 1
else
echo "*** XFail: ${NCDUMP} ${URL} failed"
fi

# Now check that with [fillmismatch], we get sucess
URL="[fillmismatch]${URL}"
rm -f ./tmp_dap4_mismatch
if ${NCDUMP} -h "${URL}" > ./tmp_dap4_mismatch ; then
echo "*** Pass: ${NCDUMP} ${URL} passed"
else
echo "*** Fail: ${NCDUMP} ${URL} failed"
exit 1
fi

# Verify result
diff -w ${srcdir}/baselineraw/$F.dmp ./tmp_dap4_mismatch
#cleanup
rm -f ./tmp_dap4_mismatch
exit
