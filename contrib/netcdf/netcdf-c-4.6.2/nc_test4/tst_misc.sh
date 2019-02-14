#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e

echo "*** Testing phony dimension creation on pure h5 file"
rm -f ./tmp
if $NCDUMP -L0 -K ${srcdir}/tdset.h5 >./tmp ; then
echo "*** Pass: phony dimension creation"
ECODE=0
else
echo "*** Fail: phony dimension creation"
ECODE=1
fi

rm -f tmp

exit $ECODE


