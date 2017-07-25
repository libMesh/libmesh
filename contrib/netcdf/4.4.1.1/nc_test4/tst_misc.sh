#!/bin/sh

if test "x$builddir" = "x"; then builddir=`pwd`; fi
if test "x$srcdir" = "x"; then srcdir=`dirname $0`; fi

# Make buildir absolute
cd $builddir
builddir=`pwd`

set -e

echo "*** Testing phony dimension creation on pure h5 file"
rm -f ./tmp
if ../ncdump/ncdump -L0 -K ${srcdir}/tdset.h5 >./tmp ; then
echo "*** Pass: phony dimension creation"
ECODE=0
else
echo "*** Fail: phony dimension creation"
ECODE=1
fi

rm -f tmp

exit $ECODE


