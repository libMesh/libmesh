#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e

#Constants
CDL5=tst_diskless5.cdl
FILE5=tst_diskless5.nc

echo ""
rm -f $FILE5
# Generate FILE5
${NCGEN} -3 -lb -o ${FILE5} ${srcdir}/${CDL5}

echo ""
${execdir}/tst_diskless5

# cleanup
rm -f $FILE5

exit
