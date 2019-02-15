#!/bin/sh

#export NCPATHDEBUG=1

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

F="ber-2002-10-01.nc"
EXPECTED="${srcdir}/expected3"

URL='[log][cache]file://'
URL="${URL}${srcdir}/testdata3/$F"
rm -f ./tmp_tst_ber
${NCDUMP} "${URL}" | sed 's/\\r//g' > ./tmp_tst_ber
diff -w ${EXPECTED}/$F.dmp ./tmp_tst_ber
#cleanup
rm -f ./tmp_tst_ber
exit
