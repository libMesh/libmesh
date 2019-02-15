#!/bin/sh
# This shell script tests the output several previous tests.
# $Id: tst_output.sh,v 1.17 2010/05/14 16:21:15 ed Exp $

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

echo ""
echo "*** Testing extended file format output."
set -e

# Figure our dst server
DTS=`./findtestserver dap2 dts`
if test "x$DTS" = "x" ; then
echo "cannot locate test server for dts"
exit
fi
URL="$DTS/test.03"

ECODE=0
echo "Test extended format output for a DAP2  file"
rm -f tmp_tst_formatx
${NCDUMP} -K "${URL}" >tmp_tst_formatx
if ! grep 'DAP2 mode=00000000' <tmp_tst_formatx ; then
echo "*** Fail: extended format for a DAP2 file"
ECODE=1
fi

rm tmp_tst_formatx

exit $ECODE
