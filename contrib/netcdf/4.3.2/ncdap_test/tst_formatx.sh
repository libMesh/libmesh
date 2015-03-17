#!/bin/sh
# This shell script tests the output several previous tests.
# $Id: tst_output.sh,v 1.17 2010/05/14 16:21:15 ed Exp $


echo ""
echo "*** Testing extended file format output."
set -e

# Figure our dst server
SVC=`./nctestserver dts`
if test "x$SVC" = "x" ; then
echo "cannot locate test server"
exit
fi
URL="$SVC/dts/test.03"

ECODE=0
echo "Test extended format output for a DAP2  file"
rm -f tmp
../ncdump/ncdump -K $URL >tmp
if ! fgrep 'DAP2 mode=00000000' <tmp ; then
echo "*** Fail: extended format for a DAP2 file"
ECODE=1
fi

rm tmp

exit $ECODE
