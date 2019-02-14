#!/bin/bash

if test "x$SETX" = x1 ; then set -x ; fi

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e

. ${srcdir}/tst_utils.sh

# get the list of test url targets
. ${srcdir}/tst_urls.sh

PARAMS="${PARAMS}[netcdf3]"

# Choose tests to run
if test "x$longtests" != x; then
WHICHTESTS="L1 LC1 LC2"
else # Standard test set
WHICHTESTS="S1 C1"
fi

TITLE="DAP to netCDF-3 translation using remote server"
EXPECTED="$expected3"
RESULTSDIR="remote_results"

rm -fr ${RESULTSDIR}
mkdir "${RESULTSDIR}"

echo "*** Testing $TITLE "
echo "        Base URL: ${DTS}"
echo "        Client Parameters: ${PARAMS}"
if test "$cache" = 0; then echo "        Caching: off"; else echo "        Caching: on"; fi
echo "    Note: The remote tests may be slow or even fail if the server is overloaded"

cd ${RESULTSDIR}
for i in $WHICHTESTS ; do
  computewhich $i
  doremotetests remote
done
cd ..
summarize
cleanup
doexit

