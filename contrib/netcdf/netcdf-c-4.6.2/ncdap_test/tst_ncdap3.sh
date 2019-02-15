#!/bin/sh

if test "x$SETX" = x1 ; then set -x ; fi

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh
set -e

. ${srcdir}/tst_utils.sh

# get the list of test files
. ${srcdir}/tst_filelists.sh

# Test executor
dotests() {
for x in ${FILETESTS} ; do
  url="${PARAMS}${FILEURL}/$x"
  if test "x$quiet" = "x0" ; then echo "*** Testing: ${x} ; url=$url" ; fi
  # determine if this is an xfailtest
  isxfail=0
  if test "x${XFAILTESTS}" != x ; then
    if IGNORE=`echo -n " ${XFAILTESTS} " | fgrep " ${x} "`; then isxfail=1; fi
  fi
  ok=1
  if ${NCDUMP} ${DUMPFLAGS} "${url}" | sed 's/\\r//g' > ${x}.dmp ; then ok=$ok; else ok=0; fi
  # compare with expected
  if diff -w ${EXPECTED}/${x}.dmp ${x}.dmp  ; then ok=$ok; else ok=0; fi
   processstatus
done
}

TITLE="DAP to netCDF-3 translation using files"
EXPECTED="$expected3"
RESULTSDIR="file_results"

rm -fr ${RESULTSDIR}
mkdir "${RESULTSDIR}"

echo "*** Testing $TITLE "
echo "        Base URL: ${TESTURL}"
echo "        Client Parameters: ${PARAMS}"

cd ${RESULTSDIR}
dotests file
cd ..
summarize
cleanup
doexit
