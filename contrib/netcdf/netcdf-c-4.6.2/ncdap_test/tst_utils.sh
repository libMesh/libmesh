# This is not a standalone script; it is invoked by
# one or more of the main scripts

if test "x$SETX" = x1 ; then set -x ; fi

quiet=0

PARAMS="[log]"
#PARAMS="${PARAMS}[fetch=memory]"
#PARAMS="${PARAMS}[show=fetch]"

OCLOGFILE=/dev/null

DUMPFLAGS=

# Locate directories
testdata3="${srcdir}/testdata3"
expected3="${srcdir}/expected3"

rm -f ./.dodsrc ./.ocrc ./.daprc

passcount=0
xfailcount=0
failcount=0

# Try to figure out our platform
myplatform=`uname -a | cut -d" " -f 1`
case "$myplatform" in
Darwin*) platform=osx ;;
MINGW*) platform=mingw ;;
CYGWIN*) platform=cygwin ;;
linux*) platform=linux ;;
Linux*) platform=linux ;;
*) platform=unknown ;;
esac

# How to access local files
FILEURL="file://${testdata3}"

processstatus() {
  if test "$ok" = 1 ; then
    status=0  # succeed
  elif test "x$isxfail" = "x0" ; then
    status=1  # fail
  else
    status=2  # xfail
  fi
  case "$status" in
  0)
    passcount=`expr $passcount + 1`
    if test "x$quiet" = "x" ; then echo "*** SUCCEED: ${x}"; fi
    ;;
  1)
    failcount=`expr $failcount + 1`
    echo "*** FAIL:  ${x}"
    ;;
  2)
    xfailcount=`expr $xfailcount + 1`
    echo "*** XFAIL : ${x}"
    ;;
  esac
}

summarize() {
totalcount=`expr $passcount + $failcount + $xfailcount`
okcount=`expr $passcount + $xfailcount`
if test "$failcount" -gt 0 ; then
echo "*** FAILED: ${okcount}/${totalcount} ; ${xfailcount} expected failures ; ${failcount} unexpected failures"
else
echo "*** PASSED: ${okcount}/${totalcount} ; ${xfailcount} expected failures ; ${failcount} unexpected failures"
fi
}

cleanup() {
rm -fr ${RESULTSDIR}
}

doexit() {
if test "$failcount" -gt 0
then
  exit 1
else
  exit 0
fi
}

