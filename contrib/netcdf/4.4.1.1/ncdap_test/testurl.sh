#!/bin/sh
#set -x

#NOP=1
#NOS=1
#NOB=1

#SHOW=1
#DBG=1
#GDB=1

# if this is part of a distcheck action, then this script
# will be executed in a different directory
# than the one containing it; so capture the path to this script
# as the location of the source directory.

# capture the build directory
# Do a hack to remove e.g. c: for CYGWIN
builddir=`pwd`/..
if test "x$TOPSRCDIR" != x ; then
srcdir="$TOPSRCDIR/ncdap_test"
else
srcdir=`dirname $0`
fi
# canonical
cd $srcdir
srcdir=`pwd`

# Hack for CYGWIN
if [ `uname | cut -d "_" -f 1` = "MINGW32" ]; then
    srcdir=`echo $srcdir | sed 's/\/c\//c:\//g'`
    builddir=`echo $builddir | sed 's/\/c\//c:\//g'`
fi
cd ${builddir}/ncdap_test

OCLOGFILE=stderr
if test "x$DBG" = x1 ; then
SHOW=1
fi

NCDUMP=$builddir/ncdump/ncdump

URL="http://remotetest.unidata.ucar.edu/dts/test.03"

PREFIX="[log][show=fetch]"
SUFFIX="log&show=fetch"
BOTHP="[log][show=fetch]"
BOTHS="noprefetch&fetch=disk"

locreset () {
    rm -f ./tmp ./errtmp
}

buildurl () {
  front="$1"
  back="$2"
  url="${front}${URL}"
  if test "x$back" != x ; then
    url="${url}#${back}"
  fi
}

pass=1

if test "x$GDB" = x1 ; then
NCDUMP="gdb --args $NCDUMP"
fi

# Initialize
locreset

if test "x$NOP" != x1 ; then
echo "***Testing url prefix parameters"
buildurl $PREFIX ""
# Invoke ncdump to extract the URL
echo "command: ${NCDUMP} -h $url"
${NCDUMP} -h "$url" >./tmp 2> ./errtmp
if test "x${SHOW}" = x1 ; then cat ./tmp ; fi
fi

locreset
if test "x$NOS" != x1 ; then
echo "***Testing url suffix parameters"
buildurl "" $SUFFIX
# Invoke ncdump to extract the URL
echo "command: ${NCDUMP} -h $url"
${NCDUMP} -h "$url" >./tmp  2> ./errtmp
if test "x${SHOW}" = x1 ; then cat ./tmp ; fi
fi

locreset
if test "x$NOB" != x1 ; then
echo "***Testing url prefix+suffix parameters"
buildurl $BOTHP $BOTHS
# Invoke ncdump to extract the URL
echo "command: ${NCDUMP} -h $url"
${NCDUMP} -h "$url" >./tmp 2> ./errtmp
if test "x${SHOW}" = x1 ; then cat ./tmp ; fi
fi

locreset

if test "x$pass" = x0 ; then
  echo "***FAIL"
  exit 1
fi
echo "***PASS"
exit 0


