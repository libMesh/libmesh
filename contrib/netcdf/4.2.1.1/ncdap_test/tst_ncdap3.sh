#!/bin/sh

#X="-x"
#grind="checkleaks"

# if this is part of a distcheck action, then this script
# will be executed in a different directory
# than the one containing it; so capture the path to this script
# as the location of the source directory.
srcdir=`dirname $0`

# compute the build directory
# Do a hack to remove e.g. c: for CYGWIN
cd `pwd`
builddir=`pwd`/..

# Hack for CYGWIN
cd $srcdir
srcdir=`pwd`
if [ `uname | cut -d "_" -f 1` = "MINGW32" ]; then
    srcdir=`pwd | sed 's/\/c\//c:\//g'`
    builddir=`echo $builddir | sed 's/\/c\//c:\//g'`
fi


cd ${builddir}/ncdap_test


exec sh $X ${srcdir}/tst_ncdap.sh "$srcdir" "$builddir" "file3" $grind

