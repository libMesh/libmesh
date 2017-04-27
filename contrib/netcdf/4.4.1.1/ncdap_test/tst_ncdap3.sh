#!/bin/sh

set -e

#X="-x"
#grind="checkleaks"

# if this is part of a distcheck action, then this script
# will be executed in a different directory
# than the one containing it; so capture the path to this script
# as the location of the source directory.

srcdir=`dirname $0`
#set -x
# compute the build directory
# Do a hack to remove e.g. c: for MSYS
cd `pwd`
builddir=`pwd`/..

# Hack for MSYS
cd $srcdir
srcdir=`pwd`
if [ `uname | cut -d "_" -f 1` = "MINGW32" ]; then
    srcdir=`pwd | sed 's/\/c\//c:\//g'`
    builddir=`echo $builddir | sed 's/\/c\//c:\//g'`

    srcdir=`pwd | sed 's/\/g\//g:\//g'`
    builddir=`echo $builddir | sed 's/\/g\//g:\//g'`
fi



cd ${builddir}/ncdap_test


#exec sh $X ${srcdir}/tst_ncdap.sh "$srcdir" "$builddir" "file3" $grind
exec sh $X ${srcdir}/tst_ncdap.sh "$srcdir" "$builddir" "dds3" $grind

