#!/bin/sh

#PROF=1
#DEBUG=1
#MEM=1

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -x
set -e

echo "Testing performance of nc_create and nc_open on file with large metadata"

ARGS="--treedepth=6 \
--ngroups=2 \
--ngroupattrs=100 \
--ndims=100 \
--ntypes=10 \
--nvars=100 \
--varrank=2 \
--nvarattrs=500"

echo "timing bigmeta:"
${execdir}/bigmeta $ARGS
echo "timing openbigmeta:"
${execdir}/openbigmeta
if test "x$PROF" = x1 ; then
rm -f perftest.txt
gprof openbigmeta gmon.out >perftest.txt
fi

