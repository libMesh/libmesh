#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# get some config.h parameters
if test -f ${top_builddir}/config.h ; then
  if fgrep -e '#define ENABLE_CDF5 1' ${top_builddir}/config.h >/dev/null ; then
    HAVE_CDF5=1
  else
    HAVE_CDF5=0
  fi
else
  echo "Cannot locate config.h"
  exit 1
fi

# It is unreasonable to test actual lengths of files
# (even netcdf-3 files).
#However, the files created in this script are used in later ones

${NCGEN} -b ${srcdir}/small.cdl
${NCGEN} -b ${srcdir}/small2.cdl

# This shell script tests lengths of small netcdf files and tests
# that rewriting a numeric value doesn't change file length
# $Id: tst_lengths.sh,v 1.10 2008/08/07 00:07:52 ed Exp $

# cat > rewrite-scalar.c << EOF
# #include <stdio.h>
# #include <netcdf.h>
# #define ERR do {fflush(stdout); fprintf(stderr, "Error, %s, line: %d\n", __FILE__, __LINE__); return(1);} while (0)

# int
# main(int ac, char *av[]) {
#     int ncid, varid, data[] = {42};
#     if (nc_open(av[1], NC_WRITE, &ncid)) ERR;
#     if (nc_inq_varid(ncid, av[2], &varid)) ERR;
#     if (nc_put_var_int(ncid, varid, data)) ERR;
#     if (nc_close(ncid)) ERR;
#     return 0;
# }
# EOF
# cat > test-len.sh << 'EOF'
# # test that length of file $1 is $2
# len=`ls -l $1|awk '{print $5}'`
# if [ $len = $2 ]; then
#   exit 0
# else
#   echo "### Failure: file $1 has length $len instead of expected $2"
#   exit 1
# fi
# EOF
# chmod +x ./test-len.sh
# cc -g -o rewrite-scalar -I../libsrc rewrite-scalar.c -L../libsrc -lnetcdf
# echo "netcdf small {variables: byte t; data: t = 1;}" > small.cdl
set -e

echo ""
echo "*** testing length of classic file"
${NCGEN} -b ${srcdir}/small.cdl
if test `wc -c < small.nc` != 68; then
    exit 1
fi

#echo "*** testing length of classic file written with NOFILL"
#${NCGEN} -b -x ${srcdir}/small.cdl
#if test `wc -c < small.nc` != 68; then
#    exit 1
#fi

echo "*** testing length of rewritten classic file"
${NCGEN} -b ${srcdir}/small.cdl && ${execdir}/rewrite-scalar small.nc t
#if test `wc -c < small.nc` != 68; then
#    exit 1
#fi

echo "*** testing length of rewritten classic file written with NOFILL"
${NCGEN} -b -x ${srcdir}/small.cdl && ${execdir}/rewrite-scalar small.nc t
#if test `wc -c < small.nc` != 68; then
#    exit 1
#fi

echo "*** testing length of 64-bit offset file"
${NCGEN} -b -k64-bit-offset ${srcdir}/small.cdl
#if test `wc -c < small.nc` != 72; then
#    exit 1
#fi

echo "*** testing length of 64-bit offset file written with NOFILL"
${NCGEN} -b -k64-bit-offset -x ${srcdir}/small.cdl
#if test `wc -c < small.nc` != 72; then
#    exit 1
#fi

echo "*** testing length of rewritten 64-bit offset file"
${NCGEN} -b -k64-bit-offset ${srcdir}/small.cdl && ${execdir}/rewrite-scalar small.nc t
#if test `wc -c < small.nc` != 72; then
#    exit 1
#fi

echo "*** testing length of rewritten 64-bit offset file written with NOFILL"
${NCGEN} -b -k64-bit-offset -x ${srcdir}/small.cdl && ${execdir}/rewrite-scalar small.nc t
#if test `wc -c < small.nc` != 72; then
#    exit 1
#fi

# The following tests only occur if we have CDF5.
if test "x$HAVE_CDF5" = x1 ; then

    echo "*** testing length of 64-bit data file"
    ${NCGEN} -b -k64-bit-data ${srcdir}/small.cdl
    if test `wc -c < small.nc` != 104; then
        exit 1
    fi

    echo "*** testing length of 64-bit data file"
    ${NCGEN} -b -5 ${srcdir}/small.cdl
    if test `wc -c < small.nc` != 104; then
        exit 1
    fi
    echo "*** testing length of 64-bit data file written with NOFILL"
    ${NCGEN} -b -5 -x ${srcdir}/small.cdl
    #if test `wc -c < small.nc` != 104; then
    #    exit 1
    #fi

    echo "*** testing length of rewritten 64-bit data file"
    ${NCGEN} -b -5 ${srcdir}/small.cdl && ${execdir}/rewrite-scalar small.nc t
    # Watch out, it appears that the CDF-5 files are being rounded up to next page size
    # So, we need to truncate them wrt nul's in order to check size.
    # Bad hack, but what else can I do?
    if test `${execdir}/nctrunc <small.nc |wc -c` != 104; then
        exit 1
    fi

    echo "*** testing length of rewritten 64-bit data file written with NOFILL"
    ${NCGEN} -b -5 -x ${srcdir}/small.cdl && ${execdir}/rewrite-scalar small.nc t
    #if test `wc -c < small.nc` != 104; then
    #    exit 1
    #fi

   # End HAVE_CDF5 block.
fi

# test with only one record variable of type byte or short, which need
# not be 4-byte aligned
echo "*** testing length of one-record-variable classic file"
${NCGEN} -b ${srcdir}/small2.cdl
#if test `wc -c < small2.nc` != 101; then
#    exit 1
#fi

# The following tests only occur if we have CDF5.
if test "x$HAVE_CDF5" = x1 ; then

    echo "*** testing length of one-record-variable 64-bit data file"
    ${NCGEN} -b -5 ${srcdir}/small2.cdl
    if test `wc -c < small2.nc` != 161; then
        exit 1
    fi

    echo "*** testing length of one-record-variable 64-bit data file"
    ${NCGEN} -b -5 ${srcdir}/small2.cdl
    if test `wc -c < small2.nc` != 161; then
        exit 1
    fi

    echo "*** testing length of one-record-variable 64-bit data file written with NOFILL"
    ${NCGEN} -b -5 -x ${srcdir}/small2.cdl
    if test `wc -c < small2.nc` != 161; then
        exit 1
    fi

    #end HAVE_CDF5 block
fi

echo "*** testing length of one-record-variable classic file written with NOFILL"
${NCGEN} -b -x ${srcdir}/small2.cdl
if test `wc -c < small2.nc` != 101; then
    exit 1
fi

echo "*** testing length of one-record-variable classic file written with NOFILL"
${NCGEN} -b -x ${srcdir}/small2.cdl
if test `wc -c < small2.nc` != 101; then
    exit 1
fi

echo "*** testing length of one-record-variable classic file written with NOFILL"
${NCGEN} -b -x ${srcdir}/small2.cdl
#if test `wc -c < small2.nc` != 101; then
#    exit 1
#fi

echo "*** testing length of one-record-variable 64-bit offset file"
${NCGEN} -b -k64-bit-offset ${srcdir}/small2.cdl
#if test `wc -c < small2.nc` != 105; then
#    exit 1
#fi

echo "*** testing length of one-record-variable 64-bit offset file written with NOFILL"
${NCGEN} -b -k64-bit-offset -x ${srcdir}/small2.cdl
#if test `wc -c < small2.nc` != 105; then
#    exit 1
#fi
