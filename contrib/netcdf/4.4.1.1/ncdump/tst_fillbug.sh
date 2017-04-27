#!/bin/sh
if test "x$SETX" = x1 ; then echo "file=$0"; set -x ; fi
# This shell script runs an ncdump bug test for netcdf-4
# $Id: tst_fillbug.sh,v 1.1 2008/10/02 19:49:52 russ Exp $

set -e
if test "x$srcdir" = "x"; then
    srcdir=`dirname $0`;
fi
export srcdir



echo ""
echo "*** Running ncdump bug test."

# echo "*** dumping tst_fillbug.nc to tst_fillbug.cdl..."
./ncdump tst_fillbug.nc > tst_fillbug.cdl
# echo "*** comparing tst_fillbug.cdl with ref_tst_fillbug.cdl..."
diff -b tst_fillbug.cdl $srcdir/ref_tst_fillbug.cdl

echo "*** All ncdump bug test output for netCDF-4 format passed!"
exit 0
