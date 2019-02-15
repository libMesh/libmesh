#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# This shell script tests a file which doesn't adhere to the
# documented null-byte padding for header information.
# See https://github.com/Unidata/netcdf-c/issues/657 for more info.

set -e

echo ""
echo "*** Testing compatibility with non-null-byte padded test file."
${NCDUMP} $srcdir/ref_null_byte_padding_test.nc

echo "Passed."
exit 0
