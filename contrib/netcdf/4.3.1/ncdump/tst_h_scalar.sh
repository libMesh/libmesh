#!/bin/sh
# This shell script runs ncdump to verify scalar attribute and variable output

set -e
echo ""
echo "*** Running ncdump scalar test."

# echo "*** dumping tst_h_scalar.nc to tst_h_scalar.cdl..."
./ncdump tst_h_scalar.nc > tst_h_scalar.cdl
# echo "*** comparing tst_h_scalar.cdl with ref_tst_h_scalar.cdl..."
diff -b tst_h_scalar.cdl $srcdir/cdl4/ref_tst_h_scalar.cdl

echo "*** All ncdump scalar test output for netCDF-4 format passed!"
exit 0

