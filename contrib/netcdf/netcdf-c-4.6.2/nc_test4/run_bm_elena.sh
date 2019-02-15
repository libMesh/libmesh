#!/bin/sh

# This shell runs some benchmarks that Elena ran as described here:
# http://hdfeos.org/workshops/ws06/presentations/Pourmal/HDF5_IO_Perf.pdf

# Ed Hartnett

# Load common values for netCDF shell script tests.
if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# Run benchmarks.
echo ""
echo "*** Testing the benchmarking program bm_file for simple float file, no compression..."
${execdir}/bm_file -h -d -f 3 -o  tst_elena_out.nc -c 0:-1:0:1024:16:256 tst_elena_int_3D.nc
${execdir}/bm_file -d -f 3 -o  tst_elena_out.nc -c 0:-1:0:1024:256:256 tst_elena_int_3D.nc
${execdir}/bm_file -d -f 3 -o  tst_elena_out.nc -c 0:-1:0:512:64:256 tst_elena_int_3D.nc
${execdir}/bm_file -d -f 3 -o  tst_elena_out.nc -c 0:-1:0:512:256:256 tst_elena_int_3D.nc
${execdir}/bm_file -d -f 3 -o  tst_elena_out.nc -c 0:-1:0:256:64:256 tst_elena_int_3D.nc
${execdir}/bm_file -d -f 3 -o  tst_elena_out.nc -c 0:-1:0:256:256:256 tst_elena_int_3D.nc
echo '*** SUCCESS!!!'

exit 0
