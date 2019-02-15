#!/bin/sh

# This shell gets some files from the netCDF ftp site for testing,
# then runs the tst_knmi benchmarking program.
# Ed Hartnett

# Load common values for netCDF shell script tests.
if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# Get files if needed.
echo ""
file_list="MSGCPP_CWP_NC3.nc  MSGCPP_CWP_NC4.nc"
echo "Getting KNMI test files $file_list"
for f1 in $file_list
do
    if ! test -f $f1; then
	wget ftp://ftp.unidata.ucar.edu/pub/netcdf/sample_data/$f1.gz
	gunzip $f1.gz
    fi
done

# Run the C program on these files.
${execdir}/tst_knmi

echo "SUCCESS!!!"

exit 0
