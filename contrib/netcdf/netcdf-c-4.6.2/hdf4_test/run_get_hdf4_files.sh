#!/bin/sh

#Q=-q

# This shell gets some sample HDF4 files from the netCDF ftp site for
# testing. Then it runs program tst_interops3 on the test file to
# check that HDF4 reading works.

# Ed Hartnett

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# Get a file from the ftp site; retry several times
getfile() {
   FTPFILE="ftp://ftp.unidata.ucar.edu/pub/netcdf/sample_data/hdf4/$1.gz"

   for try in 1 2 3 4 ; do # try 4 times

     # signal success/failure
     if wget -c $Q --passive-ftp $FTPFILE ; then
       return 0 # got it
     fi
     echo "wget failed: try $try"
     sleep 5 # seconds
   done
   return 1 # did not get it
}

set -e
echo ""
echo "Getting HDF4 sample files from Unidata FTP site..."

file_list="AMSR_E_L2_Rain_V10_200905312326_A.hdf AMSR_E_L3_DailyLand_V06_20020619.hdf \
    MYD29.A2009152.0000.005.2009153124331.hdf MYD29.A2002185.0000.005.2007160150627.hdf \
    MOD29.A2000055.0005.005.2006267200024.hdf"
echo "Getting HDF4 test files $file_list"

for f1 in $file_list
do
    if ! test -f $f1; then
	if getfile $f1 ; then
  	  gunzip $f1.gz
	else
          echo Could not ftp $f1.gz
          return 1
	fi
    fi
done


echo ""
echo "Running test program to check HDF4 sample files..."
${execdir}/tst_interops3

echo "SUCCESS!!!"

exit 0
