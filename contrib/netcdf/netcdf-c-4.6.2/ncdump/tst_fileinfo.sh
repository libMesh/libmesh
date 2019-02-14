#!/bin/sh
# Author: Dennis Heimbigner

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e
echo ""

EXIT=0

NCF="./nc4_fileinfo.nc"
HDF="./hdf5_fileinfo.hdf"

NF="${top_srcdir}/ncdump/ref_tst_compounds4.nc"
NPV1="${top_srcdir}/ncdump/ref_provenance_v1.nc"

# Create various files
${execdir}/tst_fileinfo

# Do a false negative test
rm -f ./tst_fileinfo.tmp
if $NCDUMP -s $NF | fgrep '_IsNetcdf4 = 0' > ./tst_fileinfo.tmp ; then
   echo "Pass: False negative for file: $NF"
else
   echo "FAIL: False negative for file: $NF"
   EXIT=1
fi

rm -f ./tst_fileinfo.tmp
if test -e $NCF ; then
   # look at the _IsNetcdf4 flag
   N_IS=`${NCDUMP} -s $NCF | fgrep '_IsNetcdf4' | tr -d ' ;'`
   N_IS=`echo $N_IS | cut -d= -f2`
   H_IS=`${NCDUMP} -s $HDF | fgrep '_IsNetcdf4' | tr -d ' ;'`
   H_IS=`echo $H_IS | cut -d= -f2`
   if test "x$N_IS" = 'x0' ;then
     echo "FAIL: $NCF is marked as not netcdf-4"
     EXIT=1
   fi
   if test "x$H_IS" = 'x1' ;then
     echo "FAIL: $HDF is marked as netcdf-4"
     EXIT=1
   fi
else
  echo "FAIL: tst_fileinfo: $NCF does not exist"
  EXIT=1
fi

# Test what happens when we read a file that used provenance version 1
rm -f ./tst_fileinfo.tmp ./tst_fileinfo2.tmp
$NCDUMP -hs $NPV1 >tst_fileinfo2.tmp
fgrep '_NCProperties' <tst_fileinfo2.tmp > ./tst_fileinfo.tmp
if ! fgrep 'version=1' tst_fileinfo.tmp ; then
  echo "FAIL: $NPV1 is not marked as version=1"
  EXIT=1
fi

rm -f $NCF
rm -f $HDF
rm -f *.tmp

if test "x$EXIT" = x0 ; then
echo "*** Pass all tests"
fi
exit $EXIT
