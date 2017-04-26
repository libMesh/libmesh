#!/bin/sh
if test "x$SETX" = x1 ; then echo "file=$0"; set -x ; fi
set -e
echo ""

if test "x$builddir" = "x"; then builddir=`pwd`; fi
if test "x$srcdir" = "x"; then srcdir=`dirname $0`; fi

# Make buildir absolute
cd $builddir
builddir=`pwd`

# Make srcdir be absolute
cd $srcdir
srcdir=`pwd`

cd $builddir

export verbose
export srcdir
export builddir

EXIT=0

NCF=$srcdir/nc4_fileinfo.nc
HDF=$srcdir/hdf5_fileinfo.hdf
NF=ref_tst_compounds4.nc

# Do a false negative test 
rm -f ./tmp
if $builddir/ncdump -s $builddir/$NF | fgrep '_IsNetcdf4 = 0' > ./tmp ; then
   echo "Pass: False negative for file: $NF"
else
   echo "FAIL: False negative for file: $NF"
fi
rm -f ./tmp

if ./tst_fileinfo > /dev/null ; then
   # look at the _IsNetcdf4 flag
   N_IS=`${builddir}/ncdump -s $NCF | fgrep '_IsNetcdf4' | tr -d ' ;'`
   N_IS=`echo $N_IS | cut -d= -f2`
   H_IS=`${builddir}/ncdump -s $HDF | fgrep '_IsNetcdf4' | tr -d ' ;'`
   H_IS=`echo $H_IS | cut -d= -f2`
   if test "x$N_IS" = 'x0' ;then
     echo "FAIL: $NCF is marked as not netcdf-4"
   fi
   if test "x$H_IS" = 'x1' ;then
     echo "FAIL: $HDF is marked as netcdf-4"
   fi
else
echo "FAIL: tst_fileinfo"
EXIT=1
fi

rm -f $NCF
rm -f $HDF

exit $EXIT

