#!/bin/sh

# This shell script runs tst_interops2 to create a HDF4 file, and read
# it with netCDF. Then the script runs ncdump on the HDF4 file.

# Ed Hartnett

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

FILE=tst_interops2.h4

ECODE=0

echo ""
echo "*** Testing extended file format output."
set -e

echo "Creating HDF4 file"
${execdir}/tst_interops2

echo "Test extended format output for a HDF4 file"
rm -f tmp_tst_formatx_hdf4
${NCDUMP} -K $FILE >tmp_tst_formatx_hdf4
if ! fgrep 'HDF4 mode=00001000' <tmp_tst_formatx_hdf4 ; then
TMP=`cat tmp_tst_formatx_hdf4`
echo "*** Fail: extended format for an HDF4 file: result=" $TMP
ECODE=1
fi

rm -f tmp_tst_formatx_hdf4

# Exit if there was a failure.
if test $ECODE = 1 ; then
    exit $ECODE
fi

echo ""
echo "*** Testing reading an individual variable from an HDF4 file."

${NCDUMP} -v hdf4_dataset_type_0 $FILE
${NCDUMP} -v hdf4_dataset_type_1 $FILE
${NCDUMP} -v hdf4_dataset_type_2 $FILE
${NCDUMP} -v hdf4_dataset_type_3 $FILE
${NCDUMP} -v hdf4_dataset_type_4 $FILE
${NCDUMP} -v hdf4_dataset_type_5 $FILE
${NCDUMP} -v hdf4_dataset_type_6 $FILE
${NCDUMP} -v hdf4_dataset_type_7 $FILE

echo "*** Success."

