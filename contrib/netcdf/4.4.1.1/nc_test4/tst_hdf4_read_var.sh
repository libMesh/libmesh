#!/bin/bash
# This shell script tests that an hdf4 file can be read a
# variable at a time.
#
# this was added in support of https://github.com/Unidata/netcdf-c/issues/264

FILE=tst_interops2.h4

set -e

echo ""
echo "*** Testing reading an individual variable from an HDF4 file."

../ncdump/ncdump -v hdf4_dataset_type_0 $FILE
../ncdump/ncdump -v hdf4_dataset_type_1 $FILE
../ncdump/ncdump -v hdf4_dataset_type_2 $FILE
../ncdump/ncdump -v hdf4_dataset_type_3 $FILE
../ncdump/ncdump -v hdf4_dataset_type_4 $FILE
../ncdump/ncdump -v hdf4_dataset_type_5 $FILE
../ncdump/ncdump -v hdf4_dataset_type_6 $FILE
../ncdump/ncdump -v hdf4_dataset_type_7 $FILE

echo "*** Success."
