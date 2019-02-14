#!/bin/sh
#
# Added in support of https://github.com/Unidata/netcdf-c/gh425 and
# https://github.com/Unidata/netcdf-c/gh469
#
# The output isn't validated, but the regression it is fixing fails on nccopy.
#

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# For a netCDF-3 build, test nccopy on netCDF files in this directory

set -e
echo ""

INFILE=${TOPSRCDIR}/ncdump/ref_nccopy3_subset.nc
OUTFILE=nccopy3_subset_out.nc

echo "*** Testing netCDF-3 nccopy -v/-V flags on $IN"

echo "*** One Dimensional Tests"
echo "*** Testing nccopy -v"

${NCCOPY} -v lat ${INFILE} ${OUTFILE}

echo "*** Testing nccopy -V"

${NCCOPY} -V lat ${INFILE} ${OUTFILE}

echo "*** Two Dimensional Tests"
echo "*** Testing nccopy -v"

${NCCOPY} -v lat_2D_rct ${INFILE} ${OUTFILE}

echo "*** Testing nccopy -V"

${NCCOPY} -V lat_2D_rct ${INFILE} ${OUTFILE}

echo "nccopy passed!"
exit 0
