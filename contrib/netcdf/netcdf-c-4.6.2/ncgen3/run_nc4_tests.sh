#!/bin/sh
# This shell script runs the ncdump tests.
# $Id: run_nc4_tests.sh,v 1.3 2009/09/24 18:19:10 dmh Exp $

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh


echo "*** Testing ncgen3 for netCDF-4."
set -e
echo "*** creating netCDF-4 file c0_4.nc from c0.cdl..."
${NCGEN3} -k3 -b -o c0_4.nc ${ncgen3c0}
echo "*** creating netCDF-4 classic model file c0_4c.nc from c0.cdl..."
${NCGEN3} -k4 -b -o c0_4c.nc ${ncgen3c0}

echo "*** Test successful!"
exit 0
