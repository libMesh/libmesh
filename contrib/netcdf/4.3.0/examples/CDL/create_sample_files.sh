#!/bin/sh
# This shell script creates the netCDF example files from CDL scipts.
# $Id: create_sample_files.sh,v 1.2 2006/07/14 18:39:45 ed Exp $

set -e
echo ""
echo "*** Creating example data files from CDL scripts."
echo "*** creating simple_xy.nc..."
../../ncgen/ncgen -b -o simple_xy.nc $srcdir/simple_xy.cdl

echo "*** checking sfc_pres_temp.nc..."
../../ncgen/ncgen -b -o sfc_pres_temp.nc $srcdir/sfc_pres_temp.cdl

echo "*** checking pres_temp_4D.nc..."
../../ncgen/ncgen -b -o pres_temp_4D.nc $srcdir/pres_temp_4D.cdl

echo "*** All example creations worked!"
exit 0
