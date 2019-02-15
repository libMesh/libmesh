#!/bin/sh
# This shell script runs the cmp test on the example programs.
# $Id: do_comps.sh,v 1.1 2006/06/27 17:44:54 ed Exp $

set -e

##
# This stanza was originally in create_sample_files.sh.
# Moved here.
##
echo ""
echo "*** Creating example data files from CDL scripts."
echo "*** creating simple_xy.nc..."
../../ncgen/ncgen -b -o simple_xy.nc $srcdir/simple_xy.cdl

echo "*** checking sfc_pres_temp.nc..."
../../ncgen/ncgen -b -o sfc_pres_temp.nc $srcdir/sfc_pres_temp.cdl

echo "*** checking pres_temp_4D.nc..."
../../ncgen/ncgen -b -o pres_temp_4D.nc $srcdir/pres_temp_4D.cdl

echo "*** All example creations worked!"

##
# End create_sample_files.
##

echo ""
echo "*** Testing that the CDL examples produced same files as C examples."
echo "*** checking simple_xy.nc..."
cmp simple_xy.nc ../C/simple_xy.nc

echo "*** checking sfc_pres_temp.nc..."
cmp sfc_pres_temp.nc ../C/sfc_pres_temp.nc

echo "*** checking pres_temp_4D.nc..."
cmp pres_temp_4D.nc ../C/pres_temp_4D.nc

echo "*** All CDL example comparisons worked!"
exit 0
