#!/bin/sh
# This shell script runs the cmp test on the example programs.
# $Id: do_comps.sh,v 1.1 2006/06/27 17:44:54 ed Exp $

set -e
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
