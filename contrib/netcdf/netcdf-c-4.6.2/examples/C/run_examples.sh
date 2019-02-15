#!/bin/sh
# This shell script runs the examples.
# Ed Hartnett

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../../test_common.sh

echo "*** Running examples."
set -e

echo "*** running simple_xy examples..."
${execdir}/simple_xy_wr
${execdir}/simple_xy_rd

echo "*** running sfc_pres_temp examples..."
${execdir}/sfc_pres_temp_wr
${execdir}/sfc_pres_temp_rd

echo "*** running pres_temp_4D examples..."
${execdir}/pres_temp_4D_wr
${execdir}/pres_temp_4D_rd

echo "*** Examples successful!"
exit 0
