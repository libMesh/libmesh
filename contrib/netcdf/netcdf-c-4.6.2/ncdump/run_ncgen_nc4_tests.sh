#!/bin/sh
# This shell script runs the ncdump tests.
# Dennis Heimbigner

if test "x$srcdir" = x ; then srcdir="."; fi
. ../test_common.sh

##
# Function to test a netCDF CDL file.
# First generate binary nc. Then use ncdump to compare against
# original CDL file.
# Input: CDL file name, minus the suffix, output filename
# Other input: arguments.
#
# Example:
#     $ validateNC compound_datasize_test -k nc4
##
validateNC() {
    ORIGNAME=$1
    BASENAME=tst_$1_run_ncgen_nc4_tests
    INFILE=$top_srcdir/ncgen/$1.cdl
    TMPFILE=tst_$2.cdl
    shift
    shift
    ARGS=$@

    echo "*** generating $BASENAME.nc ***"
    ${NCGEN} $ARGS -o $BASENAME.nc $INFILE
    ${NCDUMP} -n $ORIGNAME $BASENAME.nc | sed 's/e+0/e+/g' > $TMPFILE
    echo "*** comparing binary against source CDL file *** "
    diff -b -w $INFILE $TMPFILE
}

echo "*** Testing ncgen for netCDF-4."
set -e

echo "*** creating netCDF-4 file c0_4.nc from c0_4.cdl..."
validateNC "c0_4" "c0_4" -k nc4 -b -o c0_4.nc

echo "*** creating netCDF-4 classic model file c0_4c.nc from c0.cdl..."
validateNC "c0" "c0_4c" -k nc7 -b

echo "*** creating C code for CAM file ref_camrun.cdl..."
${NCGEN} -lc $top_srcdir/ncgen/ref_camrun.cdl > camrun.c

echo "*** test for jira NCF-199 bug"
validateNC "ncf199" "ncf199" -k nc4

if test "x$NC_VLEN_NOTEST" = x ; then
echo "*** creating binary files for github issue 323..."
echo "*** github issue 323 test 1"
validateNC "compound_datasize_test" "compound_datasize_test" -k nc4
fi

echo "*** github issue 323 test 2"
validateNC "compound_datasize_test2" "compound_datasize_test2"  -k nc4

echo "*** Global Attribute with Enum type"
validateNC "tst_gattenum" "tst_gattenum"  -k nc4

echo "*** Integer constant with just 'u' suffix"
validateNC "tst_usuffix" "tst_usuffix"  -k nc4

echo "*** Test successful!"
exit 0
