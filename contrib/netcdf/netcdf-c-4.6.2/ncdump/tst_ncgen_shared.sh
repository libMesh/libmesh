#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

# To add a new test,
# 1. put the .cdl file in the 'cdl' directory
# 2. put the result of running ncgen then ncdump
#    into the directory 'expected' as .dmp
# 3. Modify the file tst_ncgen_shared.sh to add
#    the test to the end of the TESTS4 variable
#    or CLASSIC variable.
# 4. Add the new files into cdl/Makfile.am
#    and expected/Makefile.am 

set -e
RESULTSDIR="./results_$$"
#SHOWXFAILS=1

# Locate the cdl and expected directory
cdl="${srcdir}/cdl"
expected="${srcdir}/expected"

case "x${KFLAG}" in
x1) CLASSIC=1; MODE=3;;
x2) CLASSIC=1; MODE=3;;
x3) CLASSIC=0; MODE=4;;
x4) CLASSIC=1; MODE=3;;
*) echo "illegal KFLAG" ; exit 1;;
esac

# Define the set of tests that can be
# processed with either the -k nc3 or -k nc4 or -k nc7 flag

# The netcdf-3 tests are divided into two parts
# These test can be run when --enable-netcdf-4 is false
CLASSIC3="\
nc_enddef \
ref_tst_unicode \
ref_tst_utf8 \
simple_xy \
small \
nc_sync \
ref_tst_small \
small2 \
tst_ncml
n3time \
ref_tst_chardata \
ref_tst_nul3 \
ref_tst_long_charconst \
tst_chararray \
unlimtest1 \
ref_keyword"

NONCLASSIC3="\
test0 \
sfc_pres_temp \
fills \
c0 \
example_good \
pres_temp_4D \
ref_nctst \
ref_nctst_64bit_offset \
ref_ctest1_nc4 \
ref_ctest1_nc4c \
ref_nctst_netcdf4 \
ref_nctst_netcdf4_classic \
ref_tst_unlim2 \
ref_tst_names \
"

if test "${CLASSIC}" = "1" ; then
TESTS3="${CLASSIC3}"
else
TESTS3="${CLASSIC3} ${NONCLASSIC3}"
fi

# Define the set of tests that must be
# processed with the -k nc4 flag

TESTS4="\
ref_dimscope \
ref_typescope \
ref_tst_string_data \
ref_tst_comp \
ref_tst_comp2 \
ref_tst_comp3 \
ref_tst_group_data \
ref_tst_opaque_data \
ref_tst_solar_1 \
ref_tst_solar_2 \
ref_tst_enum_data \
ref_tst_special_atts \
ref_tst_nans \
ref_solar \
unlimtest2 \
ref_niltest \
ref_tst_h_scalar \
ref_tst_nul4 \
"

if test "x$NC_VLEN_NOTEST" = x ; then
TESTS4="$TESTS4 ref_tst_vlen_data ref_tst_vlen_data2"
fi

SPECIALTESTS3="ref_tst_special_atts3"

SPECIALTESTS="${SPECIALTESTS3} ref_tst_special_atts"

XFAILTESTS=""
# Fails because ncdump does not output multiple unlim char types correctly
XFAILTESTS="ref_tst_unlim2 $XFAILTESTS"
# Fails because ?
XFAILTESTS="ref_const_test $XFAILTESTS"
# Fails because ?
XFAILTESTS="ref_tst_chardata $XFAILTESTS"
# Fails because ncdump is crashing
#XFAILTESTS="ref_tst_econst $XFAILTESTS"

# Following are generally not run
# Because of the size of their output
BIGTESTS3="\
bigf1 \
bigf2 \
bigr1 \
bigr2"

# Following deliberately produces a very large
# file: too large for netcdf to handle
# Currently not used because of space and time
# constraints
XFAILBIG="bigf3"

BIGTESTS4="ref_tst_solar_1"

# This test is both big and slow
# File was too large to reasonably include
# so I removed it
#BIGBIG3="gfs1"

BIGTESTS="${BIGTESTS3} ${BIGTESTS4} ${BIGBIG3}"

failcount=0
passcount=0
xfailcount=0

rm -fr $RESULTSDIR
