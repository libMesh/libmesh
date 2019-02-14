#!/bin/sh
# This shell script tests ncdump for netcdf-4
# Ed Hartnett, Dennis Heimbigner, Ward Fisher

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e

# Remove the version information from _NCProperties
cleanncprops() {
  src="$1"
  dst="$2"
  rm -f $dst
  cat $src \
  | sed -e 's/_SuperblockVersion = 1/_SuperblockVersion = 0/' \
  | sed -e 's/\(netcdflibversion\|netcdf\)=.*|/\1=NNNN|/' \
  | sed -e 's/\(hdf5libversion\|hdf5\)=.*"/\1=HHHH"/' \
  | grep -v '_NCProperties' \
  | cat >$dst
}

ERR() {
    RES=$?
    if [ $RES -ne 0 ]; then
        echo "Error found: $RES"
        exit $RES
    fi
}

echo ""
echo "*** Testing ncgen and ncdump for netCDF-4 format."
${NCGEN} -k nc4 -b -o tst_netcdf4_c0_4.nc ${ncgenc04} ;ERR
${NCDUMP} -n c1 tst_netcdf4_c0_4.nc | sed 's/e+0/e+/g' > tst_netcdf4_c1_4.cdl ; ERR
diff -b tst_netcdf4_c1_4.cdl $srcdir/ref_ctest1_nc4.cdl ; ERR

echo "*** Creating test output tst_netcdf4_c0.nc."
${NCGEN} -k nc7 -b -o tst_netcdf4_c0.nc ${ncgenc0} ; ERR

echo "*** Testing that program tst_h_rdc0 can read tst_netcdf4_c0.nc."
${execdir}/tst_h_rdc0 ; ERR

echo "*** Running tst_create_files.c to create test files."
${execdir}/tst_create_files ; ERR
echo "*** Testing tst_create_files output for netCDF-4 features."
${NCDUMP} tst_solar_1.nc | sed 's/e+0/e+/g' > tst_solar_1.cdl ; ERR
diff -b tst_solar_1.cdl $srcdir/ref_tst_solar_1.cdl ; ERR
${NCDUMP} tst_solar_2.nc | sed 's/e+0/e+/g' > tst_solar_2.cdl ; ERR
diff -b tst_solar_2.cdl $srcdir/ref_tst_solar_2.cdl ; ERR

echo "*** Running tst_group_data.c to create test files."
${execdir}/tst_group_data ; ERR
${NCDUMP} tst_group_data.nc | sed 's/e+0/e+/g' > tst_group_data.cdl ; ERR
diff -b tst_group_data.cdl $srcdir/ref_tst_group_data.cdl ; ERR

echo "*** Testing -v option with absolute name and groups..."
${NCDUMP} -v g2/g3/var tst_group_data.nc | sed 's/e+0/e+/g' > tst_group_data.cdl ; ERR
diff -b tst_group_data.cdl $srcdir/ref_tst_group_data_v23.cdl ; ERR

echo "*** Testing -v option with relative name and groups..."
${NCDUMP} -v var,var2 tst_group_data.nc | sed 's/e+0/e+/g' > tst_group_data.cdl ; ERR
diff -b tst_group_data.cdl $srcdir/ref_tst_group_data.cdl ; ERR

echo "*** Running tst_enum_data.c to create test files."
${execdir}/tst_enum_data ; ERR
${NCDUMP} tst_enum_data.nc | sed 's/e+0/e+/g' > tst_enum_data.cdl ; ERR
diff -b tst_enum_data.cdl $srcdir/ref_tst_enum_data.cdl ; ERR

echo "*** Running tst_opaque_data.c to create test files."
${execdir}/tst_opaque_data ; ERR
${NCDUMP} tst_opaque_data.nc | sed 's/e+0/e+/g' > tst_opaque_data.cdl ; ERR
diff -b tst_opaque_data.cdl $srcdir/ref_tst_opaque_data.cdl ; ERR

if test "x$NC_VLEN_NOTEST" = x ; then
echo "*** Running tst_vlen_data.c to create test files."
${execdir}/tst_vlen_data ; ERR
${NCDUMP} tst_vlen_data.nc | sed 's/e+0/e+/g' > tst_vlen_data.cdl ; ERR
diff -b tst_vlen_data.cdl $srcdir/ref_tst_vlen_data.cdl ; ERR
fi

echo "*** Running tst_comp.c to create test files."
${execdir}/tst_comp ; ERR
${NCDUMP} tst_comp.nc | sed 's/e+0/e+/g' > tst_comp.cdl ; ERR
diff -b tst_comp.cdl $srcdir/ref_tst_comp.cdl ; ERR

echo "*** Running tst_nans.c to create test files."
${execdir}/tst_nans ; ERR
${NCDUMP} tst_nans.nc | sed 's/e+0/e+/g' > tst_nans.cdl ; ERR
diff -b tst_nans.cdl $srcdir/ref_tst_nans.cdl ; ERR

# Do unicode test only if it exists => BUILD_UTF8 is true
if test -f ./tst_unicode -o -f ./tst_unicode.exe ; then
  echo "*** dumping tst_unicode.nc to tst_unicode.cdl..."
  ${execdir}/tst_unicode ; ERR
${NCDUMP} tst_unicode.nc | sed 's/e+0/e+/g' > tst_unicode.cdl ; ERR
  #echo "*** comparing tst_unicode.cdl with ref_tst_unicode.cdl..."
  #diff -b tst_unicode.cdl $srcdir/ref_tst_unicode.cdl
fi

echo "*** Running tst_special_atts.c to create test files."
${execdir}/tst_special_atts ; ERR
${NCDUMP} -c -s tst_special_atts.nc  > tst_special_atts.cdl ; ERR
cleanncprops tst_special_atts.cdl tst_special_atts.tmp
cleanncprops $srcdir/ref_tst_special_atts.cdl ref_tst_special_atts.tmp
echo "*** comparing tst_special_atts.cdl with ref_tst_special_atts.cdl..."
diff -b tst_special_atts.tmp ref_tst_special_atts.tmp ; ERR

#echo ""
#echo "*** Testing ncdump on file with corrupted header "
#rm -f ./ignore_tst_netcdf4
#if ${NCDUMP} ${srcdir}/ref_test_corrupt_magic.nc > ./ignore_tst_netcdf4 2>&1 ; then
#echo "***Fail: ncdump should have failed on ref_test_corrupt_magic.nc"
#else
#echo "***XFail: ncdump properly failed on ref_test_corrupt_magic.nc"
#fi
#rm -fr ./ignore_tst_netcdf4

# This should work, but does not. See github issue 982.
#echo "*** creating tst_output_irish_rover.cdl from ref_tst_irish_rover.nc..."
#${NCDUMP} ref_tst_irish_rover.nc > tst_output_irish_rover.cdl

rm -f *.tmp
echo "*** All ncgen and ncdump test output for netCDF-4 format passed!"
exit 0
