#!/bin/sh
# This shell script runs extra tests ncdump for netcdf-4
# Dennis Heimbigner, Ward Fisher

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

echo ""
echo "*** Running extra netcdf-4 tests."

#
# In windows, these tests fail because srcdir is prepended.
# e.g., Instead of 'ncdump ref_tst_compounds2' the file would
# contain:
#        'ncdump ./ref_tst_compounds2'. This causes the test to fail.
# But, 'srcdir' is necessary for make distcheck.
#
# Short term solution, use sed when on windows/MSYS to
# remove the './','../../ncdump'.
#
# I am undoing this because libdispatch/dwinpath.c
# should be taking care of this. If not, then that is
# what we need to fix. Alternatively, we can use top_srcdir,
# which is an absolute path

echo "*** running tst_string_data to create test files..."
${execdir}/tst_string_data

echo "*** dumping tst_string_data.nc to tst_string_data.cdl..."
${NCDUMP} tst_string_data.nc > tst_string_data.cdl
cleanncprops tst_string_data.cdl tst_string_data.tmp
cleanncprops ${srcdir}/ref_tst_string_data.cdl ref_tst_string_data.tmp
echo "*** comparing tst_string_data.cdl with ref_tst_string_data.cdl..."
diff -b tst_string_data.tmp ref_tst_string_data.tmp

#echo '*** testing non-coordinate variable of same name as dimension...'
#${NCGEN} -v4 -b -o tst_noncoord.nc ${top_srcdir}/ncdump/ref_tst_noncoord.cdl

echo '*** testing reference file ref_tst_compounds2.nc...'
${NCDUMP} ${top_srcdir}/ncdump/ref_tst_compounds2.nc > tst_compounds2.cdl
diff -b tst_compounds2.cdl ${top_srcdir}/ncdump/ref_tst_compounds2.cdl

echo '*** testing reference file ref_tst_compounds3.nc...'
${NCDUMP} ${top_srcdir}/ncdump/ref_tst_compounds3.nc > tst_compounds3.cdl
diff -b tst_compounds3.cdl ${top_srcdir}/ncdump/ref_tst_compounds3.cdl

echo '*** testing reference file ref_tst_compounds4.nc...'
${NCDUMP} ${top_srcdir}/ncdump/ref_tst_compounds4.nc > tst_compounds4.cdl
diff -b tst_compounds4.cdl ${top_srcdir}/ncdump/ref_tst_compounds4.cdl

# Exercise Jira NCF-213 bug fix
#    rm -f tst_ncf213.cdl tst_ncf213.nc
# Remove specific _NCProperties values
${NCGEN} -b -o tst_ncf213.nc $srcdir/ref_tst_ncf213.cdl
${NCDUMP} -s -h tst_ncf213.nc > tst_ncf213.cdl
cleanncprops tst_ncf213.cdl tst_ncf213.tmp
cleanncprops ${srcdir}/ref_tst_ncf213.cdl ref_tst_ncf213.tmp
# Now compare
ok=1;
if diff -b tst_ncf213.tmp ref_tst_ncf213.tmp ; then ok=1; else ok=0; fi

# cleanup
#rm -f tst_ncf213.cdl tst_ncf213.nc
rm -f *.tmp

if test $ok = 0 ; then
  echo "*** FAIL: NCF-213 Bug Fix test"
  exit 1
fi

echo "*** All ncgen and ncdump extra test output for netCDF-4 format passed!"
exit 0
