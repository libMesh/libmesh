# This shell script runs ncdump tests relating to the new UTF8 name stuff.
# Russ Rew, Dennis Heimbigner, Ward Fisher, Ed Hartnett
#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh
set -e

echo ""
echo "*** Testing ncgen and ncdump for UTF8 support..."

# Run tst_utf8.c to produce test file tst_utf8.nc.
${execdir}/tst_utf8

# if test "x$builddir" = x ; then
# builddir="."
# fi
# if test "x$srcdir" = x ; then
# srcdir="."
# fi

rm -f tst_utf8_nc3.nc tst_utf8_nc3.cdl
echo "*** creating classic offset file with utf8 characters..."
${NCGEN} -b -o tst_utf8_nc3.nc ${srcdir}/ref_tst_utf8.cdl
echo "*** dump and compare utf8 output..."
${NCDUMP} -n utf8 tst_utf8_nc3.nc > tst_utf8_nc3.cdl
diff -b -w tst_utf8_nc3.cdl ${srcdir}/ref_tst_utf8.cdl

rm -f tst_utf8_64.nc tst_utf8_64.cdl
echo "*** creating 64-bit offset file with utf8 characters..."
${NCGEN} -k 64-bit-offset -b -o tst_utf8_64.nc ${srcdir}/ref_tst_utf8.cdl
echo "*** (64 bit) dump and compare utf8 output..."
${NCDUMP} -n utf8 tst_utf8_64.nc > tst_utf8_64.cdl
diff -b -w tst_utf8_64.cdl ${srcdir}/ref_tst_utf8.cdl

echo "*** dumping tst_utf8.nc to tst_utf8.cdl..."
rm -f tst_tst8.cdl tst_utf8_tst8.cdl
sed -e 's/^netcdf tst_unicode/netcdf tst_utf8/' <${srcdir}/ref_tst_unicode.cdl > tst_tst8.cdl
${NCDUMP} tst_utf8.nc > tst_utf8_tst8.cdl
echo "*** comparing tst_utf8_tst8.cdl with tst_tst8.cdl..."
diff -b -w tst_utf8_tst8.cdl tst_tst8.cdl
echo "*** All utf8 tests of ncgen and ncdump passed!"
exit 0
