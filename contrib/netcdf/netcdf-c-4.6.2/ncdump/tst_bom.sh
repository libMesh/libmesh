#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh

# This shell script tests BOM support in ncgen

set -e

# add hack for sunos
export srcdir;

echo ""

rm -f tst_bom.cdl tmp.cdl tst_bom8.* tst_bom16.*

cat <<EOF >>tst_bom.cdl
netcdf tst_bom {
variables:
  float f;
data:

  f = 1;
}
EOF

echo "*** Generate a cdl file with leading UTF-8 BOM."
${execdir}/bom 8 >tst_bom8.cdl
cat tst_bom.cdl >> tst_bom8.cdl

echo "*** Verify .nc file"
${NCGEN} -k nc3 -o tst_bom8.nc tst_bom8.cdl
${NCDUMP} -n tst_bom tst_bom8.nc > tmp.cdl
diff -w tst_bom.cdl tmp.cdl

# Do it again but with Big-Endian 16; should fail

rm -f tmp.cdl tst_bom8.* tst_bom16.*

echo "*** Generate a cdl file with leading UTF-16 BOM."
${execdir}/bom 16 >tst_bom16.cdl
cat tst_bom.cdl >> tst_bom16.cdl

echo "*** Verify UTF-16 file fails"
if ${NCGEN} -k nc3 -o tst_bom16.nc tst_bom16.cdl ; then
echo 'BOM Big Endian 16 succeeded, but should not'
exit 1
else
echo '***XFAIL : BOM Big Endian 16'
fi

# Cleanup
rm -f tst_bom.cdl tmp.cdl tst_bom8.* tst_bom16.*

exit 0
