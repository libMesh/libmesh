#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# Get the target OS and CPU
CPU=`uname -p`
OS=`uname`

# Test diskless on a reasonably large file size

# Try a large in-memory file; two instances at once may be needed when doing realloc
#SIZE=1000000000
SIZE=500000000

FILE4=tst_diskless4.nc

# Uncomment to get timing
#TIME=time

# Create the reference ncdump output for tst_diskless4
rm -fr ref_tst_diskless4.cdl
cat >ref_tst_diskless4.cdl <<EOF
netcdf tst_diskless4 {
dimensions:
	dim = 500000000 ;
variables:
	byte var0(dim) ;
}
EOF

echo ""
rm -f $FILE4
$TIME ./tst_diskless4 $SIZE create
# Validate it
${NCDUMP} -h $FILE4 |diff -w - ref_tst_diskless4.cdl

echo ""
rm -f $FILE4
$TIME ./tst_diskless4 $SIZE creatediskless
# Validate it
${NCDUMP} -h $FILE4 |diff -w - ref_tst_diskless4.cdl

echo ""
$TIME ./tst_diskless4 $SIZE open

echo ""
$TIME ./tst_diskless4 $SIZE opendiskless

# cleanup
rm -f $FILE4 tst_diskless4.cdl ref_tst_diskless4.cdl

exit
