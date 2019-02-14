#!/bin/sh

export SETX=1
if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# Which test cases to exercise
API=1
NG=1
NCP=1
UNK=1
NGC=1
MISC=1

# Load the findplugins function
. ${builddir}/findplugin.sh
echo "findplugin.sh loaded"

# Function to remove selected -s attributes from file;
# These attributes might be platform dependent
sclean() {
cat $1 \
  | sed -e '/var:_Endianness/d' \
  | sed -e '/_NCProperties/d' \
  | sed -e '/_SuperblockVersion/d' \
  | sed -e '/_IsNetcdf4/d' \
  | cat > $2
}

# Function to extract _Filter attribute from a file
# These attributes might be platform dependent
getfilterattr() {
sed -e '/var:_Filter/p' -ed <$1 >$2
}

trimleft() {
sed -e 's/[ 	]*\([^ 	].*\)/\1/' <$1 >$2
}

# Locate the plugin path and the library names; argument order is critical
# Find bzip2 and capture
findplugin bzip2
BZIP2PATH="${HDF5_PLUGIN_PATH}/${HDF5_PLUGIN_LIB}"
# Find misc and capture
findplugin misc
MISCPATH="${HDF5_PLUGIN_PATH}/${HDF5_PLUGIN_LIB}"

echo "final HDF5_PLUGIN_PATH=${HDF5_PLUGIN_PATH}"
export HDF5_PLUGIN_PATH

# Verify
if ! test -f ${BZIP2PATH} ; then echo "Unable to locate ${BZIP2PATH}"; exit 1; fi
if ! test -f ${MISCPATH} ; then echo "Unable to locate ${MISCPATH}"; exit 1; fi

# Execute the specified tests

if test "x$API" = x1 ; then
echo "*** Testing dynamic filters using API"
rm -f ./bzip2.nc ./bzip2.dump ./tst_filter.txt
${execdir}/test_filter
${NCDUMP} -s bzip2.nc > ./tst_filter.txt
# Remove irrelevant -s output
sclean ./tst_filter.txt ./bzip2.dump
diff -b -w ${srcdir}/bzip2.cdl ./bzip2.dump
echo "*** Pass: API dynamic filter"
fi

if test "x$MISC" = x1 ; then
echo
echo "*** Testing dynamic filters parameter passing"
rm -f ./testmisc.nc tst_filter.txt tst_filter2.txt
${execdir}/test_filter_misc
# Verify the parameters via ncdump
${NCDUMP} -s testmisc.nc > ./tst_filter.txt
# Extract the parameters
getfilterattr ./tst_filter.txt ./tst_filter2.txt
rm -f ./tst_filter.txt
trimleft ./tst_filter2.txt ./tst_filter.txt
rm -f ./tst_filter2.txt
cat >./tst_filter2.txt <<EOF
var:_Filter = "32768,1,239,23,65511,27,77,93,1145389056,3287505826,1097305129,1,2147483648,4294967295,4294967295" ;
EOF
diff -b -w ./tst_filter.txt ./tst_filter2.txt
echo "*** Pass: parameter passing"
fi

if test "x$NG" = x1 ; then
echo "*** Testing dynamic filters using ncgen"
rm -f ./bzip2.nc ./bzip2.dump ./tst_filter.txt
${NCGEN} -lb -4 -o bzip2.nc ${srcdir}/bzip2.cdl
${NCDUMP} -s bzip2.nc > ./tst_filter.txt
# Remove irrelevant -s output
sclean ./tst_filter.txt ./bzip2.dump
diff -b -w ${srcdir}/bzip2.cdl ./bzip2.dump
echo "*** Pass: ncgen dynamic filter"
fi

if test "x$NCP" = x1 ; then
echo "*** Testing dynamic filters using nccopy"
rm -f ./unfiltered.nc ./filtered.nc ./tmp.nc ./filtered.dump ./tst_filter.txt
${NCGEN} -4 -lb -o unfiltered.nc ${srcdir}/unfiltered.cdl
echo "	*** Testing simple filter application"
${NCCOPY} -M0 -F "/g/var,307,9,4" unfiltered.nc filtered.nc
${NCDUMP} -s filtered.nc > ./tst_filter.txt
# Remove irrelevant -s output
sclean ./tst_filter.txt ./filtered.dump
diff -b -w ${srcdir}/filtered.cdl ./filtered.dump
echo "	*** Pass: nccopy simple filter"

echo "	*** Testing pass-thru of filters"
rm -f ./tst_filter.txt tst_filter2.txt ./tst_filter2.nc
# Prevent failure by allowing any chunk size
${NCCOPY} -M0 ./filtered.nc ./tst_filter2.nc
${NCDUMP} -s tst_filter2.nc > ./tst_filter.txt
sed -e '/_Filter/p' -e d < ./tst_filter.txt >tst_filter2.txt
test -s tst_filter2.txt
echo "	*** Pass: pass-thru of filters"

echo "	*** Testing -F none"
rm -f ./tst_filter.txt ./tst_filter2.txt ./tst_filter.nc
${NCCOPY} -M0 -F none ./filtered.nc ./tst_filter.nc
${NCDUMP} -s tst_filter.nc > ./tst_filter.txt
sed -e '/_Filter/p' -e d < ./tst_filter.txt >./tst_filter2.txt
test ! -s tst_filter2.txt
echo "	*** Pass: -F none"

echo "	*** Testing -F var,none "
rm -f ./tst_filter.txt ./tst_filter.nc
${NCCOPY} -M0 -F "/g/var,none" ./filtered.nc ./tst_filter.nc
${NCDUMP} -s tst_filter.nc > ./tst_filter.txt
sed -e '/_Filter/p' -e d < ./tst_filter.txt >tst_filter2.txt
test ! -s tst_filter2.txt
echo "	*** Pass: -F var,none"

echo "*** Pass: all nccopy filter tests"
fi

if test "x$UNK" = x1 ; then
echo "*** Testing access to filter info when filter dll is not available"
rm -f bzip2.nc ./tst_filter.txt
# build bzip2.nc
${NCGEN} -lb -4 -o bzip2.nc ${srcdir}/bzip2.cdl
# dump and clean bzip2.nc header only when filter is avail
${NCDUMP} -hs bzip2.nc > ./tst_filter.txt
# Remove irrelevant -s output
sclean ./tst_filter.txt bzip2.dump
# Now hide the filter code
mv ${BZIP2PATH} ${BZIP2PATH}.save
# dump and clean bzip2.nc header only when filter is not avail
rm -f ./tst_filter.txt
${NCDUMP} -hs bzip2.nc > ./tst_filter.txt
# Remove irrelevant -s output
sclean ./tst_filter.txt bzip2x.dump
# Restore the filter code
mv ${BZIP2PATH}.save ${BZIP2PATH}
diff -b -w ./bzip2.dump ./bzip2x.dump
echo "*** Pass: ncgen dynamic filter"
fi

if test "x$NGC" = x1 ; then
rm -f ./test_bzip2.c
echo "*** Testing dynamic filters using ncgen with -lc"
${NCGEN} -lc -4 ${srcdir}/bzip2.cdl > test_bzip2.c
diff -b -w ${srcdir}/ref_bzip2.c ./test_bzip2.c
echo "*** Pass: ncgen dynamic filter"
fi

#cleanup
rm -f ./bzip*.nc ./unfiltered.nc ./filtered.nc ./tst_filter.txt ./tst_filter2.txt *.dump bzip*hdr.*
rm -f ./test_bzip2.c
rm -f ./testmisc.nc
rm -f ./tst_filter2.nc
echo "*** Pass: all selected tests passed"

exit 0
