#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# For a netCDF-4 build, test nccopy on netCDF files in this directory

if test -f tst_comp2${ext} ; then ${execdir}/tst_comp2 ; fi


set -e
echo ""

# These files are actually in $srcdir in distcheck builds, so they
# need to be handled differently.
# ref_tst_compounds2 ref_tst_compounds3 ref_tst_compounds4
TESTFILES='tst_comp tst_comp2 tst_enum_data tst_fillbug
 tst_group_data tst_nans tst_opaque_data tst_solar_1 tst_solar_2
 tst_solar_cmp tst_special_atts tst_string_data tst_unicode'

if test "x$NC_VLEN_NOTEST" = x ; then
TESFILES="$TESTFILES tst_vlen_data"
fi

echo "*** Testing netCDF-4 features of nccopy on ncdump/*.nc files"
for i in $TESTFILES ; do
    echo "*** Test nccopy $i.nc copy_of_$i.nc ..."
    ${NCCOPY} $i.nc copy_of_$i.nc
${NCDUMP} -n copy_of_$i $i.nc > tmp.cdl
${NCDUMP} copy_of_$i.nc > copy_of_$i.cdl
#    echo "*** compare " with copy_of_$i.cdl
    diff copy_of_$i.cdl tmp.cdl
    rm copy_of_$i.nc copy_of_$i.cdl tmp.cdl
done
# echo "*** Testing compression of deflatable files ..."
./tst_compress
echo "*** Test nccopy -d1 can compress a classic format file ..."
${NCCOPY} -d1 tst_inflated.nc tst_deflated.nc
if test `wc -c < tst_deflated.nc` -ge  `wc -c < tst_inflated.nc`; then
    exit 1
fi
echo "*** Test nccopy -d1 can compress a netCDF-4 format file ..."
${NCCOPY} -d1 tst_inflated4.nc tst_deflated.nc
if test `wc -c < tst_deflated.nc` -ge  `wc -c < tst_inflated4.nc`; then
    exit 1
fi
echo "*** Test nccopy -d1 -s can compress a classic model netCDF-4 file even more ..."
${NCCOPY} -d1 -s tst_inflated.nc tmp.nc
if test `wc -c < tmp.nc` -ge  `wc -c < tst_inflated.nc`; then
    exit 1
fi
echo "*** Test nccopy -d1 -s can compress a netCDF-4 file even more ..."
${NCCOPY} -d1 -s tst_inflated4.nc tmp.nc
if test `wc -c < tmp.nc` -ge  `wc -c < tst_inflated4.nc`; then
    exit 1
fi
echo "*** Test nccopy -d0 turns off compression, shuffling of compressed, shuffled file ..."
${NCCOPY} -d0 tst_inflated4.nc tmp.nc
${NCDUMP} -sh tmp.nc > tmp.cdl
if fgrep '_DeflateLevel' < tmp.cdl ; then
    exit 1
fi
if fgrep '_Shuffle' < tmp.cdl ; then
    exit 1
fi
rm tst_deflated.nc tst_inflated.nc tst_inflated4.nc tmp.nc tmp.cdl

echo "*** Testing nccopy -d1 -s on ncdump/*.nc files"
for i in $TESTFILES ; do
    echo "*** Test nccopy -d1 -s $i.nc copy_of_$i.nc ..."
    ${NCCOPY} -d1 -s $i.nc copy_of_$i.nc
${NCDUMP} -n copy_of_$i $i.nc > tmp.cdl
${NCDUMP} copy_of_$i.nc > copy_of_$i.cdl
#    echo "*** compare " with copy_of_$i.cdl
    diff copy_of_$i.cdl tmp.cdl
    rm copy_of_$i.nc copy_of_$i.cdl tmp.cdl
done
./tst_chunking
echo "*** Test that nccopy -c can chunk and unchunk files"
${NCCOPY} -M0 tst_chunking.nc tmp.nc
${NCDUMP} tmp.nc > tmp.cdl
${NCCOPY} -c dim0/,dim1/1,dim2/,dim3/1,dim4/,dim5/1,dim6/ tst_chunking.nc tmp-chunked.nc
${NCDUMP} -n tmp tmp-chunked.nc > tmp-chunked.cdl
diff tmp.cdl tmp-chunked.cdl
${NCCOPY} -c dim0/,dim1/,dim2/,dim3/,dim4/,dim5/,dim6/ tmp-chunked.nc tmp-unchunked.nc
${NCDUMP} -n tmp tmp-unchunked.nc > tmp-unchunked.cdl
diff tmp.cdl tmp-unchunked.cdl
${NCCOPY} -c // tmp-chunked.nc tmp-unchunked.nc
${NCDUMP} -n tmp tmp-unchunked.nc > tmp-unchunked.cdl
diff tmp.cdl tmp-unchunked.cdl
echo "*** Test that nccopy -c works as intended for record dimension default (1)"
${NCGEN} -b -o tst_bug321.nc $srcdir/tst_bug321.cdl
${NCCOPY} -k nc7 -c"lat/2,lon/2" tst_bug321.nc tmp.nc
${NCDUMP} -n tst_bug321 tmp.nc > tmp.cdl
diff -b $srcdir/tst_bug321.cdl tmp.cdl
# echo "*** Test that nccopy compression with chunking can improve compression"
rm tst_chunking.nc tmp.nc tmp.cdl tmp-chunked.nc tmp-chunked.cdl tmp-unchunked.nc tmp-unchunked.cdl

echo "*** All nccopy tests passed!"
exit 0
