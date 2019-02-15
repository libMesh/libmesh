#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# If we want to run valgrind
#NCCOPY="valgrind --leak-check=full ${NCCOPY}"

# Choose tests to run
T1=1
T2=1
T3=1
T4=1
T5=1

# For a netCDF-4 build, test nccopy chunking rules

set -e
echo ""

# Trim off leading and trailing whitespace
# Also remove any <cr>
# usage: trim <line>
# Leaves result in variable TRIMMED
trim() {
    # trim leading whitespace and remove <cr>
    TMP=`echo "$1" |tr -d '\r' | sed -e 's/^[ 	]*//'`
    # trim trailing whitespace
    TRIMMED=`echo "$TMP" | sed -e 's/[ 	]*$//'`
}

# usage: checkfvar <file>
checkfvar() {
  # Make sure that fvar was not chunked
  C5FVAR=`sed -e '/fvar:_ChunkSizes/p' -e d <$1`
  if test "x$C5FVAR" != x ; then
      echo "***Fail: fvar was chunked"
      exit 1
  fi
}

# usage: checkivar <file>
checkivar() {
  # Make sure that ivar was not chunked
  C5IVAR=`sed -e '/ivar:_ChunkSizes/p' -e d <$1`
  if test "x$C5IVAR" != x ; then
      echo "***Fail: ivar was chunked"
      exit 1
  fi
}

# usage: verifychunkline line1 line2
verifychunkline() {
    # trim leading whitespace
    trim "$1"; L1="$TRIMMED"
    trim "$2"; L2="$TRIMMED"
    if test "x$L1" != "x$L2" ; then
	echo "chunk line mismatch |$L1| |$L2|"
	exit 1;
    fi
}

# Remove any temporary files
cleanup() {
    rm -f tmp_nc5.nc tmp_nc5a.nc
    rm -f tmp_nc5.cdl tmp_nc5a.cdl tmp_nc5b.cdl
    rm -f tmp_nc5_omit.nc tmp_nc5_omit.cdl
}

# remove all created files
reset() {
    cleanup
    rm -fr tst_nc5.nc tst_nc5.cdl
    rm -f tst_nc5_omit.nc tst_nc5_omit.cdl
}

reset

if test "x$T1" = x1 ; then

# Create a simple classic input file 
./tst_chunking tst_nc5.nc

# Save a .cdl version
${NCDUMP} tst_nc5.nc > tst_nc5.cdl

echo "*** Test nccopy -c with per-variable chunking; classic->enhanced"
# This should produce same as -c dim0/,dim1/1,dim2/,dim3/1,dim4/,dim5/1,dim6/
${NCCOPY} -c ivar:7,1,2,1,5,1,9 tst_nc5.nc tmp_nc5.nc
${NCDUMP} -n tst_nc5 tmp_nc5.nc > tmp_nc5.cdl
# Verify that the core cdl is the same
diff tst_nc5.cdl tmp_nc5.cdl

# Look at the output chunking of ivar
rm -f tmp_nc5a.cdl # reuse
${NCDUMP} -hs -n tst_nc5 tmp_nc5.nc > tmp_nc5.cdl
# extract the chunking line
TESTLINE=`sed -e '/ivar:_ChunkSizes/p' -e d <tmp_nc5.cdl`
# track line to match
BASELINE='ivar:_ChunkSizes = 7, 1, 2, 1, 5, 1, 9 ;'
verifychunkline "$TESTLINE" "$BASELINE"

# Make sure that fvar was not chunked
checkfvar tmp_nc5.cdl

fi # T1

if test "x$T2" = x1 ; then

echo "*** Test nccopy -c with per-variable chunking; enhanced->enhanced"
reset
./tst_chunking tst_nc5.nc deflate
${NCDUMP} -n tst_nc5 tst_nc5.nc > tst_nc5.cdl
${NCCOPY} -c ivar:4,1,2,1,5,2,3 tst_nc5.nc tmp_nc5.nc
${NCDUMP} -n tst_nc5 tmp_nc5.nc > tmp_nc5.cdl
diff tst_nc5.cdl tmp_nc5.cdl

# Look at the output chunking
rm -f tmp_nc5.cdl # reuse
${NCDUMP} -hs -n tst_nc5 tmp_nc5.nc > tmp_nc5.cdl
# extract the chunking line
TESTLINE=`sed -e '/ivar:_ChunkSizes/p' -e d <tmp_nc5.cdl`
# track line to match
BASELINE='ivar:_ChunkSizes = 4, 1, 2, 1, 5, 2, 3 ;'
verifychunkline "$TESTLINE" "$BASELINE"

# Make sure that fvar was not chunked
checkfvar tmp_nc5.cdl

fi # T2

if test "x$T3" = x1 ; then

echo "*** Test nccopy -c with FQN var name; enhanced ->enhanced"
reset
./tst_chunking tst_nc5.nc group
${NCDUMP} -n tst_nc5 tst_nc5.nc > tst_nc5.cdl
${NCCOPY} -c /g/ivar:4,1,2,1,5,2,3 tst_nc5.nc tmp_nc5.nc
${NCDUMP} -n tst_nc5 tmp_nc5.nc > tmp_nc5.cdl
diff tst_nc5.cdl tmp_nc5.cdl

# Verify chunking
${NCDUMP} -hs -n tst_nc5 tmp_nc5.nc > tmp_nc5.cdl
# extract the chunking line
TESTLINE=`sed -e '/ivar:_ChunkSizes/p' -e d <tmp_nc5.cdl`
# track line to match
BASELINE='ivar:_ChunkSizes = 4, 1, 2, 1, 5, 2, 3 ;'
verifychunkline "$TESTLINE" "$BASELINE"

fi #T3

if test "x$T4" = x1 ; then

echo "*** Test nccopy -c with unlimited dimension; classic ->enhanced"
reset
./tst_chunking tst_nc5.nc unlimited
${NCDUMP} -n tst_nc5 tst_nc5.nc > tst_nc5.cdl
${NCCOPY} -c ivar:5,3 tst_nc5.nc tmp_nc5.nc
${NCDUMP} -n tst_nc5 tmp_nc5.nc > tmp_nc5.cdl
diff tst_nc5.cdl tmp_nc5.cdl

# Verify chunking
${NCDUMP} -hs -n tst_nc5 tmp_nc5.nc > tmp_nc5.cdl
# extract the chunking line
TESTLINE=`sed -e '/ivar:_ChunkSizes/p' -e d <tmp_nc5.cdl`
# track line to match
BASELINE='   ivar:_ChunkSizes = 5, 3 ;   '
verifychunkline "$TESTLINE" "$BASELINE"

# Make sure that fvar was not chunked
checkfvar tmp_nc5.cdl

fi # T4

if test "x$T5" = x1 ; then

echo "*** Test nccopy -c fvar: to suppress chunking; classic ->enhanced"
reset
./tst_chunking tst_nc5_omit.nc
${NCDUMP} -n tst_nc5_omit tst_nc5_omit.nc > tst_nc5_omit.cdl
${NCCOPY} -c ivar:7,1,2,1,5,1,9 -c fvar: tst_nc5_omit.nc tmp_nc5_omit.nc
${NCDUMP} -n tst_nc5_omit tmp_nc5_omit.nc > tmp_nc5_omit.cdl
diff tst_nc5_omit.cdl tmp_nc5_omit.cdl

# Verify chunking of ivar
${NCDUMP} -hs -n tst_nc5_omit tmp_nc5_omit.nc > tmp_nc5_omit.cdl
# extract the chunking line
TESTLINE=`sed -e '/ivar:_ChunkSizes/p' -e d <tmp_nc5_omit.cdl`
# track line to match
BASELINE='   ivar:_ChunkSizes = 7, 1, 2, 1, 5, 1, 9 ;   '
verifychunkline "$TESTLINE" "$BASELINE"

# Make sure that fvar was not chunked
checkfvar tmp_nc5_omit.cdl

fi # T5

# Cleanup all created files
reset

echo "*** All nccopy tests passed!"
exit 0

