#!/bin/sh
# Run test_chunk_hdf4 passing ${src_dir}
set -x
CHUNKED=chunked.hdf4
CONTIG=contiguous.hdf4

echo ""
echo "*** Testing hdf4 chunking..."

if test "x${src_dir}" = "x" ; then
src_dir="."
fi

# Move the data sets into place
ISDISTCHECK=0
if test -f ./${CHUNKED} ; then
ISDISTCHECK=0
else
ISDISTCHECK=1
cp ${src_dir}/${CHUNKED} .
cp ${src_dir}/${CONTIG} .
fi

if ./tst_chunk_hdf4 ; then
  echo "***SUCCESS!! tst_chunk_hdf4"
else
  echo "***FAIL: tst_chunk_hdf4"
fi

if test "x${ISDISTCHECK}" = "x1" ; then
echo rm -f ./${CHUNKED} ./${CONTIG}
fi

exit 0
