#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

. ${srcdir}/d4test_common.sh

set -e

echo "test_remote.sh:"

#BIG=1
#NOCSUM=1

F="\
test_atomic_array.nc \
test_atomic_types.nc \
test_enum.nc \
test_enum_2.nc \
test_enum_array.nc \
test_fill.nc \
test_groups1.nc \
test_misc1.nc \
test_one_var.nc \
test_one_vararray.nc \
test_opaque.nc \
test_opaque_array.nc \
test_struct1.nc \
test_struct_array.nc \
test_struct_nested.nc \
test_struct_nested3.nc \
test_struct_type.nc \
test_utf8.nc \
test_vlen1.nc \
test_vlen2.nc \
test_vlen3.nc \
test_vlen4.nc \
test_vlen5.nc \
test_vlen6.nc \
test_vlen7.nc \
test_vlen8.nc \
test_vlen9.nc \
test_vlen10.nc \
test_vlen11.nc \
tst_fills.nc \
test_struct_nested.hdf5 \
test_struct_nested3.hdf5 \
test_vlen3.hdf5 \
test_vlen4.hdf5 \
test_vlen5.hdf5 \
test_anon_dim.syn \
test_atomic_array.syn \
test_atomic_types.syn \
test_sequence_1.syn \
test_sequence_2.syn \
test_struct_array.syn \
"

setresultdir results_test_remote

TESTSERVER=`${execdir}/findtestserver4 dap4 d4ts`
if test "x$TESTSERVER" = x ; then
echo "***XFAIL: Cannot find d4ts testserver"
exit 1
fi

if test "x${RESET}" = x1 ; then rm -fr ${BASELINER}/*.dmp ; fi
for f in $F ; do
    URL="[log][show=fetch][dap4]${TESTSERVER}/testfiles/${f}"
    if test "x$BIG" = x1; then
	URL="[ucar.littleendian=0]${URL}"
    fi
    if test "x$NOCSUM" = x1; then
	URL="[ucar.checksummode=none]${URL}"
    fi
    if ! ${NCDUMP} ${DUMPFLAGS} "${URL}" > ${builddir}/results_test_remote/${f}.dmp; then
        failure "${URL}"
    fi
    if test "x${TEST}" = x1 ; then
	if ! diff -wBb "${BASELINEREM}/${f}.dmp" "${builddir}/results_test_remote/${f}.dmp" ; then
	    failure "diff ${f}.dmp"
	fi
    elif test "x${RESET}" = x1 ; then
	echo "${f}:" 
	cp "${builddir}/results_test_remote/${f}.dmp" "${BASELINEREM}/${f}.dmp"
    fi
done

rm -fr ${builddir}/results_test_remote

finish

