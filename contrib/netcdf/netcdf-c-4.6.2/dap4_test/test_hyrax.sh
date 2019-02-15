#!/bin/sh

if test "x$srcdir" = "x"; then srcdir=`dirname $0`; fi
export srcdir;

. ${srcdir}/../test_common.sh

FRAG="#checksummode=ignore"

F="\
nc4_nc_classic_comp.nc \
nc4_nc_classic_no_comp.nc \
nc4_strings.nc \
nc4_strings_comp.nc \
nc4_unsigned_types.nc \
nc4_unsigned_types_comp.nc \
ref_tst_compounds.nc \
"

failure() {
      echo "*** Fail: $1"
      exit 1
}

setresultdir results_test_hyrax

if test "x${RESET}" = x1 ; then rm -fr ${BASELINEH}/*.dmp ; fi
for f in $F ; do
    URL="dap4://test.opendap.org:8080/opendap/nc4_test_files/${f}${FRAG}"
    echo "testing: $URL"
    if ! ${NCDUMP} "${URL}" > ./results_test_hyrax/${f}.hyrax; then
        failure "${URL}"
    fi
    if test "x${TEST}" = x1 ; then
	if ! diff -wBb ${BASELINEREM}/${f}.hyrax ./results_test_hyrax/${f}.hyrax ; then
	    failure "diff ${f}.hyrax"
	fi
    elif test "x${RESET}" = x1 ; then
	echo "${f}:" 
	cp ./results_test_hyrax/${f}.hyrax ${BASELINEH}/${f}.hyrax
    fi
done

rm -rf ./results_test_hyrax

echo "*** Pass"
exit 0

