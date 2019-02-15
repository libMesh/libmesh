#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

. ${srcdir}/d4test_common.sh

echo "test_parse.sh:"

cd ${DMRTESTFILES}
F=`ls -1 *.dmr | sed -e 's/[.]dmr//' |tr '\r\n' '  '`
cd $WD

setresultdir results_test_parse

if test "x${RESET}" = x1 ; then rm -fr ${BASELINE}/*.d4p ; fi
for f in $F ; do
    echo "testing: $f"
    if ! ${VG} ${execdir}/test_parse ${DMRTESTFILES}/${f}.dmr > ./results_test_parse/${f}.d4p ; then
	failure "${f}"
    fi
    if test "x${TEST}" = x1 ; then
	if ! diff -wBb ${BASELINE}/${f}.d4p ./results_test_parse/${f}.d4p ; then
	    failure "${f}"
	fi
    elif test "x${DIFF}" = x1 ; then
	echo "diff -wBb ${DMRTESTFILES}/${f}.dmr ./results_test_parse/${f}.d4p"
	rm -f ./tmp
	cat ./results_test_parse/${f}.d4p \
	| sed -e '/<Dimensions>/d' -e '/<Types>'/d -e '/<Variables>'/d -e '/<Groups>'/d \
	| sed -e '/<\/Dimensions>/d' -e '/<\/Types>'/d -e '/<\/Variables>'/d -e '/<\/Groups>'/d  \
	| sed -e '/_edu.ucar.opaque.size/,+2d' \
	| cat > ./tmp
	if ! diff -wBb ${DMRTESTFILES}/${f}.dmr ./tmp ; then
	    failure "${f}" 
	fi
    elif test "x${RESET}" = x1 ; then
	echo "${f}:" 
	cp ./results_test_parse/${f}.d4p ${BASELINE}/${f}.d4p	
    fi
done
rm -rf ./results_test_parse

finish

exit 0


