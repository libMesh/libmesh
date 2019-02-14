#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

. ${srcdir}/d4test_common.sh

echo "test_meta.sh:"

set -e

cd ${DMRTESTFILES}
F=`ls -1 *.dmr | sed -e 's/[.]dmr//g' | tr '\r\n' '  '`
cd $WD

CDL=
for f in ${F} ; do
STEM=`echo $f | cut -d. -f 1`
if test -e ${CDLTESTFILES}/${STEM}.cdl ; then
  CDL="${CDL} ${STEM}"
else
  echo "Not found: ${CDLTESTFILES}/${STEM}.cdl; ignored"
fi
done

if test "x${RESET}" = x1 ; then rm -fr ${BASELINE}/*.d4m ; fi

setresultdir results_test_meta

for f in ${F} ; do
    echo "checking: $f"
    if ! ${VG} ${execdir}/test_meta ${DMRTESTFILES}/${f}.dmr ./results_test_meta/${f} ; then
        failure "${execdir}/test_meta ${DMRTESTFILES}/${f}.dmr ./results_test_meta/${f}"
    fi
    ${NCDUMP} ${DUMPFLAGS} -h ./results_test_meta/${f} > ./results_test_meta/${f}.d4m
    if test "x${TEST}" = x1 ; then
	if ! diff -wBb ${BASELINE}/${f}.d4m ./results_test_meta/${f}.d4m ; then
	    failure "diff -wBb ${BASELINE}/${f}.ncdump ./results_test_meta/${f}.d4m"
	fi
    elif test "x${RESET}" = x1 ; then
	echo "${f}:" 
	cp ./results_test_meta/${f}.d4m ${BASELINE}/${f}.d4m
    fi
done

if test "x${CDLDIFF}" = x1 ; then
  for f in $CDL ; do
    echo "diff -wBb ${CDLTESTFILES}/${f}.cdl ./results_test_meta/${f}.d4m"
    rm -f ./tmp
    cat ${CDLTESTFILES}/${f}.cdl \
    cat >./tmp
    echo diff -wBbu ./tmp ./results_test_meta/${f}.d4m
    if ! diff -wBbu ./tmp ./results_test_meta/${f}.d4m ; then
	failure "${f}" 
    fi
  done
fi
rm -rf ./results_test_meta
finish


