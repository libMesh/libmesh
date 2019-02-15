#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

set -e
. ${srcdir}/d4test_common.sh

echo "test_raw.sh:"

# Compute the set of testfiles
cd ${srcdir}/daptestfiles
F=`ls -1d *.dap`
cd -
F=`echo $F | tr '\r\n' '  '`
F=`echo $F | sed -e s/.dap//g`

# Do cleanup on the baseline file
baseclean() {
  if test $# != 2 ; then
    echo "simplify: too few args"
  else
    rm -f $2
    while read -r iline; do
      oline=`echo $iline | tr "'" '"'`
      echo "$oline" >> $2
    done < $1
  fi
}

# Do cleanup on the result file
resultclean() {
  if test $# != 2 ; then
    echo "simplify: too few args"
  else
    rm -f $2
    while read -r iline; do
      oline=`echo $iline | sed -e 's|^\(netcdf.*\)[.]nc\(.*\)$|\\1\\2|'`
      echo "$oline" >> $2
    done < $1
  fi
}

setresultdir results_test_raw

if test "x${RESET}" = x1 ; then rm -fr ${BASELINERAW}/*.dmp ; fi
for f in $F ; do
    echo "testing: $f"
    URL="[log][dap4]file://${DAPTESTFILES}/${f}"
    if ! ${NCDUMP} "${URL}" > ${builddir}/results_test_raw/${f}.dmp; then
        failure "${URL}"
    fi
    if test "x${TEST}" = x1 ; then
	if ! diff -wBb ${BASELINERAW}/${f}.dmp ${builddir}/results_test_raw/${f}.dmp ; then
	    failure "diff ${f}.dmp"
	fi
    elif test "x${RESET}" = x1 ; then
	echo "${f}:"
	cp ${builddir}/results_test_raw/${f}.dmp ${BASELINERAW}/${f}.dmp
    elif test "x${DIFF}" = x1 ; then
	echo "hdrtest: ${f}"
	baseclean
        if ! diff -wBb ${BASELINERAW}/${f}.dmp ${BASELINE}/${f}.ncdump ; then
          failure diff -wBb ${BASELINERAW}/${f}.dmp ${BASELINE}/${f}.ncdump
	fi
    fi
done
rm -rf ${builddir}/results_test_raw

finish
