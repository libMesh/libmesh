#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

. ${srcdir}/d4test_common.sh

echo "test_data.sh:"

cd ${DAPTESTFILES}
F=`ls -1 *.dap | sed -e 's/[.]dap//g' | tr '\r\n' '  '`
cd $WD

setresultdir results_test_data

if test "x${RESET}" = x1 ; then rm -fr ${BASELINE}/*.d4d ; fi
for f in $F ; do
    echo "testing: ${f}"
    if ! ${VG} ${execdir}/test_data ${DAPTESTFILES}/${f} ./results_test_data/${f}.nc ; then
        failure "${execdir}/test_data ${DAPTESTFILES}/${f} ./results_test_data/${f}.nc"
    fi
    ${NCDUMP} ./results_test_data/${f}.nc > ./results_test_data/${f}.d4d
    if test "x${TEST}" = x1 ; then
	if ! diff -wBb ${BASELINE}/${f}.d4d ./results_test_data/${f}.d4d ; then
	    failure "diff -wBb ${BASELINE}/${f}.d4d ./results_test_data/${f}.d4d"
	fi
    elif test "x${RESET}" = x1 ; then
	echo "${f}:"
	cp ./results_test_data/${f}.d4d ${BASELINE}/${f}.d4d
    fi
done

# Remove empty lines and trim lines in a cdl file
trim() {
  if test $# != 2 ; then
    echo "simplify: too few args"
  else
    rm -f $2
    while read -r iline; do
      oline=`echo $iline | sed -e 's/^[\t ]*\([^\t ]*\)[\t ]*$/\\1/'`
      if test "x$oline" = x ; then continue ; fi
      echo "$oline" >> $2
    done < $1
  fi
}

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

if test "x${CDLDIFF}" = x1 ; then
  for f in $F ; do
    STEM=`echo $f | cut -d. -f 1`
    if ! test -e ${CDLTESTFILES}/${STEM}.cdl ; then
      echo "Not found: ${CDLTESTFILES}/${STEM}.cdl"
      continue
    fi
    echo "diff -wBb ${CDLTESTFILES}/${STEM}.cdl ./results_test_data/${f}.d4d"
    rm -f ./b1 ./b2 ./r1 ./r2
    trim ${CDLTESTFILES}/${STEM}.cdl ./b1
    trim ./results_test_data/${f}.d4d ./r1
    baseclean b1 b2
    resultclean r1 r2
    if ! diff -wBb ./b2 ./r2 ; then
	failure "${f}"
    fi
  done
fi
rm -rf ./results_test_data

finish
