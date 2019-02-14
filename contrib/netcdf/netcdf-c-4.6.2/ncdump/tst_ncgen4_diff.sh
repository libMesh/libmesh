#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi 
. ../test_common.sh


. ${srcdir}/tst_ncgen_shared.sh

ALLXFAIL="${XFAILTESTS} ${SPECIALTESTS}"

if test "${MODE}" = 3 ; then
  if test "${KFLAG}" = 4 ; then
    TESTSET="${TESTS3} ${SPECIALTESTS3}"
  else
    TESTSET="${TESTS3}"
  fi
else
TESTSET="${TESTS3} ${TESTS4} ${ALLXFAIL}"
fi

echo "*** Testing ncgen with -k${KFLAG}"

mkdir ${RESULTSDIR}
cd ${RESULTSDIR}
for x in ${TESTSET} ; do
  if test $verbose = 1 ; then echo "*** Testing: ${x}" ; fi
	# determine if we need the specflag set
	specflag=
	headflag=
	for s in ${SPECIALTESTS} ; do
	if test "x${s}" = "x${x}" ; then specflag="-s"; headflag="-h"; fi
	done
	# determine if this is an xfailtest
	isxfail=
	for t in ${ALLXFAIL} ; do
	if test "x${t}" = "x${x}" ; then isxfail=1; fi
	done
  rm -f ${x}.nc ${x}.dmp
  ${NCGEN} -b -k${KFLAG} -o ${x}_$$.nc ${cdl}/${x}.cdl
  # dump .nc file
  # if windows, we need to remove any leading 0's in exponents.
  ${NCDUMP} ${headflag} ${specflag} -n ${x} ${x}_$$.nc | sed 's/e+0/e+/g' > ${x}.dmp
  # compare the expected (silently if XFAIL)
  if test "x$isxfail" = "x1" -a "x$SHOWXFAILS" = "x" ; then
    if diff -b -bw ${expected}/${x}.dmp ${x}.dmp >/dev/null 2>&1; then ok=1; else ok=0; fi
  else
    if diff -b -w ${expected}/${x}.dmp ${x}.dmp ; then ok=1; else ok=0; fi
  fi
  if test "x$ok" = "x1" ; then
    test $verbose = 1 && echo "*** SUCCEED: ${x}"
    passcount=`expr $passcount + 1`
  elif test "x${isxfail}" = "x1" ; then
    echo "*** XFAIL : ${x}"
    xfailcount=`expr $xfailcount + 1`
  else
    echo "*** FAIL: ${x}"
    failcount=`expr $failcount + 1`
  fi
done
cd ..

totalcount=`expr $passcount + $failcount + $xfailcount`
okcount=`expr $passcount + $xfailcount`

echo "*** PASSED: ${okcount}/${totalcount} ; ${xfailcount} expected failures ; ${failcount} unexpected failures"
rm -rf ${RESULTSDIR}

if test $failcount -gt 0 ; then
  exit 1
else
  exit 0
fi
