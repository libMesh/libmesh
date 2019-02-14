if test $# = 0 ; then
TEST=1
else
for arg in "$@"; do
  case "${arg}" in
  test) TEST=1 ;;
  reset) RESET=1 ;;
  diff) DIFF=1 ;;
  cdl) CDLDIFF=1 ;;
  log) ;;
  *) echo unknown argument $arg ;;
  esac
done
fi

# Define input paths
WD=`pwd`
cd ${srcdir}/daptestfiles; DAPTESTFILES=`pwd` ; cd ${WD}
cd ${srcdir}/dmrtestfiles; DMRTESTFILES=`pwd` ; cd ${WD}
cd ${srcdir}/cdltestfiles; CDLTESTFILES=`pwd` ; cd ${WD}
cd ${srcdir}/baseline; BASELINE=`pwd` ; cd ${WD}
cd ${srcdir}/baselineraw; BASELINERAW=`pwd` ; cd ${WD}
cd ${srcdir}/baselineremote; BASELINEREM=`pwd` ; cd ${WD}
BASELINEH=${BASELINEREM}

setresultdir() {
rm -fr ${builddir}/$1
mkdir -p ${builddir}/$1
}

FAILURES=
failure() {
  echo "*** Fail: $1"
  FAILURES=1
  if test "x$2" != x ; then
    exit 1
  fi
}

PUSHD() {
pushd $1 >>/dev/null
}
POPD() {
popd >>/dev/null
}

filesexist() {
    for x in "$@" ; do
	if ! test -e $x ; then
	  failure "missing file: $x"
	fi
    done
}

finish() {
if test "x$FAILURES" = x1 ; then
echo "*** Fail"
exit 1
else
echo "*** Pass"
exit 0
fi
}

# Return the complete set of arguments minus any found in $SUPPRESS
suppress() {
  local args="$@"
  suppress=
  for a in $args; do
  for s in $SUPPRESS; do if test "x$a" != "x$s" ; then suppress="$suppress $a"; fi;  done
  done        
}

VG="valgrind --leak-check=full --error-exitcode=1 --num-callers=100"
if test "x$USEVG" = x ; then VG=; fi

DUMPFLAGS=

