##################################################
# Remote test info
##################################################

# This is not a standalone script; it is invoked by
# one or more of the main scripts

# Define various sets of test targets

# Figure our dst server; if none, then just stop
DTS=`${execdir}/findtestserver dap2 dts`
if test "x$DTS" = "x" ; then
echo "WARNING: Cannot locate test server for dts"
exit 1
fi

if test "x$timing" = "x1" ; then TIMECMD="time"; else TIMECMD=""; fi
expected3="${srcdir}/expectremote3"

# For special testing
X="test.03"
XC="test.03;1;s0,s1"

# Set to run a specific test
S0="test.07"

# These shorter tests are always run
S1="\
test.01 test.02 test.04 test.05 test.07a test.07 \
test.21 \
test.50 test.53 test.55 test.56 test.57 \
test.66 test.67 test.68 test.69"

# Server is failing on some tests ; investigate why
S1FAIL="test.06a test.22 test.23 test.31"

# This fails because of the duplicate name problem and the way thredds dap handles grids"
S2FAIL="test.06"

# These longer tests are optional
SL1="\
test.03 \
b31 b31a D1 Drifters EOSDB \
ingrid nestedDAS NestedSeq NestedSeq2 OverideExample \
SimpleDrdsExample test.an1 \
test.dfp1 test.gr1 \
test.gr2 test.gr3 test.gr4 test.gr5 \
test.sds1 test.sds2 test.sds3 \
test.sds4 test.sds5 \
test.vs1 test.vs2 test.vs3 test.vs4 test.vs5 \
whoi"

# Anything larger than about 100k will not be in the distribution
TOOBIGL1="parserBug0001 test.satimage Sat_Images test.32"

# Following contain %XX escapes which I cannot handle yet
ESCAPEDFAIL="test.dfr1 test.dfr2 test.dfr3 test.GridFile test.PointFile test.SwathFile test.sds6 test.sds7"

# Following tests are to check constraint handling
C1="\
test.01;1;f64 \
test.02;1;b[1:2:10] \
test.03;1;i32[0:1][1:2][0:2] \
test.04;1;types.i32 \
test.05;1;types.floats.f32 \
test.07;1;person.age \
test.07;3;person \
test.07;4;types.f32"

# See S2FAIL above
SR1FAIL="test.06;1;ThreeD"

# Constrained long tests
LC1="test.03;2;s1"

# Unknown problem: test.07;2;&age>2
IGNORE="test.07.2" 

# Columbia hack test Not sure if operative
COLUMBIA="http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NASA-GMAO/.MONTHLY/.sst"

# Known to fail

XFAILTESTS=
# Suppress some tests if not windows platform.
if test "x$platform" == xmingw ; then
    XFAILTESTS="$XFAILTESTS test.67"
fi

# Following tests must be run as not cached
NOCACHETESTS="test.07"

# Following tests must be run as not prefetch
NOPREFETCHTESTS="test.07"

computewhich() { # set REMOTETESTS and constrained
  case "$1" in
  S0) REMOTETESTS="$S0" ; constrained=0 ;;
  S1) REMOTETESTS="$S1" ; constrained=0 ;;
  S2) REMOTETESTS="$S2" ; constrained=0 ;;
  L1) REMOTETESTS="$L1" ; constrained=0 ;;
  L2) REMOTETESTS="$L2" ; constrained=0 ;;
  C1) REMOTETESTS="$C1" ; constrained=1 ;;
  C2) REMOTETESTS="$C2" ; constrained=1 ;;
  C3) REMOTETESTS="$C3" ; constrained=1 ;;
  LC1) REMOTETESTS="$LC1" ; constrained=1 ;;
  LC2) REMOTETESTS="$LC2" ; constrained=1 ;;
  X) REMOTETESTS="$X" ; constrained=0 ;;
  XC) REMOTETESTS="$XC" ; constrained=1 ;;
  *) echo "Unknown which test: $1" ; exit 1;;
  esac
}

constrain() {
  T="$1;;" # add semicolons to fake out the cut command
  # see if we are using constraints will set testname and ce and testno and contrained
  testname=`echo -n $T | cut "-d;" -f1`
  testno=`echo -n $T | cut "-d;" -f2`
  ce=`echo -n $T | cut "-d;" -f3`
  if test "x$ce" = x ; then constrained=0; else constrained=1; fi
}

setcache() {
  CACHE="[cache]"
  if test "x${NOCACHETESTS}" != x ; then
    if IGNORE=`echo -n " ${NOCACHETESTS} " | fgrep " $1 "`; then CACHE=; fi
  fi
  PARAMS="${PARAMS}${CACHE}"
}

setprefetch() {
  PREFETCH="[prefetch]"
  if test "x${NOPREFETCHTESTS}" != x ; then
    if IGNORE=`echo -n " ${NOPREFETCHTESTS} " | fgrep " $1 "`; then PREFETCH=; fi
  fi
  PARAMS="${PARAMS}${PREFETCH}"
}

# Use specialized test executor
doremotetests() {
for x in ${REMOTETESTS} ; do
  PARAMS=
  setcache $x
  setprefetch $x
  constrain $x
  if test "x$constrained" = "x1" ; then
    name="${testname}.${testno}"
    url="${PARAMS}${DTS}/$testname?${ce}"
  else
    name="${testname}"
    url="${PARAMS}${DTS}/$testname"
  fi
  if test "x$quiet" = "x0" ; then echo "*** Testing: ${name} ; url=$url" ; fi
  # determine if this is an xfailtest
  isxfail=0
  if test "x${XFAILTESTS}" != x ; then
    if IGNORE=`echo -n " ${XFAILTESTS} " | fgrep " ${name} "`; then isxfail=1; fi
  fi
  ok=1
  if ${NCDUMP} ${FLAGS} "${url}" | sed 's/\\r//g' > ${name}.dmp ; then ok=$ok; else ok=0; fi
  # compare with expected
  if diff -w ${EXPECTED}/${name}.dmp ${name}.dmp  ; then ok=$ok; else ok=0; fi
   processstatus
done
}
