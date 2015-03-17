#!/bin/sh
#set -x

#NOEMBED=1
#NOLOCAL=1
#NOHOME=1

SHOW=1
#DBG=1
#GDB=1

CMP=1

WD=`pwd`

# if this is part of a distcheck action, then this script
# will be executed in a different directory
# than the one containing it; so capture the path to this script
# as the location of the source directory.

# capture the build directory
builddir=`pwd`/..
if test "x$TOPSRCDIR" != x ; then
srcdir="$TOPSRCDIR/ncdap_test"
else
srcdir=`dirname $0`
fi
# canonical
cd $srcdir
srcdir=`pwd`

# Hack for CYGWIN
if [ `uname | cut -d "_" -f 1` = "MINGW32" ]; then
    srcdir=`echo $srcdir | sed 's/\/c\//c:\//g'`
    builddir=`echo $builddir | sed 's/\/c\//c:\//g'`
fi
cd ${builddir}/ncdap_test

NETRCFILE=$WD/test_auth_netrc
# This is the control variable
NETRC=$NETRCFILE
COOKIES="${WD}/test_auth_cookies"
RC=.ocrc
EXPECT=$srcdir/expected3/testData.nc.dmp
OCLOGFILE=stderr


# Major parameters

BASICCOMBO="tiggeUser:tigge"
URLSERVER="remotetest.unidata.ucar.edu"
URLPATH="thredds/dodsC/restrict/testData.nc"

# See if we need to override
if test "x$URS" != "x" ; then
#https://54.86.135.31/opendap/data/nc/fnoc1.nc.dds
URLSERVER="54.86.135.31"
URLPATH="opendap/data/nc/fnoc1.nc"
BASICCOMBO="$URS"
NOEMBED=1
NETRC=$NETRCFILE
CMP=
else
NETRC=
fi

if test "x$CMP" = x1 ; then SHOW= ; else SHOW=1; fi
if test "x$DBG" = x1 ; then SHOW=1; fi

# Split the combo
BASICUSER=`echo $BASICCOMBO | cut -d: -f1`
BASICPWD=`echo $BASICCOMBO | cut -d: -f2`

NCDUMP=$builddir/ncdump/ncdump

LOCALRC=./$RC
HOMERC=${HOME}/$RC
HOMERC=`echo "$HOMERC" | sed -e "s|//|/|g"`
SPECRC=$TEMP
if test "x$TEMP" = x ; then
  SPECRC="/tmp"
fi
SPECRC="$SPECRC/temprc"

function createrc {
  if test "x$1" = x ; then
    echo "createrc: no rc file specified"
  else
  rm -f $1
  RCP=$1
  NRC=$2
  echo "Creating rc file $RCP"
  if test "x${DBG}" != x ; then
    echo "HTTP.VERBOSE=1" >>$RCP
  fi	
  echo "HTTP.COOKIEJAR=${COOKIES}" >>$RCP
  if test "x${URS}" = x ; then
    echo "HTTP.CREDENTIALS.USERPASSWORD=${BASICCOMBO}" >>$RCP
  fi
  if test "x${NRC}" != x ; then
    echo "HTTP.NETRC=${NRC}" >>$RCP
  fi
  fi
}

function createnetrc {
  if test "x$1" != x ; then
    rm -f $1
    echo "Creating netrc file $1"
    echo "machine uat.urs.earthdata.nasa.gov login $BASICUSER password $BASICPWD" >>$1
#   echo "machine 54.86.135.31 login $BASICUSER password $BASICPWD" >>$1
  fi
}

function reset {
  for f in ./$RC $HOME/$RC $SPECRC $COOKIES $NETRC ; do
    if test "x$DBG" = x1 ; then echo "Deleting $f"; fi
    rm -f ${f}
    if test -f ${f}.save ; then
        if test "x$DBG" = x1 ; then echo "restoring old ${f}" ; fi
      cp ${f}.save ${f}
    fi      
  done      
  # unconditional
  rm -f ./tmp
}

function save {
  for f in ./$RC $HOME/$RC $SPECRC $COOKIES $NETRC ; do
    if test -f $f ; then
      if test -f ${f}.save ; then
        ignore=1
      else
        if test "x$DBG" = x1 ; then echo "saving $f"; fi
        cp ${f} ${f}.save
      fi
    fi      
  done      
}

function compare {
  if test "x$CMP" = x1 ; then
    if diff -w ./tmp $EXPECT ; then
      echo "***Pass"
    else
      echo "***FAIL"
      pass=1
    fi
  fi
}

# Invoke ncdump to extract the URL
function doit {
  reset
  createnetrc $2
  createrc $1 $2
  echo "command: ${NCDUMP} -h $URL"
  ${NCDUMP} -h "$URL" >./tmp
  if test "x$CMP" = x1 ; then compare ; fi
  if test "x${SHOW}" = x1 ; then cat ./tmp ; fi
}

pass=0

# Assemble the ncdump command
if test "x$GDB" = x1 ; then
  NCDUMP="gdb --args $NCDUMP"
fi

# Initialize
save
reset

if test "x$NOEMBED" != x1 ; then
echo "***Testing rc file with embedded user:pwd"
URL="https://${BASICCOMBO}@${URLSERVER}/$URLPATH"
# Invoke ncdump to extract the URL
doit
fi

URL="https://${URLSERVER}/$URLPATH"
if test "x$NOLOCAL" != x1 ; then
  URL="https://${URLSERVER}/$URLPATH"
  echo "***Testing rc file in local directory"
  doit $LOCALRC $NETRC
fi

if test "x$NOHOME" != x1 ; then
  URL="https://${URLSERVER}/$URLPATH"
  echo "***Testing rc file in home directory"
  doit $HOMERC $NETRC
fi

#cleanup
reset

exit $pass
