#!/bin/sh

if test "x$srcdir" = x ; then srcdir=`pwd`; fi
. ../test_common.sh

# Enable if using localhost
LOCAL=1

RCEMBED=1
RCLOCAL=1
RCHOME=1
RCENV=1
RCPREC=1

# Not currently testable in netcdf
#RCSPEC=1

#SHOW=1
#DBG=1

# Choose at most 1
#GDB=1
#VG=1

NFL=1

WD=`pwd`

NETRCFILE=$WD/test_auth_netrc
# This is the control variable; set when needed
unset NETRC

COOKIES="${WD}/test_auth_cookies"

RC=.daprc

OCLOGFILE=stderr
if test "x$DBG" = x1 ; then
SHOW=1
fi

# Major parameters

BASICCOMBO="tiggeUser:tigge"
BADCOMBO="tiggeUser:xxxxx"
URLPATH="thredds/dodsC/testRestrictedDataset/testData2.nc"
PROTO=http
if test "x$LOCAL" = x ; then
URLSERVER="remotetest.unidata.ucar.edu"
else
URLSERVER="localhost:8081"
fi

# See if we need to override
if test "x$URS" != "x" ; then
#https://54.86.135.31/opendap/data/nc/fnoc1.nc.dds
URLSERVER="54.86.135.31"
URLPATH="opendap/data/nc/fnoc1.nc"
BASICCOMBO="$URS"
RCEMBED=0
NETRC=$NETRCFILE
PROTO=https
fi

if test "x$DBG" = x1 ; then
URLPATH="${URLPATH}#log&show=fetch"
fi

# Split the combo
BASICUSER=`echo $BASICCOMBO | cut -d: -f1`
BASICPWD=`echo $BASICCOMBO | cut -d: -f2`

OUTPUT="./.output"

if test "x$TEMP" = x ; then
  TEMP="/tmp"
fi
TEMP=`echo "$TEMP" | sed -e "s|/$||"`

LOCALRC=./$RC
HOMERC=${HOME}/$RC
HOMERC=`echo "$HOMERC" | sed -e "s|//|/|g"`
SPECRC="$TEMP/temprc"
ENVRC="$WD/envrc"

createrc() {
  RCP="$1" ; shift
  unset NOPWD
  unset BADPWD
  while [[ $# > 0 ]] ; do
    case "$1" in
    nopwd) NOPWD=1 ;;
    badpwd) BADPWD=1 ;;
    *) ;;
    esac
    shift
  done
  if test "x$RCP" != x ; then
    rm -f $RCP
    echo "Creating rc file $RCP"
  else
    echo "createrc: no rc specified"
    exit 1
  fi
  if test "x${DBG}" != x ; then
    echo "HTTP.VERBOSE=1" >>$RCP
  fi	
  echo "HTTP.COOKIEJAR=${COOKIES}" >>$RCP
  if test "x${URS}" = x ; then
    if test "x${NOPWD}" = x ; then
      if test "x${BADPWD}" = x ; then
        echo "HTTP.CREDENTIALS.USERPASSWORD=${BASICCOMBO}" >>$RCP
      else
        echo "HTTP.CREDENTIALS.USERPASSWORD=${BADCOMBO}" >>$RCP
      fi
    fi
  fi
  if test "x${NETRC}" != x && test "x$NFL" = x ; then
    echo "HTTP.NETRC=${NETRC}" >>$RCP
  fi
}

createnetrc() {
  NCP="$1" ; shift
  unset NOPWD
  unset BADPWD
  while [[ $# > 0 ]] ; do
    case "$1" in
    nopwd) NOPWD=1 ;;
    badpwd) BADPWD=1 ;;
    *) ;;
    esac
    shift
  done
  if test "x$NCP" != x ; then
    rm -f $NCP
    echo "Creating netrc file $NCP"
  else
    echo "createnetrc: no rc specified"
    exit 1
  fi
  if test "x$URS" != x ; then
    echo "machine uat.urs.earthdata.nasa.gov login $BASICUSER password $BASICPWD" >>$NCP
    #echo "machine 54.86.135.31 login $BASICUSER password $BASICPWD" >>$1
  else
    echo -n "${PROTO}://$URLSERVER/$URLPATH" >>$NCP
    if test "x$NOPWD" = x ; then
      if test "x$BADPWD" = x ; then
        echo -n " login $BASICUSER password $BASICPWD" >>$NCP
      else
        echo -n " login $BASICUSER password xxxxxx" >>$NCP
      fi
    fi
    echo "" >>$NCP
  fi
}

reset() {
  for f in ./$RC $HOME/$RC $SPECRC $ENVRC $COOKIES $NETRC $OUTPUT ; do
    rm -f ${f}
  done      
  unset DAPRCFILE
}

restore() {
  reset
  for f in ./$RC $HOME/$RC $SPECRC $ENVRC $COOKIES $NETRC ; do
    if test -f ${f}.save ; then
      echo "restoring old ${f}"
      cp ${f}.save ${f}
    fi      
  done      
}

save() {
  for f in ./$RC $HOME/$RC $SPECRC $ENVRC $COOKIES $NETRC ; do
    if test -f $f ; then
      if test -f ${f}.save ; then
        ignore=1
      else
        echo "saving $f"
        cp ${f} ${f}.save
      fi
    fi      
  done      
}

show() {
  if test "x$SHOW" = x1 ; then cat $OUTPUT; fi
  if test "x$OUTPUT" != "x"; then rm -f $OUTPUT; fi
}

# Assemble the ncdump command
if test "x$DBG" = x1; then
NCDUMP="$NCDUMP -D1"
fi

if test "x$GDB" = x1 ; then
  NCDUMP="gdb --args $NCDUMP"
fi
if test "x$VG" = x1 ; then
NCDUMP="valgrind --leak-check=full $NCDUMP"
fi

# Initialize
save
reset

if test "x$RCEMBED" = x1 ; then
  echo "***Testing rc file with embedded user:pwd"
  URL="${PROTO}://${BASICCOMBO}@${URLSERVER}/$URLPATH"
  unset NETRC
  # Invoke ncdump to extract a file the URL
  echo "command: ${NCDUMP} -h ${URL} > $OUTPUT"
  ${NCDUMP} -h "$URL" > $OUTPUT
  show
fi

# Rest of tests assume these defaults
URL="${PROTO}://${URLSERVER}/$URLPATH"
NETRC=$NETRCFILE

if test "x$RCLOCAL" = x1 ; then
  echo "***Testing rc file in local directory"
  # Create the rc file and (optional) netrc fil in ./
  reset
  createnetrc $NETRC
  createrc $LOCALRC

  # Invoke ncdump to extract a file using the URL
  echo "command: ${NCDUMP} -h ${URL} > $OUTPUT"
  ${NCDUMP} -h "$URL" > $OUTPUT
  show
fi

if test "x$RCHOME" = x1 ; then
  echo "***Testing rc file in home directory"
  # Create the rc file and (optional) netrc file in ./
  reset
  createnetrc $NETRC
  createrc $HOMERC

  # Invoke ncdump to extract a file the URL
  echo "command: ${NCDUMP} -h ${URL} > $OUTPUT"
  ${NCDUMP} -h "$URL" > $OUTPUT
  show
fi

if test "x$RCSPEC" == x1 ; then
  echo "*** Testing rc file in specified directory"
  # Create the rc file and (optional) netrc file
  reset
  createnetrc $NETRC
  createrc $SPECRC

  # Invoke ncdump to extract a file the URL
  echo "command: ${NCDUMP} -h ${URL} > $OUTPUT"
  ${NCDUMP} -h "$URL" > $OUTPUT
  show
fi

if test "x$RCENV" = x1 ; then
  echo "*** Testing rc file using env variable"
  # Create the rc file and (optional) netrc file
  reset
  createnetrc $NETRC
  echo "ENV: export DAPRCFILE=$ENVRC"
  export DAPRCFILE=$ENVRC
  createrc $DAPRCFILE

  # Invoke ncdump to extract a file the URL
  echo "command: ${NCDUMP} -h ${URL} > $OUTPUT"
  ${NCDUMP} -h "$URL" > $OUTPUT
  show
  export DAPRCFILE=
fi

# Test that .daprc overrides netcrc for password
URL="${PROTO}://${URLSERVER}/$URLPATH"
NETRC=$NETRCFILE
if test "x$RCPREC" = x1 ; then
  echo "***Testing rc vs netrc file precedence"
  # Create the rc file and (optional) netrc file in ./
  reset
  createnetrc $NETRC badpwd
  createrc $LOCALRC

  # Invoke ncdump to extract a file using the URL
  echo "command: ${NCDUMP} -h ${URL} > $OUTPUT"
  ${NCDUMP} -h "$URL" > $OUTPUT
  ${NCDUMP} -h "$URL"
  show
fi

reset
restore

exit

