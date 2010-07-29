dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl SLEPc
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_SLEPC],
[
  dnl Test to see if SLEPC_DIR set by user.  If not set, then
  dnl try to autodetect in a default directory
  if test "x$SLEPC_DIR" = x ; then
    SLEPC_DIR=/usr/lib/slepc
  fi

  dnl Test to see if SLEPC_DIR set by user or autodetection.  If not set, then
  dnl disable slepc, print a message.
  if test "x$SLEPC_DIR" = x ; then
    enableslepc=no
    AC_MSG_RESULT(<<< SLEPc disabled.  Please set your "\$SLEPC_DIR" environment variable correctly. >>>)

  else
    AC_CHECK_FILE($SLEPC_DIR/include/slepc.h,
                  SLEPC_H_PATH=$SLEPC_DIR/include/slepc.h)
                                                                                         
    if (test -r $SLEPC_DIR/include/slepc.h) ; then
      AC_SUBST(SLEPC_DIR)
  
      dnl Similar to petsc, we need the slepc version number.
      dnl Note slepc will most likely only work with the corresponding version of petsc
      slepcmajor=`grep "define SLEPC_VERSION_MAJOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_MAJOR[ ]*//g"`
      slepcminor=`grep "define SLEPC_VERSION_MINOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_MINOR[ ]*//g"`
      slepcsubminor=`grep "define SLEPC_VERSION_SUBMINOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_SUBMINOR[ ]*//g"`
      slepcversion=$slepcmajor.$slepcminor.$slepcsubminor
      AC_SUBST(slepcversion)
  
      if (test $slepcversion != $petscversion) ; then
        AC_MSG_RESULT(WARNING:)
        AC_MSG_RESULT(>>> Different version numbers for SLEPc and PETSc <<<)
        AC_MSG_RESULT(Will skip SLEPc support)
        enableslepc=no
      else	
        AC_DEFINE(HAVE_SLEPC, 1,
                  [Flag indicating whether or not SLEPc is available])
        AC_MSG_RESULT(<<< Configuring library with SLEPc version $slepcversion support >>>)
      fi
    else
      enableslepc=no
    fi
  fi

  AC_SUBST(enableslepc)
])
