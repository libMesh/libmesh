dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl SLEPc
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_SLEPC],
[
  dnl Although PETSc 3.x does not require PETSC_ARCH, it appears that
  dnl SLEPc 3.x may (based on patches sent in by Jed Brown)?  
  dnl If that's the case, we should probably uncomment 
  dnl this check and disable SLEPc when PETSC_ARCH is not set...
  dnl if test "x$PETSC_ARCH" = x ; then
  dnl   enableslepc=no
  dnl   AC_MSG_RESULT([<<< SLEPc disabled.  Please set your "\$PETSC_ARCH" environment variable correctly. >>>])
  dnl fi

  dnl Test to see if SLEPC_DIR set by user.  If not set, then
  dnl try to autodetect in a default directory
dnl  if test "x$SLEPC_DIR" = x ; then
dnl    SLEPC_DIR=/usr/lib/slepc
dnl  fi

  dnl Test to see if SLEPC_DIR and petscversion were set by user or
  dnl autodetection.  If not set, then disable slepc, print a message.
  if test "x$SLEPC_DIR" = x || test "x$petscversion" = x; then
    enableslepc=no
    AC_MSG_RESULT(<<< SLEPc disabled.  Please set your "\$SLEPC_DIR" environment variable correctly to enable SLEPc. >>>)

  else
    AC_CHECK_FILE($SLEPC_DIR/include/slepc.h,
                  SLEPC_H_PATH=$SLEPC_DIR/include/slepc.h)
                                                                                         
    if (test -r $SLEPC_DIR/include/slepc.h) ; then
  
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
      fi
dnl        AC_MSG_RESULT(Will skip SLEPc support)
dnl        enableslepc=no
dnl      else	
        AC_DEFINE(HAVE_SLEPC, 1,
                  [Flag indicating whether or not SLEPc is available])
        AC_MSG_RESULT(<<< Configuring library with SLEPc version $slepcversion support >>>)
dnl      fi
    else
      AC_MSG_RESULT(<<< Invalid "\$SLEPC_DIR" detected (slepc.h not found). SLEPc disabled. >>>)
      enableslepc=no
    fi
  fi

  AC_SUBST(SLEPC_DIR)
  AC_SUBST(enableslepc)
])
