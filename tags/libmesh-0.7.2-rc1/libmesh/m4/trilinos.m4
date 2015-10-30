dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl Trilinos 10
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRILINOS_10], 
[
  if test "x$TRILINOS_DIR" = "x"; then
    TRILINOS_DIR=no
  fi  	

  AC_ARG_WITH(trilinos,
              AC_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  if test "$withtrilinosdir" != no ; then
    AC_CHECK_FILE($withtrilinosdir/include/Makefile.export.Trilinos,
                  TRILINOS_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.Trilinos,
    AC_CHECK_FILE($withtrilinosdir/Makefile.export.Trilinos,
                  TRILINOS_MAKEFILE_EXPORT=$withtrilinosdir/Makefile.export.Trilinos,
	 	enabletrilinos10=no))

    if test "$enabletrilinos10" != no ; then
       enabletrilinos10=yes
       AC_DEFINE(HAVE_TRILINOS, 1,
                 [Flag indicating whether the library shall be compiled to use the Trilinos solver collection])
       AC_MSG_RESULT(<<< Configuring library with Trilinos 10 support >>>)
       
       dnl ------------------------------------------------------
       dnl AztecOO
       dnl ------------------------------------------------------
       AC_CHECK_FILE($withtrilinosdir/include/AztecOO_config.h,
                     enableaztecoo=yes,
       AC_CHECK_FILE($withtrilinosdir/packages/aztecoo/src/AztecOO_config.h,
                     enableaztecoo=yes,
                     enableaztecoo=no))
                     
       if test "$enableaztecoo" != no ; then
          AC_DEFINE(HAVE_AZTECOO, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos AztecOO linear solver])
          AC_MSG_RESULT(<<< Configuring library with AztecOO support >>>)
       fi
       
       dnl ------------------------------------------------------
       dnl NOX
       dnl ------------------------------------------------------
       AC_CHECK_FILE($withtrilinosdir/include/NOX_Config.h,
                     enablenox=yes,
       AC_CHECK_FILE($withtrilinosdir/packages/nox/src/NOX_Config.h,
                     enablenox=yes,
                     enablenox=no))
                     
       if test "$enablenox" != no ; then
          AC_DEFINE(HAVE_NOX, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos NOX nonlinear solver])
          AC_MSG_RESULT(<<< Configuring library with NOX support >>>)
       fi
       
       dnl ------------------------------------------------------
       dnl ML
       dnl ------------------------------------------------------
       AC_CHECK_FILE($withtrilinosdir/include/ml_include.h,
                     enableml=yes,
       AC_CHECK_FILE($withtrilinosdir/packages/ml/src/Include/ml_include.h,
                     enableml=yes,
                     enableml=no))
                     
       if test "$enableml" != no ; then
          AC_DEFINE(HAVE_ML, 1,
                    [Flag indicating whether the library shall be compiled to use the Trilinos ML package])
          AC_MSG_RESULT(<<< Configuring library with ML support >>>)
       fi
    fi
  else
    enabletrilinos10=no
  fi

  AC_SUBST(TRILINOS_MAKEFILE_EXPORT)
  AC_SUBST(enabletrilinos10)
])


dnl -------------------------------------------------------------
dnl Trilinos 9
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRILINOS_9], 
[
  if test "x$TRILINOS_DIR" = "x"; then
    TRILINOS_DIR=no
  fi  	

  AC_ARG_WITH(trilinos,
              AC_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  dnl AztecOO
  if test "$withtrilinosdir" != no ; then
    AC_CHECK_FILE($withtrilinosdir/include/Makefile.export.aztecoo,
                  AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.aztecoo,
    AC_CHECK_FILE($withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo,
                  AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo,
	                enableaztecoo=no))

    if test "$enableaztecoo" != no ; then
       enableaztecoo=yes
       AC_DEFINE(HAVE_AZTECOO, 1,
                 [Flag indicating whether the library shall be compiled to use the Trilinos solver collection])
       AC_DEFINE(HAVE_TRILINOS, 1,
                 [])
       AC_MSG_RESULT(<<< Configuring library with Trilinos 9 support >>>)
    fi
  else
    enableaztecoo=no
  fi

  dnl Nox
  AC_ARG_WITH(nox,
              AC_HELP_STRING([--with-nox=PATH],[Specify the path to Nox installation]),
              withnoxdir=$withval,
              withnoxdir=$TRILINOS_DIR)
  if test "$withnoxdir" != no ; then
    AC_CHECK_FILE($withnoxdir/include/Makefile.export.nox,
                  NOX_MAKEFILE_EXPORT=$withnoxdir/include/Makefile.export.nox,
    AC_CHECK_FILE($withnoxdir/packages/nox/Makefile.export.nox,
                  NOX_MAKEFILE_EXPORT=$withnoxdir/packages/nox/Makefile.export.nox,
	 	              enablenox=no))

    if test "$enablenox" != no ; then
       enablenox=yes
       AC_DEFINE(HAVE_NOX, 1,
                 [Flag indicating whether the library shall be compiled to use the Nox solver collection])
       AC_MSG_RESULT(<<< Configuring library with Nox support >>>)
    fi
  else
    enablenox=no
  fi
  
  dnl ML
  AC_ARG_WITH(ml,
              AC_HELP_STRING([--with-ml=PATH],[Specify the path to ML installation]),
              withmldir=$withval,
              withmldir=$TRILINOS_DIR)
  if test "$withmldir" != no ; then
    AC_CHECK_FILE($withmldir/include/Makefile.export.ml,
                  ML_MAKEFILE_EXPORT=$withmldir/include/Makefile.export.ml,
    AC_CHECK_FILE($withmldir/packages/nox/Makefile.export.ml,
                  ML_MAKEFILE_EXPORT=$withmldir/packages/nox/Makefile.export.ml,
	 	              enableml=no))

    if test "$enableml" != no ; then
       enableml=yes
       AC_DEFINE(HAVE_ML, 1,
                 [Flag indicating whether the library shall be compiled to use the ML package])
       AC_MSG_RESULT(<<< Configuring library with ML support >>>)
    fi
  else
    enableml=no
  fi

  AC_SUBST(AZTECOO_MAKEFILE_EXPORT)
  AC_SUBST(enableaztecoo)

  AC_SUBST(NOX_MAKEFILE_EXPORT)
  AC_SUBST(enablenox)

  AC_SUBST(ML_MAKEFILE_EXPORT)
  AC_SUBST(enableml)

])
