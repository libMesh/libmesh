dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl Nox
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NOX], 
[
  if test "x$TRILINOS_DIR" = "x"; then
    TRILINOS_DIR=no
  fi  	

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
       AC_DEFINE(HAVE_NOX, 1,
                 [Flag indicating whether the library shall be compiled to use the Nox solver collection])
       AC_MSG_RESULT(<<< Configuring library with Nox support >>>)
    fi
  else
    enablenox=no
  fi

  AC_SUBST(NOX_MAKEFILE_EXPORT)
  AC_SUBST(enablenox)

])
