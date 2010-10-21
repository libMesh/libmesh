dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl Trilinos
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRILINOS], 
[
  if test "x$TRILINOS_DIR" = "x"; then
    TRILINOS_DIR=no
  fi  	

  AC_ARG_WITH(trilinos,
              AC_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  if test "$withtrilinosdir" != no ; then
    AC_CHECK_FILE($withtrilinosdir/include/Makefile.export.aztecoo,
                  AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.aztecoo,
                  AC_CHECK_FILE($withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo,
                                AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/packages/aztecoo/Makefile.export.aztecoo,
	 	                enabletrilinos=no))

    if test "$enabletrilinos" != no ; then
       AC_DEFINE(HAVE_TRILINOS, 1,
                 [Flag indicating whether the library shall be compiled to use the Trilinos solver collection])
       AC_MSG_RESULT(<<< Configuring library with Trilinos support >>>)
    fi
  else
    enabletrilinos=no
  fi

  AC_SUBST(AZTECOO_MAKEFILE_EXPORT)
  AC_SUBST(enabletrilinos)

])
