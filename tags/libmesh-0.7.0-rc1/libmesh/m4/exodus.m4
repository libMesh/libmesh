dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl ExodusII
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_EXODUS], 
[
dnl Exodus is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enablenetcdf = yes -a $enableexodus = yes); then
     EXODUS_INCLUDE="-I$PWD/contrib/exodusii/Lib/include"
     EXODUS_LIBRARY="\$(EXTERNAL_LIBDIR)/libexodusii\$(libext)"
     AC_DEFINE(HAVE_EXODUS_API, 1, [Flag indicating whether the library will be compiled with Exodus support])
     AC_MSG_RESULT(<<< Configuring library with Exodus API support >>>)
  else
     EXODUS_INCLUDE=""
     EXODUS_LIBRARY=""
     enableexodus=no
  fi

  AC_SUBST(EXODUS_INCLUDE)
  AC_SUBST(EXODUS_LIBRARY)	
  AC_SUBST(enableexodus)
])
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl ExodusII 
dnl -------------------------------------------------------------
dnl AC_DEFUN([CONFIGURE_EXODUS],
dnl [
dnl   AC_CHECK_FILE(./contrib/exodus/lib/$host/libexoIIv2c.a,
dnl 		EXODUS_LIB=$PWD/contrib/exodus/lib/$host/libexoIIv2c.a)
dnl   AC_CHECK_FILE(./contrib/exodus/include/exodusII.h,
dnl 		EXODUS_INCLUDE_PATH=$PWD/contrib/exodus/include)

dnl   if (test -r $EXODUS_INCLUDE_PATH/exodusII.h -a "x$EXODUS_LIB" != x) ; then
dnl     EXODUS_INCLUDE=-I$EXODUS_INCLUDE_PATH
dnl     AC_SUBST(EXODUS_LIB)
dnl     AC_SUBST(EXODUS_INCLUDE)
dnl     AC_DEFINE(HAVE_EXODUS_API, 1,
dnl 	      [Flag indicating whether the library shall be compiled to use the Exodus interface])
dnl     AC_MSG_RESULT(<<< Configuring library with Exodus API support >>>)
dnl   fi
dnl ])
