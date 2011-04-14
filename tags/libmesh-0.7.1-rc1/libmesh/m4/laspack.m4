dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl LASPACK Iterative Solvers
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_LASPACK], 
[

dnl Initialize variables
LASPACK_INCLUDE=""
LASPACK_LIB=""

dnl Sanity check: make sure the user really has the contrib directory
if (test $enablelaspack = yes); then
  AC_CHECK_FILE(./contrib/laspack/lastypes.h, [enablelaspack=yes], [enablelaspack=no])
fi


if (test $enablelaspack = yes); then

  LASPACK_INCLUDE="-I$PWD/contrib/laspack"
  LASPACK_LIB="\$(EXTERNAL_LIBDIR)/liblaspack\$(libext)"
  AC_DEFINE(HAVE_LASPACK, 1, [Flag indicating whether or not LASPACK iterative solvers are available])
  laspack_version=`grep "define LASPACK_VERSION " $PWD/contrib/laspack/version.h | sed -e "s/[[^0-9.]]*//g"`
  AC_MSG_RESULT(<<< Configuring library with LASPACK version $laspack_version support >>>)

fi

AC_SUBST(LASPACK_INCLUDE)
AC_SUBST(LASPACK_LIB)	
AC_SUBST(enablelaspack)

])
