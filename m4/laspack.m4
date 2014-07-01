dnl -------------------------------------------------------------
dnl LASPACK Iterative Solvers
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_LASPACK],
[
  AC_ARG_ENABLE(laspack,
                AS_HELP_STRING([--disable-laspack],
                               [build without LASPACK iterative solver suppport]),
		[case "${enableval}" in
		  yes)  enablelaspack=yes ;;
		   no)  enablelaspack=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-laspack) ;;
		 esac],
		 [enablelaspack=$enableoptional])



  dnl The LASPACK API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablelaspack = yes); then
     LASPACK_INCLUDE="-I\$(top_srcdir)/contrib/laspack"
     LASPACK_LIB="\$(EXTERNAL_LIBDIR)/liblaspack\$(libext)"
     AC_DEFINE(HAVE_LASPACK, 1, [Flag indicating whether the library will be compiled with LASPACK support])
     AC_MSG_RESULT(<<< Configuring library with Laspack support >>>)
  else
     LASPACK_INCLUDE=""
     LASPACK_LIB=""
     enablelaspack=no
  fi

  AC_SUBST(LASPACK_INCLUDE)
  AC_SUBST(LASPACK_LIB)
])


# AC_DEFUN([CONFIGURE_LASPACK],
# [

# dnl Initialize variables
# LASPACK_INCLUDE=""
# LASPACK_LIB=""

# dnl Sanity check: make sure the user really has the contrib directory
# if (test $enablelaspack = yes); then
#   AC_CHECK_HEADER(./contrib/laspack/lastypes.h, [enablelaspack=yes], [enablelaspack=no])
# fi


# if (test $enablelaspack = yes); then

#   LASPACK_INCLUDE="-I\$(top_srcdir)/contrib/laspack"
#   LASPACK_LIB="\$(EXTERNAL_LIBDIR)/liblaspack\$(libext)"
#   AC_DEFINE(HAVE_LASPACK, 1, [Flag indicating whether or not LASPACK iterative solvers are available])
#   laspack_version=`grep "define LASPACK_VERSION " $PWD/contrib/laspack/version.h | sed -e "s/[[^0-9.]]*//g"`
#   AC_MSG_RESULT(<<< Configuring library with LASPACK version $laspack_version support >>>)

# fi

# AC_SUBST(LASPACK_INCLUDE)
# AC_SUBST(LASPACK_LIB)
# AC_SUBST(enablelaspack)

# ])
