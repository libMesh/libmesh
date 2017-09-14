dnl -------------------------------------------------------------
dnl LASPACK Iterative Solvers
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_LASPACK],
[
  AC_ARG_ENABLE(laspack,
                AS_HELP_STRING([--disable-laspack],
                               [build without LASPACK iterative solver support]),
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
