dnl -------------------------------------------------------------
dnl Space Filling Curves
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_SFC],
[
  AC_ARG_ENABLE(sfc,
                AS_HELP_STRING([--disable-sfc],
                               [build without space-filling curves support]),
                [case "${enableval}" in
                  yes)  enablesfc=yes ;;
                  no)  enablesfc=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-sfc) ;;
                esac],
                [enablesfc=$enableoptional])

  dnl The SFC API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablesfc = yes); then
     SFC_INCLUDE="-I\$(top_srcdir)/contrib/sfcurves"
     SFC_LIB="\$(EXTERNAL_LIBDIR)/libsfcurves\$(libext)"
     AC_DEFINE(HAVE_SFCURVES, 1, [Flag indicating whether the library will be compiled with SFC support])
     AC_MSG_RESULT(<<< Configuring library with SFC support >>>)
  else
     SFC_INCLUDE=""
     SFC_LIB=""
     enablesfc=no
  fi

  AC_SUBST(SFC_INCLUDE)
  AC_SUBST(SFC_LIB)
])
