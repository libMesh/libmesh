dnl -------------------------------------------------------------
dnl fparser
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_FPARSER],
[
dnl fparser is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enablefparser = yes); then
     FPARSER_INCLUDE="-I$PWD/contrib/fparser"
     FPARSER_LIBRARY="\$(EXTERNAL_LIBDIR)/libfparser\$(libext)"
     AC_DEFINE(HAVE_FPARSER, 1, [Flag indicating whether the library will be compiled with fparser support])
     AC_MSG_RESULT(<<< Configuring library with fparser support >>>)
  else
     FPARSER_INCLUDE=""
     FPARSER_LIBRARY=""
     enablefparser=no
  fi

  AC_SUBST(FPARSER_INCLUDE)
  AC_SUBST(FPARSER_LIBRARY)	
  AC_SUBST(enablefparser)
])
