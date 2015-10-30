dnl -------------------------------------------------------------
dnl Nemesis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NEMESIS], 
[
dnl Nemesis is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enablenetcdf = yes -a $enableexodus = yes -a $enablenemesis = yes); then
     NEMESIS_INCLUDE="-I$PWD/contrib/nemesis/Lib"
     NEMESIS_LIBRARY="\$(EXTERNAL_LIBDIR)/libnemesis\$(libext)"
     AC_DEFINE(HAVE_NEMESIS_API, 1, [Flag indicating whether the library will be compiled with Nemesis support])
     AC_MSG_RESULT(<<< Configuring library with Nemesis API support >>>)
  else
     NEMESIS_INCLUDE=""
     NEMESIS_LIBRARY=""
     enablenemesis=no
  fi

  AC_SUBST(NEMESIS_INCLUDE)
  AC_SUBST(NEMESIS_LIBRARY)	
  AC_SUBST(enablenemesis)
])
