dnl -------------------------------------------------------------
dnl libHilbert
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_LIBHILBERT],
[
dnl libHilbert is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enablelibhilbert = yes); then
     LIBHILBERT_INCLUDE="-I$PWD/contrib/libHilbert/include"
     LIBHILBERT_LIBRARY="\$(EXTERNAL_LIBDIR)/libHilbert\$(libext)"
     AC_DEFINE(HAVE_LIBHILBERT, 1, [Flag indicating whether the library will be compiled with libHilbert support])
     AC_MSG_RESULT(<<< Configuring library with libHilbert support >>>)
  else
     LIBHILBERT_INCLUDE=""
     LIBHILBERT_LIBRARY=""
     enablelibhilbert=no
  fi

  AC_SUBST(LIBHILBERT_INCLUDE)
  AC_SUBST(LIBHILBERT_LIBRARY)	
  AC_SUBST(enablelibhilbert)
])
