dnl -------------------------------------------------------------
dnl libHilbert
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_LIBHILBERT],
[
  AC_ARG_ENABLE(libHilbert,
                AS_HELP_STRING([--disable-libHilbert],
                               [build without Chris Hamilton's libHilbert]),
		[case "${enableval}" in
		  yes)  enablelibhilbert=yes ;;
		   no)  enablelibhilbert=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-libHilbert) ;;
		 esac],
		 [enablelibhilbert=$enableoptional])



  dnl The LIBHILBERT API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablelibhilbert = yes); then
     LIBHILBERT_INCLUDE="-I\$(top_srcdir)/contrib/libHilbert/include"
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
])
