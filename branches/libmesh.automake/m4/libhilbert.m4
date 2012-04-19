dnl -------------------------------------------------------------
dnl libHilbert
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_LIBHILBERT], 
[
  AC_ARG_ENABLE(libhilbert,
                AC_HELP_STRING([--enable-libHilbert],
                               [build with libHilbert, from Chris Hamilton]),
		[case "${enableval}" in
		  yes)  enablelibhilbert=yes ;;
		   no)  enablelibhilbert=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-libhilbert) ;;
		 esac],
		 [enablelibhilbert=$enableoptional])

		 
		
  dnl The LIBHILBERT API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablelibhilbert = yes); then
     LIBHILBERT_INCLUDE="-I\$(top_srcdir)/contrib/libHilbert/include"
     #LIBHILBERT_LIBRARY="\$(EXTERNAL_LIBDIR)/liblibhilbert\$(libext)"
     AC_DEFINE(HAVE_LIBHILBERT, 1, [Flag indicating whether the library will be compiled with libHilbert support])
     AC_MSG_RESULT(<<< Configuring library with libHilbert support >>>)
     libmesh_contrib_INCLUDES="$LIBHILBERT_INCLUDE $libmesh_contrib_INCLUDES"
  else
     LIBHILBERT_INCLUDE=""
     #LIBHILBERT_LIBRARY=""
     enablelibhilbert=no
  fi

  AC_SUBST(LIBHILBERT_INCLUDE)
  #AC_SUBST(LIBHILBERT_LIBRARY)	
  AC_SUBST(enablelibhilbert)

  AM_CONDITIONAL(LIBMESH_ENABLE_LIBHILBERT, test x$enablelibhilbert = xyes)
])
