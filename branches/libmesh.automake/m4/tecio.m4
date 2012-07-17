dnl -------------------------------------------------------------
dnl tecio
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TECIO], 
[
  AC_ARG_ENABLE(tecio,
                AC_HELP_STRING([--enable-tecio],
                               [build with Tecplot TecIO API support (from source)]),
		[case "${enableval}" in
		  yes)  enabletecio=yes ;;
		   no)  enabletecio=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-tecio) ;;
		 esac],
		 [enabletecio=no]) #disabled by default, do more portability.


		
  dnl The TECIO API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enabletecio = yes); then
  
    # tecio platform-specific compiler flags
    TECIO_CPPFLAGS=""
    case "${host_os}" in
      *linux*)
	TECIO_CPPFLAGS="-DLINUX $TECIO_CPPFLAGS"
	AC_CHECK_SIZEOF([void *])
	if (test $ac_cv_sizeof_void_p = 8); then
	  TECIO_CPPFLAGS="-DLINUX64 $TECIO_CPPFLAGS"
	fi			 
	;;
    
      *darwin*)
	TECIO_CPPFLAGS="-DDARWIN -DLONGIS64 $TECIO_CPPFLAGS"
        ;;
	
        *)
	AC_MSG_RESULT([>>> Unrecognized TecIO platform, see contrib/tecplot/tecio/Runmake for hints on how to extend <<<])
	;;
    esac	

  
     TECIO_INCLUDE="-I\$(top_srcdir)/contrib/tecplot/tecio/tecsrc"
     TECIO_LIBRARY="\$(EXTERNAL_LIBDIR)/libtecio\$(libext)"
     AC_DEFINE(HAVE_TECPLOT_API, 1, [Flag indicating whether the library will be compiled with Tecplot TecIO API support])
     AC_MSG_RESULT(<<< Configuring library with Tecplot TecIO support >>>)
     libmesh_contrib_INCLUDES="$TECIO_INCLUDE $libmesh_contrib_INCLUDES"
     have_tecio=yes
  else
     TECIO_INCLUDE=""
     TECIO_LIBRARY=""
     enabletecio=no
     have_tecio=no
  fi

  AC_SUBST(TECIO_INCLUDE)
  AC_SUBST(TECIO_CPPFLAGS)
  AC_SUBST(enabletecio)
		 
  AM_CONDITIONAL(LIBMESH_ENABLE_TECIO, test x$enabletecio = xyes)
])
