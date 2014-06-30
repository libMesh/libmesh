dnl -------------------------------------------------------------
dnl Triangle Delaunay triangulation library by J.R. Shewchuk
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TRIANGLE],
[
  AC_ARG_ENABLE(triangle,
                AS_HELP_STRING([--disable-triangle],
                               [build without Triangle Delaunay triangulation library support]),
		[case "${enableval}" in
		  yes)  enabletriangle=yes ;;
		   no)  enabletriangle=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-triangle) ;;
		 esac],
		 [enabletriangle=$enableoptional])



  dnl The TRIANGLE API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enabletriangle = yes); then
     TRIANGLE_INCLUDE="-I\$(top_srcdir)/contrib/triangle"
     TRIANGLE_LIBRARY="\$(EXTERNAL_LIBDIR)/libtriangle\$(libext)"
     AC_DEFINE(HAVE_TRIANGLE, 1, [Flag indicating whether the library will be compiled with Triangle support])
     AC_MSG_RESULT(<<< Configuring library with Triangle support >>>)
  else
     TRIANGLE_INCLUDE=""
     TRIANGLE_LIBRARY=""
     enabletriangle=no
  fi

  AC_SUBST(TRIANGLE_INCLUDE)
  AC_SUBST(TRIANGLE_LIBRARY)
])
