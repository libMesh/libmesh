dnl -------------------------------------------------------------
dnl TetGen tetrahedrization library
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TETGEN], 
[
  AC_ARG_ENABLE(tetgen,
                AC_HELP_STRING([--enable-tetgen],
                               [build with TetGen tetrahedrization library support]),
		[case "${enableval}" in
		  yes)  enabletetgen=yes ;;
		   no)  enabletetgen=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-tetgen) ;;
		 esac],
		 [enabletetgen=$enableoptional])


		 		 
  dnl The TETGEN API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enabletetgen = yes); then
     TETGEN_INCLUDE="-I\$(top_srcdir)/contrib/tetgen"
     #TETGEN_LIBRARY="\$(EXTERNAL_LIBDIR)/libtetgen\$(libext)"
     AC_DEFINE(HAVE_TETGEN, 1, [Flag indicating whether the library will be compiled with Tetgen support])
     AC_MSG_RESULT(<<< Configuring library with Tetgen support >>>)
     libmesh_contrib_INCLUDES="$TETGEN_INCLUDE $libmesh_contrib_INCLUDES"
  else
     TETGEN_INCLUDE=""
     #TETGEN_LIBRARY=""
     enabletetgen=no
  fi

  AC_SUBST(TETGEN_INCLUDE)
  #AC_SUBST(TETGEN_LIBRARY)	
  AC_SUBST(enabletetgen)

  AM_CONDITIONAL(LIBMESH_ENABLE_TETGEN, test x$enabletetgen = xyes)
])
