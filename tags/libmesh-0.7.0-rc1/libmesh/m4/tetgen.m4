dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl TetGen tetrahedrization library
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TETGEN], 
[
dnl if TetGen is enabled we need the header path and the lib

  if (test $enabletetgen = yes) ; then
     TETGEN_INCLUDE="-I$PWD/contrib/tetgen"
     TETGEN_LIBRARY="\$(EXTERNAL_LIBDIR)/libtetgen\$(libext)"
     AC_DEFINE(HAVE_TETGEN, 1, [Flag indicating whether the library will be compiled with TetGen support])
     AC_MSG_RESULT(<<< Configuring library with TetGen support >>>)
  else
     TETGEN_INCLUDE=""
     TETGEN_LIBRARY=""
     enabletetgen=no
   fi

  dnl TetGen
  AC_SUBST(TETGEN_INCLUDE)
  AC_SUBST(TETGEN_LIBRARY)	
  AC_SUBST(enabletetgen)
])
