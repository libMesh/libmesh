dnl -------------------------------------------------------------
dnl Tecplot
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TECPLOT],
[
  AC_ARG_WITH(tecplot,
              AC_HELP_STRING([--with-tecplot=PATH],[Specify the path where Tecplot is installed]),
              withtecplot=$withval,
              withtecplot=no)

  if test "$withtecplot" = no ; then
    AC_CHECK_FILE(./contrib/tecplot/lib/$host/tecio.a,
	  	  TECPLOT_LIBRARY_PATH=$PWD/contrib/tecplot/lib/$host)
    AC_CHECK_FILE(./contrib/tecplot/include/TECIO.h,
 	  	  TECPLOT_INCLUDE_PATH=$PWD/contrib/tecplot/include)
  else
    AC_CHECK_FILE($withtecplot/lib/tecio.a,
	  	  TECPLOT_LIBRARY_PATH=$withtecplot/lib)
    AC_CHECK_FILE($withtecplot/include/TECIO.h,
 	  	  TECPLOT_INCLUDE_PATH=$withtecplot/include)
  fi

  if (test -r $TECPLOT_LIBRARY_PATH/tecio.a -a -r $TECPLOT_INCLUDE_PATH/TECIO.h) ; then

    dnl--------------------------------------------------------------------------
    dnl OK, the library and header are there, how about linking with the library?
    dnl--------------------------------------------------------------------------
    save_CPPFLAGS=$CPPFLAGS
    save_LIBS=$LIBS

    CPPFLAGS="-I$TECPLOT_INCLUDE_PATH $CPPFLAGS"
    LIBS="$TECPLOT_LIBRARY_PATH/tecio.a $LIBS"

    AC_TRY_LINK([#include <TECIO.h>],
                [ int ierr = TECEND ();],
                [
                   TECPLOT_LIBRARY=$TECPLOT_LIBRARY_PATH/tecio.a
                   TECPLOT_INCLUDE=-I$TECPLOT_INCLUDE_PATH
                   AC_SUBST(TECPLOT_LIBRARY)
                   AC_SUBST(TECPLOT_INCLUDE)
                   AC_DEFINE(HAVE_TECPLOT_API, 1,
                             [Flag indicating whether the library shall be compiled to use the Tecplot interface])
                   AC_MSG_RESULT(<<< Configuring library with Tecplot API support >>>)
                ],
                [
                   AC_MSG_RESULT( [WARNING: Found $TECPLOT_LIBRARY_PATH/tecio.a but cannot link with it!] )
                ] )
    LIBS=$save_LIBS
    CPPFLAGS=$save_CPPFLAGS
  fi
])
