# -------------------------------------------------------------
# Tecplot
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TECPLOT],
[
  AC_ARG_ENABLE(tecplot,
                AC_HELP_STRING([--enable-tecplot],
                               [build with Tecplot binary file I/O support]),
		[case "${enableval}" in
		  yes)  enabletecplot=yes ;;
		   no)  enabletecplot=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-tecplot) ;;
		 esac],
		 [enabletecplot=$enableoptional])



  # The Tecplot API is distributed with libmesh, but we want to support external installations
  # on platforms we may not have the binaries for...
  if (test $enabletecplot = yes); then
    AC_ARG_WITH(tecplot,
                AC_HELP_STRING([--with-tecplot=PATH],[Specify the path where Tecplot is installed]),
                withtecplot=$withval,
                withtecplot=no)

    # unspecified - look in contrib  
    if test "$withtecplot" = no ; then
      AC_CHECK_FILE($top_srcdir/contrib/tecplot/lib/$host/tecio.a,
  	  	  TECPLOT_LIBRARY_PATH=$top_srcdir/contrib/tecplot/lib/$host)
      AC_CHECK_FILE($top_srcdir/contrib/tecplot/include/TECIO.h,
   	  	  TECPLOT_INCLUDE_PATH=$top_srcdir/contrib/tecplot/include)

    # specified - look there
    else
      AC_CHECK_FILE($withtecplot/lib/tecio.a,
  	  	  TECPLOT_LIBRARY_PATH=$withtecplot/lib)
      AC_CHECK_FILE($withtecplot/include/TECIO.h,
   	  	  TECPLOT_INCLUDE_PATH=$withtecplot/include)
    fi
  
    if (test -r $TECPLOT_LIBRARY_PATH/tecio.a -a -r $TECPLOT_INCLUDE_PATH/TECIO.h) ; then
  
      #--------------------------------------------------------------------------
      # OK, the library and header are there, how about linking with the library?
      #--------------------------------------------------------------------------
      save_CPPFLAGS=$CPPFLAGS
      save_LIBS=$LIBS
  
      CPPFLAGS="-I$TECPLOT_INCLUDE_PATH $CPPFLAGS"
      LIBS="$TECPLOT_LIBRARY_PATH/tecio.a $LIBS"
  
      AC_LINK_IFELSE(
                  [
                     AC_LANG_PROGRAM([#include <TECIO.h>], 
                                     [int ierr = TECEND112 ();])
                  ],
                  [
                     TECPLOT_LIBRARY=$TECPLOT_LIBRARY_PATH/tecio.a
                     TECPLOT_INCLUDE=-I$TECPLOT_INCLUDE_PATH
                     AC_SUBST(TECPLOT_LIBRARY)
                     AC_SUBST(TECPLOT_INCLUDE)
                     AC_DEFINE(HAVE_TECPLOT_API, 1,
                               [Flag indicating whether the library shall be compiled to use the Tecplot interface])
                     AC_DEFINE(HAVE_TECPLOT_API_112, 1,
                               [Flag indicating tecplot API understands newer features])
                     AC_MSG_RESULT(<<< Configuring library with Tecplot API support (v11.2) >>>)
		     libmesh_contrib_INCLUDES="$TECPLOT_INCLUDE $libmesh_contrib_INCLUDES"
                  ],
		  [		  
                     AC_LINK_IFELSE(
                                 [
                                    AC_LANG_PROGRAM([#include <TECIO.h>], 
                                                    [int ierr = TECEND ();])
                                 ],
                                 [
                                    TECPLOT_LIBRARY=$TECPLOT_LIBRARY_PATH/tecio.a
                                    TECPLOT_INCLUDE=-I$TECPLOT_INCLUDE_PATH
                                    AC_SUBST(TECPLOT_LIBRARY)
                                    AC_SUBST(TECPLOT_INCLUDE)
                                    AC_DEFINE(HAVE_TECPLOT_API, 1,
                                              [Flag indicating whether the library shall be compiled to use the Tecplot interface])
                                    AC_MSG_RESULT(<<< Configuring library with legacy Tecplot API support >>>)
               		            libmesh_contrib_INCLUDES="$TECPLOT_INCLUDE $libmesh_contrib_INCLUDES"
                                 ],
                                 [
                                    AC_MSG_RESULT( [WARNING: Found $TECPLOT_LIBRARY_PATH/tecio.a but cannot link with it!] )
               		            enabletecplot=no
                                 ] )
	          ] )
				 
      LIBS=$save_LIBS
      CPPFLAGS=$save_CPPFLAGS
    else
      enabletecplot=no
    fi
  fi



  AC_SUBST(enabletecplot)

  AM_CONDITIONAL(LIBMESH_ENABLE_TECPLOT, test x$enabletecplot = xyes)		 
])
