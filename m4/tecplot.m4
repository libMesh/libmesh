# -------------------------------------------------------------
# Tecplot
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TECPLOT],
[
  AC_ARG_ENABLE(tecplot,
                AS_HELP_STRING([--enable-tecplot],
                               [build with Tecplot binary file I/O support (using distributed libraries)]),
                [case "${enableval}" in
                  yes)  enabletecplot=yes ;;
                  no)  enabletecplot=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-tecplot) ;;
                esac],
                [enabletecplot=no])

  # Can't support both vendor-provided libraries and building from source, and we prefer the latter
  if (test "x$enabletecplot" = "xyes" -a "x$enabletecio" = "xyes"); then
     AC_MSG_RESULT([>>> Not using vendor provided tecio libraries, deferring to source build <<<])
     enabletecplot=no
  fi


  # The Tecplot API is distributed with libmesh.
  if (test $enabletecplot = yes); then

      # We will check to see if we can actually link against the Tecplot library momentarily,
      # now we just see if the file exists, without using AC_CHECK_FILE!
      TECPLOT_LIBRARY_PATH=""
      if (test -r $top_srcdir/contrib/tecplot/binary/lib/$host/tecio.a) ; then
        TECPLOT_LIBRARY_PATH=$top_srcdir/contrib/tecplot/binary/lib/$host
      else
        AC_MSG_RESULT([>>> Configuring Tecplot failed, no tecio exists for $host <<<])
        enabletecplot=no
      fi

      # Note: AC_CHECK_HEADER seems to fail if the path to the header
      # is a relative one, i.e containing ".."  in it.  We'll work around this
      # by setting the relevant path in $CPPFLAGS.
      old_CPPFLAGS="$CPPFLAGS"
      CPPFLAGS="$CPPFLAGS -I$top_srcdir/contrib/tecplot/binary/include"

      AC_CHECK_HEADER(TECIO.h,
                      [TECPLOT_INCLUDE_PATH=$top_srcdir/contrib/tecplot/binary/include
                       TECPLOT_INCLUDE="-I\$(top_srcdir)/contrib/tecplot/binary/include"])

      # Reset CPPFLAGS
      CPPFLAGS="$old_CPPFLAGS"

      # And don't step on anybody's toes
      unset old_CPPFLAGS
  fi

  if (test $enabletecplot = yes); then
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
                     AC_LANG_PROGRAM([@%:@include <TECIO.h>],
                                     [int ierr = TECEND112 ();])
                  ],
                  [
                     AC_DEFINE(HAVE_TECPLOT_API, 1,
                               [Flag indicating whether the library will be compiled with Tecplot TecIO API support])
                     AC_DEFINE(HAVE_TECPLOT_API_112, 1,
                               [Flag indicating tecplot API understands newer features])
                     AC_MSG_RESULT(<<< Configuring library with Tecplot API support (v11.2) >>>)
                  ],
                  [
                     AC_LINK_IFELSE(
                                 [
                                    AC_LANG_PROGRAM([@%:@include <TECIO.h>],
                                                    [int ierr = TECEND ();])
                                 ],
                                 [
                                    AC_DEFINE(HAVE_TECPLOT_API, 1,
                                              [Flag indicating whether the library shall be compiled to use the Tecplot interface])
                                    AC_MSG_RESULT(<<< Configuring library with legacy Tecplot API support >>>)
                                 ],
                                 [
                                    AC_MSG_RESULT( [WARNING: Found $TECPLOT_LIBRARY_PATH/tecio.a but cannot link with it!] )
                                    enabletecplot=no
                                 ])
                  ])

      LIBS=$save_LIBS
      CPPFLAGS=$save_CPPFLAGS
    else
      AC_MSG_RESULT([>>> Configuring Tecplot failed, could not find at least one of tecio.a or TECIO.h <<<])
      enabletecplot=no
  fi
  fi
])
