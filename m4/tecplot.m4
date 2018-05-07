# -------------------------------------------------------------
# Tecplot
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TECPLOT],
[
  AC_ARG_ENABLE(tecplot,
                AS_HELP_STRING([--enable-tecplot],
                               [build with Tecplot binary file I/O support (using distributed libraries)]),
                [AS_CASE("${enableval}",
                         [yes], [enabletecplot=yes],
                         [no],  [enabletecplot=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-tecplot)])],
                [enabletecplot=no])

  dnl Can't support both vendor-provided libraries and building from source, and we prefer the latter
  AS_IF([test "x$enabletecplot" = "xyes" && test "x$enabletecio" = "xyes"],
        [
          AC_MSG_RESULT([>>> Not using vendor provided tecio libraries, deferring to source build <<<])
          enabletecplot=no
        ])

  dnl The Tecplot API is distributed with libmesh.
  AS_IF([test "x$enabletecplot" = "xyes"],
        [
          dnl We will check to see if we can actually link against the Tecplot library momentarily,
          dnl now we just see if the file exists, without using AC_CHECK_FILE!
          TECPLOT_LIBRARY_PATH=""
          AS_IF([test -r $top_srcdir/contrib/tecplot/binary/lib/$host/tecio.a],
                [TECPLOT_LIBRARY_PATH=$top_srcdir/contrib/tecplot/binary/lib/$host],
                [
                  AC_MSG_RESULT([>>> Configuring Tecplot failed, no tecio exists for $host <<<])
                  enabletecplot=no
                ])

          dnl Note: AC_CHECK_HEADER seems to fail if the path to the header
          dnl is a relative one, i.e containing ".."  in it.  We'll work around this
          dnl by setting the relevant path in $CPPFLAGS.
          old_CPPFLAGS="$CPPFLAGS"
          CPPFLAGS="$CPPFLAGS -I$top_srcdir/contrib/tecplot/binary/include"

          AC_CHECK_HEADER(TECIO.h,
                          [TECPLOT_INCLUDE_PATH=$top_srcdir/contrib/tecplot/binary/include
                           TECPLOT_INCLUDE="-I\$(top_srcdir)/contrib/tecplot/binary/include"])

          dnl Reset CPPFLAGS
          CPPFLAGS="$old_CPPFLAGS"

          dnl And don't step on anybody's toes
          unset old_CPPFLAGS
        ])

  AS_IF([test "x$enabletecplot" = "xyes"],
        [
          AS_IF([test -r $TECPLOT_LIBRARY_PATH/tecio.a && test -r $TECPLOT_INCLUDE_PATH/TECIO.h],
                [
                  dnl --------------------------------------------------------------------------
                  dnl  OK, the library and header are there, how about linking with the library?
                  dnl --------------------------------------------------------------------------
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
                ],
                [
                  AC_MSG_RESULT([>>> Configuring Tecplot failed, could not find at least one of tecio.a or TECIO.h <<<])
                  enabletecplot=no
                ])
        ])
])
