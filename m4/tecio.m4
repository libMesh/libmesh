dnl -------------------------------------------------------------
dnl tecio
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TECIO],
[
  AC_ARG_ENABLE(tecio,
                AS_HELP_STRING([--disable-tecio],
                               [build without Tecplot TecIO API support (from source)]),
                [AS_CASE("${enableval}",
                         [yes], [enabletecio=yes],
                         [no],  [enabletecio=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-tecio)])],
                [enabletecio=$enableoptional])

  AC_ARG_WITH(tecio-x11-include,
              AS_HELP_STRING([--with-tecio-x11-include=PATH],
                             [Path to X11 header files. E.g. /usr/include but _not_ /usr/include/X11]),
              withteciox11inc=$withval,
              withteciox11inc=no)

  dnl Start off with an empty TECIO_CPPFLAGS and possibly append to it.
  TECIO_CPPFLAGS=""

  dnl Look for an X11 installation
  AS_IF([test "x$enabletecio" = "xyes"],
        [
          dnl If the user specified a path to look for X11 headers in, honor
          dnl that and don't try to look elsewhere.  If the compilation fails
          dnl then they must figure out why on their own.
          AS_IF([test "x$withteciox11inc" != "xno"], [TECIO_CPPFLAGS="-I$withteciox11inc"],
                dnl The user did not specify where to look, so see if the file
                dnl exists in the usual Linux location...
                [test -r /usr/include/X11/Intrinsic.h], [TECIO_CPPFLAGS="-I/usr/include"],
                dnl ... and if not there, try the Mac (XQuartz) location.
                [test -r /opt/X11/include/X11/Intrinsic.h], [TECIO_CPPFLAGS="-I/opt/X11/include"])

          dnl Print a status message
          AS_IF([test "x$TECIO_CPPFLAGS" != "x"],
                [AC_MSG_RESULT(<<< Testing X11 headers with $TECIO_CPPFLAGS >>>)])

          dnl Run AC_CHECK_HEADER to see if we can actually compile a test
          dnl program with the TECIO_CPPFLAGS we determined.  If this doesn't
          dnl work, we assume the X11 installation is not working and disable
          dnl the tecio interface.
          saved_CPPFLAGS=$CPPFLAGS
          CPPFLAGS="$TECIO_CPPFLAGS $CPPFLAGS"

          AC_CHECK_HEADER([X11/Intrinsic.h],
                          [],
                          [enabletecio=no])

          dnl Restore original value of CPPFLAGS
          CPPFLAGS=$saved_CPPFLAGS
        ])

  dnl The TECIO API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  AS_IF([test "x$enabletecio" = "xyes"],
        [
          dnl Set tecio platform-specific compiler flags.
          AS_CASE("${host_os}",
                  [*linux*], [TECIO_CPPFLAGS="-DLINUX $TECIO_CPPFLAGS"
                              AC_CHECK_SIZEOF([void *])
                              AS_IF([test "x$ac_cv_sizeof_void_p" = "x8"], [TECIO_CPPFLAGS="-DLINUX64 $TECIO_CPPFLAGS"])],
                  [*darwin*], [TECIO_CPPFLAGS="-DDARWIN -DLONGIS64 $TECIO_CPPFLAGS"],
                  [AC_MSG_RESULT([>>> Unrecognized TecIO platform, see contrib/tecplot/tecio/Runmake for hints on how to extend <<<])])

           TECIO_INCLUDE="-I\$(top_srcdir)/contrib/tecplot/tecio/include"
           AC_DEFINE(HAVE_TECPLOT_API, 1, [Flag indicating whether the library will be compiled with Tecplot TecIO API support])
           AC_DEFINE(HAVE_TECPLOT_API_112, 1, [Flag indicating tecplot API understands newer features])
           AC_MSG_RESULT(<<< Configuring library with Tecplot TecIO support >>>)
           have_tecio=yes
        ],
        [
          TECIO_INCLUDE=""
          enabletecio=no
          have_tecio=no
        ])

  AC_SUBST(TECIO_INCLUDE)
  AC_SUBST(TECIO_CPPFLAGS)
])
