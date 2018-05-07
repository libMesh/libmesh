# -------------------------------------------------------------
# qhull
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_QHULL],
[
  AC_ARG_ENABLE(qhull,
                AS_HELP_STRING([--enable-qhull],
                               [build with Qhull API support]),
                [AS_CASE("${enableval}",
                         [yes], [enableqhull=yes],
                         [no],  [enableqhull=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-qhull)])],
                [enableqhull=$enableoptional]) # if unspecified, depend on enableoptional

  QHULL_INCLUDE=""
  QHULL_LIBS=""

  dnl fix for --disable-optional
  AS_IF([test "x$enableqhull" = "xyes"],
        [
          dnl look for require -lm libraries
          AC_SEARCH_LIBS([sqrt],  [m], [QHULL_LIBS="-lm"])
          AC_SEARCH_LIBS([trunc], [m], [QHULL_LIBS="-lm"])
          AC_SEARCH_LIBS([cos],   [m], [QHULL_LIBS="-lm"])

          dnl The QHULL API is distributed with libmesh, so we don't have to guess
          dnl where it might be installed...
          QHULL_INCLUDE="-I\$(top_srcdir)/contrib/qhull/qhull/src -I\$(top_srcdir)/contrib/qhull/qhull/src/libqhullcpp"
          AC_DEFINE(HAVE_QHULL_API, 1, [Flag indicating whether the library will be compiled with Qhull support])
          AC_MSG_RESULT(<<< Configuring library with Qhull version 2012.1 support >>>)
        ])

  AC_SUBST(QHULL_LIBS)
])
