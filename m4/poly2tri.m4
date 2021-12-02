# -------------------------------------------------------------
# poly2tri
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_POLY2TRI],
[
  AC_ARG_ENABLE(poly2tri,
                AS_HELP_STRING([--enable-poly2tri],
                               [build with poly2tri API support]),
                [AS_CASE("${enableval}",
                         [yes], [enablepoly2tri=yes],
                         [no],  [enablepoly2tri=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-poly2tri)])],
                [enablepoly2tri=$enableoptional]) # if unspecified, depend on enableoptional

  POLY2TRI_INCLUDE=""

  dnl fix for --disable-optional
  AS_IF([test "x$enablepoly2tri" = "xyes"],
        [
          dnl The poly2tri API is distributed with libmesh, so we don't have to guess
          dnl where it might be installed...
          POLY2TRI_INCLUDE="-I\$(top_builddir)/contrib/poly2tri/modified"
          AC_DEFINE(HAVE_POLY2TRI, 1, [Flag indicating whether the library will be compiled with poly2tri support])
          AC_MSG_RESULT(<<< Configuring library with poly2tri support >>>)
        ])
])
