# -------------------------------------------------------------
# nanoflann - a C++ header-only library for building KD-Trees
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NANOFLANN],
[
  AC_ARG_ENABLE(nanoflann,
                AS_HELP_STRING([--disable-nanoflann],
                               [build without nanoflann KD-tree support]),
                [AS_CASE("${enableval}",
                         [yes], [enablenanoflann=yes],
                         [no],  [enablenanoflann=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-nanoflann)])],
                [enablenanoflann=$enableoptional])

  dnl The NANOFLANN API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  AS_IF([test "x$enablenanoflann" = "xyes"],
        [
          NANOFLANN_INCLUDE="-I\$(top_srcdir)/contrib/nanoflann/include"
          AC_DEFINE(HAVE_NANOFLANN, 1, [Flag indicating whether the library will be compiled with nanoflann KD-Tree support])
          AC_MSG_RESULT(<<< Configuring library with nanoflann KDtree support >>>)
        ],
        [
          NANOFLANN_INCLUDE=""
          enablenanoflann=no
        ])

  AC_SUBST(NANOFLANN_INCLUDE)
])
