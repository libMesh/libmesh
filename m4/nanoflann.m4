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

  dnl The Nanoflann PointLocator is not used by default, you can turn it on by configuring
  dnl with --enable-nanoflann-pointlocator
  AC_ARG_ENABLE(nanoflann-pointlocator,
                AS_HELP_STRING([--enable-nanoflann-pointlocator],
                               [use Nanoflann-based PointLocator (experimental)]),
                [AS_CASE("${enableval}",
                         [yes], [enablenanoflannpointlocator=yes],
                         [no], [enablenanoflannpointlocator=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-nanoflann-pointlocator)])],
                [enablenanoflannpointlocator=no])

  dnl The NANOFLANN API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  AS_IF([test "x$enablenanoflann" = "xyes"],
        [
          NANOFLANN_INCLUDE="-I\$(top_srcdir)/contrib/nanoflann/include"
          AC_DEFINE(HAVE_NANOFLANN, 1, [Flag indicating whether the library will be compiled with nanoflann KD-Tree support])
          AC_MSG_RESULT(<<< Configuring library with nanoflann KDtree support >>>)

          AS_IF([test "x$enablenanoflannpointlocator" = "xyes"],
                [
                  AC_DEFINE(ENABLE_NANOFLANN_POINTLOCATOR, 1, [Flag indicating whether the library will use the (experimental) Nanoflann-based PointLocator])
                  AC_MSG_RESULT([<<< Configuring library with Nanoflann-based PointLocator >>>])
                ])
        ],
        [
          NANOFLANN_INCLUDE=""
          enablenanoflann=no
        ])

  AC_SUBST(NANOFLANN_INCLUDE)
])
