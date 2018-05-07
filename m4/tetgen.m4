dnl -------------------------------------------------------------
dnl TetGen tetrahedralization library
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TETGEN],
[
  AC_ARG_ENABLE(tetgen,
                AS_HELP_STRING([--disable-tetgen],
                               [build without TetGen tetrahedralization library support]),
                [AS_CASE("${enableval}",
                         [yes], [enabletetgen=yes],
                         [no],  [enabletetgen=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-tetgen)])],
                [enabletetgen=$enableoptional])

  dnl The TETGEN API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  AS_IF([test "x$enabletetgen" = "xyes"],
        [
          TETGEN_INCLUDE="-I\$(top_srcdir)/contrib/tetgen"
          TETGEN_LIBRARY="\$(EXTERNAL_LIBDIR)/libtetgen\$(libext)"
          AC_DEFINE(HAVE_TETGEN, 1, [Flag indicating whether the library will be compiled with Tetgen support])
          AC_MSG_RESULT(<<< Configuring library with Tetgen support >>>)
        ],
        [
          TETGEN_INCLUDE=""
          TETGEN_LIBRARY=""
          enabletetgen=no
        ])

  AC_SUBST(TETGEN_INCLUDE)
  AC_SUBST(TETGEN_LIBRARY)
])
