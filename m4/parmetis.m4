dnl -------------------------------------------------------------
dnl Parmetis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PARMETIS],
[
  AC_ARG_ENABLE(parmetis,
                AS_HELP_STRING([--disable-parmetis],
                               [build without Parmetis parallel graph partitioning support]),
                [AS_CASE("${enableval}",
                         [yes], [enableparmetis=yes],
                         [no],  [enableparmetis=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-parmetis)])],
                [enableparmetis=$enableoptional])

  dnl Trump --enable-parmetis with --disable-mpi
  AS_IF([test "x$enablempi" = "xno"], [enableparmetis=no])

  dnl The PARMETIS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  AS_IF([test "x$enableparmetis" = "xyes"],
        [
          dnl We always configure for Metis before ParMetis, and we will respect the value of $build_metis
          dnl determined there -- if we are using PETSc's Metis, we'll assume we are also using PETSc's ParMetis.
          AS_IF([test "x$build_metis" = "xyes"],
                [
                  PARMETIS_INCLUDE="-I\$(top_srcdir)/contrib/parmetis/include"
                  PARMETIS_LIB="\$(EXTERNAL_LIBDIR)/libparmetis\$(libext)"
                ])
          AC_DEFINE(HAVE_PARMETIS, 1, [Flag indicating whether the library will be compiled with Parmetis support])
          AC_MSG_RESULT(<<< Configuring library with Parmetis support >>>)
        ],
        [
          PARMETIS_INCLUDE=""
          PARMETIS_LIB=""
          enableparmetis=no
        ])

  AC_SUBST(PARMETIS_INCLUDE)
  AC_SUBST(PARMETIS_LIB)
])
