dnl -------------------------------------------------------------
dnl Metis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_METIS],
[
  AC_ARG_ENABLE(metis,
                AS_HELP_STRING([--disable-metis],
                               [build without Metis graph partitioning support]),
                [AS_CASE("${enableval}",
                         [yes], [enablemetis=yes],
                         [no],  [enablemetis=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-metis)])],
                [enablemetis=$enableoptional])

  AC_ARG_WITH(metis,
             AS_HELP_STRING([--with-metis=<internal,PETSc,/some/libdir>],
                            [internal: build from contrib, PETSc: rely on PETSc]),
             [AS_CASE("${withval}",
                      [internal], [build_metis=yes],
                      [PETSc],    [build_metis=petsc],
                      [METIS_LIB="-L${withval} -lmetis"
                       build_metis=no])
              enablemetis=yes],
             [build_metis=yes])


  AC_ARG_WITH(metis-include,
             AS_HELP_STRING([--with-metis-include=</some/includedir>]),
             [METIS_INCLUDE="-I${withval}"
              enablemetis=yes
              build_metis=no],
             [])

  dnl If PETSc has its own METIS, default to using that one regardless
  dnl of what the user specified (if anything) in --with-metis.
  AS_IF([test "x$petsc_have_metis" != "x" && test $petsc_have_metis -gt 0],
        [AC_MSG_RESULT(<<< Using PETSc Metis support to avoid PETSc conflict>>>)
         build_metis=petsc])

  dnl Conversely, if:
  dnl .) METIS is enabled in libmesh,
  dnl .) PETSc does not have a METIS or we aren't using PETSc, and
  dnl .) build_metis=petsc because user said --with-metis=PETSc,
  dnl then we need to make sure that libmesh builds its own METIS!
  AS_IF([test "$enablemetis" = "yes" && test "$build_metis" = "petsc"],
        [
          AS_IF([test "x$petsc_have_metis" = "x0" || test "$enablepetsc" = "no"],
                [build_metis=yes])
        ])

  dnl The METIS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  AS_IF([test "$enablemetis" = "yes"],
        [
          dnl Only put Metis contrib directories into the compiler include paths if we're building our own Metis.
          AS_IF([test "$build_metis" = "yes"],
                [
                  METIS_INCLUDE="-I\$(top_srcdir)/contrib/metis/include"
                  METIS_LIB="" # contrib Metis gets lumped into libcontrib
                ])
          AC_DEFINE(HAVE_METIS, 1, [Flag indicating whether the library will be compiled with Metis support])
          AC_MSG_RESULT(<<< Configuring library with Metis support >>>)

          dnl look for thread-local storage
          AX_TLS
        ],
        [
          METIS_INCLUDE=""
          METIS_LIB=""
          enablemetis=no
        ])
  AM_CONDITIONAL(BUILD_METIS, test x$build_metis = xyes)

  AC_SUBST(METIS_INCLUDE)
  AC_SUBST(METIS_LIB)
])
