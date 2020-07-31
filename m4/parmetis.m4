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

  AC_ARG_WITH(parmetis,
             AS_HELP_STRING([--with-parmetis=<internal,PETSc,/some/libdir>],
                            [internal: build from contrib, PETSc: rely on PETSc]),
             [AS_CASE("${withval}",
                      [internal], [build_parmetis=yes],
                      [PETSc],    [build_parmetis=petsc],
                      [PARMETIS_LIB="-L${withval} -lparmetis"
                       build_parmetis=no])
              enableparmetis=yes],
             [build_parmetis=yes])

  AC_ARG_WITH(parmetis-include,
             AS_HELP_STRING([--with-parmetis-include=</some/includedir>]),
             [PARMETIS_INCLUDE="-I${withval}"
              enableparmetis=yes
              build_parmetis=no],
             [])

  dnl Initialize $petsc_have_* to 0 if not already set. They may be
  dnl unset if the user configured with --disable-petsc
  AS_IF([test "x$petsc_have_metis" = "x"], [petsc_have_metis=0])
  AS_IF([test "x$petsc_have_parmetis" = "x"], [petsc_have_parmetis=0])

  dnl If we're using a PETSc with its own METIS, default to using
  dnl it's ParMETIS too or no ParMETIS at all, regardless of what the
  dnl user specified (if anything) in --with-parmetis.
  AS_IF([test $petsc_have_metis -gt 0],
        [
         AS_IF([test $petsc_have_parmetis -gt 0],
               [AC_MSG_RESULT(<<< Using PETSc ParMETIS support to avoid PETSc conflict>>>)
                build_parmetis=petsc],
               [AC_MSG_RESULT(<<< Disabling ParMETIS support to avoid PETSc conflict>>>)
                enableparmetis=no])
        ])

  dnl Conversely, if:
  dnl .) ParMETIS is enabled in libmesh,
  dnl .) PETSc does not have a ParMETIS or METIS, or we aren't using PETSc, and
  dnl .) build_parmetis=petsc because user said --with-parmetis=PETSc,
  dnl then we need to make sure that libmesh builds its own ParMETIS!
  AS_IF([test "$enableparmetis" = "yes" && test "$build_parmetis" = "petsc"],
        [
          AS_IF([test "x$petsc_have_metis" = "x0" &&
                 test "x$petsc_have_parmetis" = "x0"],
                [build_parmetis=yes])
          AS_IF([test "$enablepetsc" = "no"],
                [build_parmetis=yes])
        ])

  dnl Trump --enable-parmetis with --disable-mpi
  AS_IF([test "x$enablempi" = "xno"], [enableparmetis=no])

  dnl The PARMETIS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  AS_IF([test "x$enableparmetis" = "xyes"],
        [
          AC_DEFINE(HAVE_PARMETIS, 1, [Flag indicating whether the library will be compiled with Parmetis support])

          dnl Only put ParMETIS contrib directories into the compiler include paths if we're building our own.
          AS_IF([test "$build_parmetis" = "yes"],
                [
                  PARMETIS_INCLUDE="-I\$(top_srcdir)/contrib/parmetis/include"
                  PARMETIS_LIB="" # contrib Parmetis gets lumped into libcontrib
                  AC_MSG_RESULT([<<< Configuring library with internal Parmetis support >>>])
                ],
                [test "$build_parmetis" = "petsc"],
                [
                  PARMETIS_INCLUDE=""
                  PARMETIS_LIB=""
                  AC_MSG_RESULT([<<< Configuring library with PETSc Parmetis support >>>])
                ],
                [
                  AC_MSG_RESULT([<<< Configuring library with external Parmetis support >>>])
                ])
        ],
        [
          PARMETIS_INCLUDE=""
          PARMETIS_LIB=""
          enableparmetis=no
          build_parmetis=no
        ])
  AM_CONDITIONAL(BUILD_PARMETIS, test x$build_parmetis = xyes)

  AC_SUBST(PARMETIS_INCLUDE)
  AC_SUBST(PARMETIS_LIB)
])
