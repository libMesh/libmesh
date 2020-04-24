AC_DEFUN([ACSM_COMPILER_CONTROL_ARGS],
[
# -------------------------------------------------------------------
# MPI -- enabled by default.  Check for it now so we can be somewhat
#                             smart about which compilers to look for
# -------------------------------------------------------------------
AC_ARG_ENABLE(mpi,
              AS_HELP_STRING([--disable-mpi],
                             [build without MPI message passing support]),
              [AS_CASE("${enableval}",
                       [yes], [enablempi=yes],
                       [no],  [enablempi=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --enable-mpi)])],
              [enablempi=yes])

AC_ARG_WITH([mpi],
            AS_HELP_STRING([--with-mpi@<:@=DIR@:>@],
                           [Prefix where MPI is installed (default is MPIHOME and then MPI_HOME)]),
            [
              dnl We have no way of knowing whether the user explicitly enabled mpi
              dnl with --enable-mpi or whether it was set by default, so if the user
              dnl specified --with-mpi=no or --without-mpi we're just going to tell them they've given
              dnl competing options
              AS_IF([test "$enablempi" = yes && test "$withval" = no],
                    [AC_MSG_ERROR([Did you mean to disable MPI? If you really mean it, use the --disable-mpi option instead])]
                    )
              MPI="$withval"
            ],
            [
              AS_ECHO(["note: MPI library path not given..."])
              AS_IF([test x"$MPIHOME" != x],
                    [
                      AS_ECHO(["trying prefix=$MPIHOME"])
                      MPI=$MPIHOME
                    ],
                    [
                      AS_IF([test x"$MPI_HOME" != x],
                            [
                              AS_ECHO(["trying prefix=$MPI_HOME"])
                              MPI=$MPI_HOME
                            ])
                    ])
            ])

# --------------------------------------------------------------
# Allow for disable-optional
# --------------------------------------------------------------
AC_ARG_ENABLE(optional,
              AS_HELP_STRING([--disable-optional],
                             [build without most optional external libraries]),
              [AS_CASE("${enableval}",
                       [yes], [enableoptional=yes],
                       [no],  [enableoptional=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --enable-optional)])],
              [enableoptional=yes])

# ----------------------------------------------------------------------
# PETSc is our usual solver package. We place this configure option here
# because we may use its CXX for our own compiles
# ----------------------------------------------------------------------
AC_ARG_ENABLE(petsc,
              AS_HELP_STRING([--disable-petsc],
                             [build without PETSc iterative solver support]),
              [AS_CASE("${enableval}",
                       [yes], [enablepetsc=yes;enablepetsc_mpi=yes],
                       [no],  [enablepetsc=no;enablepetsc_mpi=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --enable-petsc)])],
              [enablepetsc=$enableoptional;enablepetsc_mpi=$enableoptional])
])
