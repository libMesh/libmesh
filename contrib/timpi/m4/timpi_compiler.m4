# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([TIMPI_SET_COMPILERS],
[
  # --------------------------------------------------------------
  # look for a decent C++ compiler or honor --with-cxx=...
  CXX_TRY_LIST="g++ icpc icc pgCC c++"

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

  AS_IF([test "$enablempi" != no],
        [CXX_TRY_LIST="mpicxx mpiCC mpicc $CXX_TRY_LIST"],
        AC_MSG_RESULT([>>> Disabling MPI per user request <<<]))

  AC_ARG_WITH([cxx],
              AS_HELP_STRING([--with-cxx=CXX], [C++ compiler to use]),
              [CXX="$withval"],
              [])

  # --------------------------------------------------------------
  # Determines a C++ compiler to use.  First checks if the variable CXX is
  # already set.  If not, then searches under g++, c++, and other names.
  # --------------------------------------------------------------
  AC_PROG_CXX([$CXX_TRY_LIST])
  # --------------------------------------------------------------

  # --------------------------------------------------------------
  # Figure out which version of a particular compiler, e.g. GCC 4.0,
  # we are using.
  # --------------------------------------------------------------
  DETERMINE_CXX_BRAND
])
