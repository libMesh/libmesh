# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_SET_COMPILERS],
[
  AC_REQUIRE([ACSM_COMPILER_CONTROL_ARGS])
  AC_REQUIRE([ACSM_SCRAPE_PETSC_CONFIGURE])

  # --------------------------------------------------------------
  # look for a decent C compiler or honor --with-cc=...
  CC_TRY_LIST="gcc icc pgcc cc"
  AS_IF([test "$enablempi" != no],
        [CC_TRY_LIST="mpicc $CC_TRY_LIST"])
  AC_ARG_WITH([cc],
              AS_HELP_STRING([--with-cc=CC], [C compiler to use]),
              [CC="$withval"],
              [])

  # --------------------------------------------------------------
  # Determine a C compiler to use.  If CC is not already set, checks for
  # gcc, cc, and other C compilers.  Then sets the CC variable to the result.
  # --------------------------------------------------------------
  AC_PROG_CC([$CC_TRY_LIST])
  # --------------------------------------------------------------


  # --------------------------------------------------------------
  # libMesh itself is not written in any Fortran and does not need
  # a Fortran compiler.  Many optional packages however are and
  # we need the compiler to figure out how to link those libraries
  #
  # note though than on OSX for example the XCode tools provide
  # a 'mpif77' which will be detected below but is actually an
  # empty shell script wrapper.  Then the compiler will fail to
  # make executables and we will wind up with static libraries
  # due to a bizarre chain of events.  So, add support for
  # --disable-fortran
  # --------------------------------------------------------------
  AC_ARG_ENABLE(fortran,
                AS_HELP_STRING([--disable-fortran],
                               [build without Fortran language support]),
                [AS_CASE("${enableval}",
                         [yes], [enablefortran=yes],
                         [no],  [enablefortran=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-fortran)])],
                [enablefortran=yes])

  AS_IF([test "x$enablefortran" = xyes],
        [
          dnl look for a decent F90+ compiler or honor --with-fc=...
          FC_TRY_LIST="gfortran ifort pgf90 xlf95"
          AS_IF([test "$enablempi" != no],
                [FC_TRY_LIST="mpif90 $FC_TRY_LIST"])
          AC_ARG_WITH([fc],
                      AS_HELP_STRING([--with-fc=FC], [Fortran compiler to use]),
                      [FC="$withval"],
                      [])

          dnl --------------------------------------------------------------
          dnl Determine a F90+ compiler to use.
          dnl --------------------------------------------------------------
          AC_PROG_FC([$FC_TRY_LIST])

          AS_IF([test "x$FC" = "x"],
                [AC_MSG_RESULT(>>> No valid Fortran compiler <<<)
                FC=no
                enablefortran=no])

          dnl look for a decent F77 compiler or honor --with-77=...
          F77_TRY_LIST="gfortran g77 ifort f77 xlf frt pgf77 fort77 fl32 af77 f90 xlf90 pgf90 epcf90 f95 fort xlf95 ifc efc pgf95 lf95"
          AS_IF([test "$enablempi" != no],
                [F77_TRY_LIST="mpif77 $F77_TRY_LIST"])
          AC_ARG_WITH([f77],
                      AS_HELP_STRING([--with-f77=F77], [Fortran compiler to use]),
                      [F77="$withval"],
                      [])

          dnl --------------------------------------------------------------
          dnl Determine a F77 compiler to use.
          dnl --------------------------------------------------------------
          AC_PROG_F77([$F77_TRY_LIST])

          AS_IF([test "x$F77" = "x"],
                [AC_MSG_RESULT(>>> No valid Fortran 77 compiler <<<)
                 F77=no
                 enablefortran=no])
        ],
        [
          dnl when --disable-fortran is specified, explicitly set these
          dnl to "no" to instruct libtool not to bother with them.
          AC_MSG_RESULT(>>> Disabling Fortran language support per user request <<<)
          FC=no
          F77=no
        ])

# --------------------------------------------------------------
  # look for a decent C++ compiler or honor --with-cxx=...
  CXX_TRY_LIST="g++ icpc icc pgCC c++"

  AC_ARG_WITH([cxx],
              AS_HELP_STRING([--with-cxx=CXX], [C++ compiler to use]),
              [CXX="$withval"],
              [])

  dnl If we already have CXX either from --with-cxx or from the environment, then there's no sense going
  dnl any further. Moreover if we are not enabling mpi then we don't have to query for mpi compilers
  dnl or for a compiler from PETSc
  AS_IF([test -z "$CXX" && test "$enablempi" != no],
        [
          dnl Did we get --with-mpi=DIR or was there a MPI_HOME or MPIHOME set?
          AS_IF([test x"$MPI" != x],
                [
                  dnl Inspect $MPI/bin
                  AS_IF([test -d "$MPI/bin"],
                        [
                          AC_CHECK_PROGS(LOCAL_CXX, [mpicxx mpiCC mpicc], [], ["$MPI/bin"])
                          AS_IF([test -z "$LOCAL_CXX"],
                                [AS_ECHO(["None of the wrappers we look for exist in $MPI/bin. We will not try to use mpi compiler wrappers"])],
                                [MPI_USING_WRAPPERS=1;CXX="$MPI/bin/$LOCAL_CXX"])
                        ],
                        [AS_ECHO(["An MPI directory was specified, but $MPI/bin does not exist. We will not try to use mpi compiler wrappers"])])
                ],
                [
                  dnl No MPI directory specified. If we have PETSc, let's try to snoop some
                  dnl information from there. We'll use this information further below and in
                  dnl mpi.m4
                  AS_IF([test x"$PETSC_HAVE_MPI" = x1 && test x"$PETSC_CXX" != x],
                        [],
                        dnl PETSc doesn't define a CXX so we'll just try to pull one from the environment
                        [CXX_TRY_LIST="mpicxx mpiCC mpicc $CXX_TRY_LIST"])
                ]) dnl AS_IF([test x"$MPI" != x])
        ])

  dnl See whether we are using PETSC_CXX. Unfortunately PETSC_CXX may
  dnl be prefixed with a PATH, and its not straightforward to strip it off.
  dnl If CXX is not set AC_PROG_CXX will call AC_CHECK_TOOLS which will prefix
  dnl every argument in CXX_TRY_LIST with values in $PATH, so we will
  dnl not find something like PATH/PETSC_PREFIX/mpicxx. The solution
  dnl then is just to set CXX to PETSC_CXX so that AC_CHECK_TOOLS
  dnl never gets called
  AS_IF([test -z "$CXX" && test x"$PETSC_HAVE_MPI" = x1 && test x"$PETSC_CXX" != x],
        [CXX="$PETSC_CXX"])

  dnl If we still don't have a CXX set then we will try to pick one up from CXX_TRY_LIST
  AC_PROG_CXX([$CXX_TRY_LIST])

  # --------------------------------------------------------------
  # Use new ACSM macro, but then define old vars too from it
  # --------------------------------------------------------------
  ACSM_DETERMINE_CXX_BRAND
  REAL_GXX="$ACSM_REAL_GXX"
  GXX_VERSION_STRING="$ACSM_GXX_VERSION_STRING"
  GXX_VERSION="$ACSM_GXX_VERSION"
])



# -------------------------------------------------------------
# Set compiler+linker flags to their default values. They will be
# modified/augmented according to other options in later steps of
# configuration
#
# Use ACSM for most of this, but then set the environment variables
# that our other configure scripts still expect.
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_SET_CXX_FLAGS],
[
  ACSM_SET_CXX_FLAGS

  ACSM_SET_GLIBCXX_DEBUG_FLAGS

  CXXFLAGS_OPT="$ACSM_CXXFLAGS_OPT"
  CXXFLAGS_DEVEL="$ACSM_CXXFLAGS_DEVEL"
  CXXFLAGS_DBG="$ACSM_CXXFLAGS_DBG"

  CPPFLAGS_OPT="$ACSM_CPPFLAGS_OPT"
  CPPFLAGS_DEVEL="$ACSM_CPPFLAGS_DEVEL"
  CPPFLAGS_DBG="$ACSM_CPPFLAGS_DBG"

  ASSEMBLY_FLAGS="$ACSM_ASSEMBLY_FLAGS"
  NODEPRECATEDFLAG="$ACSM_NODEPRECATEDFLAG"
  OPROFILE_FLAGS="$ACSM_OPROFILE_FLAGS"
  PARANOID_FLAGS="$ACSM_PARANOID_FLAGS"
  PROFILING_FLAGS="$ACSM_PROFILING_FLAGS"
  RPATHFLAG="$ACSM_RPATHFLAG"
  WERROR_FLAGS="$ACSM_WERROR_FLAGS"

  LIBS="$LIBS $ACSM_XDR_LIBS"
])
