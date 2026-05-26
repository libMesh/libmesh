dnl -------------------------------------------------------------
dnl Kokkos -- optional, enables the native Kokkos FE math path
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_KOKKOS],
[
  AC_ARG_WITH([kokkos],
    AS_HELP_STRING([--with-kokkos=DIR],
                   [Enable Kokkos support using the installation at DIR]),
    [KOKKOS_DIR="$withval"],
    [KOKKOS_DIR="no"])

   AC_ARG_WITH([kokkos-include],
    AS_HELP_STRING([--with-kokkos-include=DIR],
                   [Enable Kokkos support using the headers in DIR]),
    [KOKKOS_INCLUDE_DIR="$withval"],
    [KOKKOS_INCLUDE_DIR="$KOKKOS_DIR/include"])

   AC_ARG_WITH([kokkos-lib],
    AS_HELP_STRING([--with-kokkos-lib=DIR],
                   [Enable Kokkos support using the libraries in DIR]),
    [KOKKOS_LIB_DIR="$withval"],
    [KOKKOS_LIB_DIR="$KOKKOS_DIR/lib"])

  AC_ARG_WITH([kokkos-backend],
    AS_HELP_STRING([--with-kokkos-backend=BACKEND],
                   [cuda|hip|sycl|openmp|serial (default: auto-detect from KokkosCore_config.h)]),
    [KOKKOS_BACKEND="$withval"], [KOKKOS_BACKEND="auto"])

  dnl Setting --enable-kokkos-required causes an error to be emitted during
  dnl configure if Kokkos is not successfully detected.
  AC_ARG_ENABLE(kokkos-required,
                AS_HELP_STRING([--enable-kokkos-required],
                               [Error if Kokkos is not detected by configure]),
                [AS_CASE("${enableval}",
                         [yes], [kokkosrequired=yes],
                         [no],  [kokkosrequired=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-kokkos-required)])],
                     [kokkosrequired=no])

  dnl Allow the caller (e.g. MOOSE's configure_libmesh.sh) to pre-set the
  dnl Kokkos compiler and flags via environment variables.  If KOKKOS_CXX is
  dnl already set, we skip auto-detection entirely — the caller knows best.
  dnl We use AC_SUBST (not AC_ARG_VAR) so these flags stay scoped to .K
  dnl compilation rules and don't leak into the main CPPFLAGS/CXXFLAGS.

  AS_IF([test "x$KOKKOS_INCLUDE_DIR" != "xno/include" -a "x$KOKKOS_LIB_DIR" != "xno/lib"],
    [
      AC_CHECK_FILE([$KOKKOS_INCLUDE_DIR/Kokkos_Core.hpp],
        [
          enablekokkos=yes
          libmesh_optional_INCLUDES="$libmesh_optional_INCLUDES -I$KOKKOS_INCLUDE_DIR"
          libmesh_optional_LIBS="$libmesh_optional_LIBS -L$KOKKOS_LIB_DIR -lkokkoscore"

          dnl Only auto-detect if KOKKOS_CXX was not pre-set by the caller
          AS_IF([test "x$KOKKOS_CXX" = "x"],
            [
              KOKKOS_CFG="$KOKKOS_INCLUDE_DIR/KokkosCore_config.h"

              dnl Auto-detect backend
              AS_IF([test "x$KOKKOS_BACKEND" = "xauto"],
                [
                  AS_IF([test -r "$KOKKOS_CFG"],
                    [
                      AS_IF([grep -q 'KOKKOS_ENABLE_CUDA' "$KOKKOS_CFG"],
                        [KOKKOS_BACKEND=cuda],
                        [AS_IF([grep -q 'KOKKOS_ENABLE_HIP' "$KOKKOS_CFG"],
                          [KOKKOS_BACKEND=hip],
                          [AS_IF([grep -q 'KOKKOS_ENABLE_SYCL' "$KOKKOS_CFG"],
                            [KOKKOS_BACKEND=sycl],
                            [AS_IF([grep -q 'KOKKOS_ENABLE_OPENMP' "$KOKKOS_CFG"],
                              [KOKKOS_BACKEND=openmp],
                              [KOKKOS_BACKEND=serial])])])])
                    ],
                    [KOKKOS_BACKEND=serial])
                ])

              AC_MSG_RESULT([Kokkos backend: $KOKKOS_BACKEND])

              dnl Check if Kokkos was built with OpenMP
              have_kokkos_openmp=no
              AS_IF([test -r "$KOKKOS_CFG"],
                [AS_IF([grep -q 'KOKKOS_ENABLE_OPENMP' "$KOKKOS_CFG"],
                  [have_kokkos_openmp=yes])])

              case "$KOKKOS_BACKEND" in
                cuda)
                  AC_PATH_PROG([NVCC],[nvcc],[no],[$PATH])
                  AS_IF([test "x$NVCC" = "xno"],
                    [AC_MSG_ERROR([nvcc not found but Kokkos CUDA backend requested])])
                  kokkos_cuda_arch=`sed -n 's/^#define PETSC_PKG_CUDA_MIN_ARCH //p' "${PETSC_DIR}/include/petscconf.h" "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h" 2>/dev/null | head -n 1`
                  AS_IF([test "x$kokkos_cuda_arch" = "x" && test -r "$KOKKOS_CFG"],
                    [kokkos_cuda_arch=`sed -n \
                      -e 's/^#define KOKKOS_ARCH_BLACKWELL120.*/120/p' \
                      -e 's/^#define KOKKOS_ARCH_BLACKWELL100.*/100/p' \
                      -e 's/^#define KOKKOS_ARCH_HOPPER90.*/90/p' \
                      -e 's/^#define KOKKOS_ARCH_ADA89.*/89/p' \
                      -e 's/^#define KOKKOS_ARCH_AMPERE86.*/86/p' \
                      -e 's/^#define KOKKOS_ARCH_AMPERE80.*/80/p' \
                      -e 's/^#define KOKKOS_ARCH_TURING75.*/75/p' \
                      -e 's/^#define KOKKOS_ARCH_VOLTA72.*/72/p' \
                      -e 's/^#define KOKKOS_ARCH_VOLTA70.*/70/p' \
                      -e 's/^#define KOKKOS_ARCH_PASCAL61.*/61/p' \
                      -e 's/^#define KOKKOS_ARCH_PASCAL60.*/60/p' \
                      -e 's/^#define KOKKOS_ARCH_MAXWELL53.*/53/p' \
                      -e 's/^#define KOKKOS_ARCH_MAXWELL52.*/52/p' \
                      -e 's/^#define KOKKOS_ARCH_MAXWELL50.*/50/p' \
                      -e 's/^#define KOKKOS_ARCH_KEPLER37.*/37/p' \
                      -e 's/^#define KOKKOS_ARCH_KEPLER35.*/35/p' \
                      -e 's/^#define KOKKOS_ARCH_KEPLER32.*/32/p' \
                      -e 's/^#define KOKKOS_ARCH_KEPLER30.*/30/p' \
                      "$KOKKOS_CFG" 2>/dev/null | head -n 1`])
                  kokkos_cuda_arch_flag=""
                  AS_IF([test "x$kokkos_cuda_arch" != "x"],
                    [kokkos_cuda_arch_flag="-arch=sm_${kokkos_cuda_arch}"])
                  libmesh_kokkos_host_cxx=""
                  for libmesh_kokkos_cxx_word in $CXX; do
                    case "$libmesh_kokkos_cxx_word" in
                      ccache|sccache|distcc)
                        ;;
                      *)
                        libmesh_kokkos_host_cxx=$libmesh_kokkos_cxx_word
                        break
                        ;;
                    esac
                  done
                  AS_IF([test "x$libmesh_kokkos_host_cxx" = "x"],
                    [libmesh_kokkos_host_cxx=$CXX])
                  KOKKOS_CXX="$NVCC"
                  KOKKOS_CXXFLAGS="$kokkos_cuda_arch_flag --forward-unknown-to-host-compiler --extended-lambda --disable-warnings --x=cu -ccbin=$libmesh_kokkos_host_cxx"
                  KOKKOS_LDFLAGS="--forward-unknown-to-host-compiler $kokkos_cuda_arch_flag -L$KOKKOS_LIB_DIR"
                  AS_IF([test "x$have_kokkos_openmp" = "xyes"],
                    [
                      KOKKOS_CXXFLAGS="$KOKKOS_CXXFLAGS -fopenmp"
                      KOKKOS_LDFLAGS="$KOKKOS_LDFLAGS -fopenmp"
                    ])
                  ;;
                hip)
                  AC_PATH_PROG([HIPCC],[hipcc],[no],[$PATH])
                  AS_IF([test "x$HIPCC" = "xno"],
                    [AC_MSG_ERROR([hipcc not found but Kokkos HIP backend requested])])
                  KOKKOS_CXX="$HIPCC"
                  KOKKOS_LDFLAGS="-L$KOKKOS_LIB_DIR"
                  ;;
                sycl)
                  AC_PATH_PROG([ICPX],[icpx],[no],[$PATH])
                  AS_IF([test "x$ICPX" = "xno"],
                    [AC_MSG_ERROR([icpx not found but Kokkos SYCL backend requested])])
                  KOKKOS_CXX="$ICPX"
                  KOKKOS_CXXFLAGS="-fsycl"
                  KOKKOS_LDFLAGS="-fsycl -L$KOKKOS_LIB_DIR"
                  ;;
                openmp)
                  KOKKOS_CXX="${CXX}"
                  KOKKOS_CXXFLAGS="-fopenmp -x c++"
                  KOKKOS_LDFLAGS="-fopenmp -L$KOKKOS_LIB_DIR"
                  ;;
                serial|*)
                  KOKKOS_CXX="${CXX}"
                  KOKKOS_CXXFLAGS="-x c++"
                  KOKKOS_LDFLAGS="-L$KOKKOS_LIB_DIR"
                  ;;
              esac
            ],
            [AC_MSG_RESULT([Using caller-provided KOKKOS_CXX=$KOKKOS_CXX])])

          dnl Set defaults for any variables not provided by caller or auto-detect
          KOKKOS_CPPFLAGS="${KOKKOS_CPPFLAGS:--DLIBMESH_KOKKOS_COMPILATION -I$KOKKOS_INCLUDE_DIR}"
          KOKKOS_LDFLAGS="${KOKKOS_LDFLAGS:--L$KOKKOS_LIB_DIR}"
          KOKKOS_LIBS="${KOKKOS_LIBS:--lkokkoscore}"

          dnl If KOKKOS_CXX differs from the main compiler, it may not be the MPI
          dnl wrapper and thus may need the wrapper's compile flags explicitly in
          dnl order to find mpi.h.  Query the primary CXX wrapper for compile-time
          dnl flags and fall back to MPI_INCLUDES when probing is unavailable.
          KOKKOS_MPI_CPPFLAGS=""
          AS_IF([test "x$enablempi" = "xyes" && test "x$KOKKOS_CXX" != "x$CXX"],
            [
              AC_MSG_CHECKING([for MPI compile flags usable with KOKKOS_CXX])

              dnl Check for flags from OpenMPI mpicxx
              KOKKOS_MPI_CPPFLAGS=`$CXX -showme:compile 2>/dev/null`

              dnl If we found no OpenMPI results, try MPICH arguments
              AS_IF([test "x$KOKKOS_MPI_CPPFLAGS" = "x"],
                [KOKKOS_MPI_CPPFLAGS=`$CXX -cxx='' -compile_info 2>/dev/null`])

              dnl If we still have nothing, try Intel MPI arguments
              AS_IF([test "x$KOKKOS_MPI_CPPFLAGS" = "x"],
                [KOKKOS_MPI_CPPFLAGS=`$CXX -show 2>/dev/null | sed 's/^[^ ]* //'`])

              dnl Our MPI compiler might be reporting a full compiler command
              dnl rather than just preprocessor flags.  Bare words such as the
              dnl host compiler name are harmless to normal wrappers but nvcc
              dnl treats them as extra input files, which breaks non-linking
              dnl compile rules that also specify -o.  Keep only flags usable
              dnl during preprocessing/compilation here; link flags are handled
              dnl separately by the existing MPI variables.

              STRIPPED_FLAGS=""
              flag_arg_action=""
              for flag in $KOKKOS_MPI_CPPFLAGS; do
                AS_IF([test "x$flag_arg_action" = "xkeep"],
                  [
                    STRIPPED_FLAGS="$STRIPPED_FLAGS $flag"
                    flag_arg_action=""
                  ],
                  [AS_IF([test "x$flag_arg_action" = "xskip"],
                    [flag_arg_action=""],
                    [case "$flag" in
                      -O*|-std=*|-x)
                        ;;
                      -I|-D|-U|-isystem|-iquote|-idirafter|-include|-imacros)
                        STRIPPED_FLAGS="$STRIPPED_FLAGS $flag"
                        flag_arg_action="keep"
                        ;;
                      -I*|-D*|-U*|-pthread|-pthreads|-f*|-m*|-W*)
                        STRIPPED_FLAGS="$STRIPPED_FLAGS $flag"
                        ;;
                      -Xcompiler|-Xcudafe|-Xptxas)
                        STRIPPED_FLAGS="$STRIPPED_FLAGS $flag"
                        flag_arg_action="keep"
                        ;;
                      -L*|-l*|-Wl,*)
                        ;;
                      -*)
                        STRIPPED_FLAGS="$STRIPPED_FLAGS $flag"
                        ;;
                      *)
                        ;;
                    esac])])
              done
              KOKKOS_MPI_CPPFLAGS=$STRIPPED_FLAGS

              dnl If we still have nothing, fall back to the env?
              AS_IF([test "x$KOKKOS_MPI_CPPFLAGS" = "x"],
                [KOKKOS_MPI_CPPFLAGS="$MPI_INCLUDES"])

              dnl Report what we finally do or do not have
              AS_IF([test "x$KOKKOS_MPI_CPPFLAGS" = "x"],
                [AC_MSG_RESULT([not found])],
                [AC_MSG_RESULT([$KOKKOS_MPI_CPPFLAGS])])
            ])

          dnl Fail configure early if the chosen Kokkos compiler/flags/libs cannot
          dnl actually compile and link a minimal Kokkos program.
          AC_MSG_CHECKING([whether the Kokkos compiler configuration works])
          libmesh_save_CXX="$CXX"
          libmesh_save_CPPFLAGS="$CPPFLAGS"
          libmesh_save_CXXFLAGS="$CXXFLAGS"
          libmesh_save_LDFLAGS="$LDFLAGS"
          libmesh_save_LIBS="$LIBS"

          CXX="$KOKKOS_CXX"
          CPPFLAGS="$CPPFLAGS $KOKKOS_CPPFLAGS $KOKKOS_MPI_CPPFLAGS"
          CXXFLAGS="$CXXFLAGS $KOKKOS_CXXFLAGS"
          LDFLAGS="$LDFLAGS $KOKKOS_LDFLAGS"
          LIBS="$LIBS $KOKKOS_LIBS"
          AC_LANG_PUSH([C++])

          AS_IF([test "x$enablempi" = "xyes"],
            [
              LDFLAGS="$LDFLAGS $MPI_LDFLAGS"
              LIBS="$LIBS $MPI_LIBS"
              AC_LINK_IFELSE(
                [AC_LANG_SOURCE([[
  #include <mpi.h>
  #include <Kokkos_Core.hpp>
  int main(int argc, char ** argv)
  {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
  }
  ]])],
                [kokkos_config_works=yes],
                [kokkos_config_works=no])
            ],
            [
              AC_LINK_IFELSE(
                [AC_LANG_SOURCE([[
  #include <Kokkos_Core.hpp>
  int main(int argc, char ** argv)
  {
    Kokkos::initialize(argc, argv);
    Kokkos::finalize();
    return 0;
  }
  ]])],
                [kokkos_config_works=yes],
                [kokkos_config_works=no])
            ])
          AC_LANG_POP([C++])

          CXX="$libmesh_save_CXX"
          CPPFLAGS="$libmesh_save_CPPFLAGS"
          CXXFLAGS="$libmesh_save_CXXFLAGS"
          LDFLAGS="$libmesh_save_LDFLAGS"
          LIBS="$libmesh_save_LIBS"

          AS_IF([test "x$kokkos_config_works" = "xyes"],
            [AC_MSG_RESULT([yes])],
            [AC_MSG_RESULT([no])
             AC_MSG_ERROR([Kokkos compiler/flags failed to compile and link a minimal test program])])

          AC_DEFINE([HAVE_KOKKOS], [1],
                    [Define if Kokkos support is enabled in libMesh])
          AC_MSG_RESULT(<<< Configuring library with Kokkos support >>>)
        ],
        [
          AC_MSG_WARN([Kokkos not found at $KOKKOS_INCLUDE_DIR -- disabling Kokkos FE support])
          enablekokkos=no
        ])
    ],
    [AC_MSG_NOTICE(<<< Configuring library without Kokkos support >>>)
     enablekokkos=no])

  dnl If Kokkos is not enabled, but it *was* required, error out now
  dnl instead of compiling libmesh in an invalid configuration.
  AS_IF([test "$enablekokkos" = "no" && test "$kokkosrequired" = "yes"],
        dnl We return error code 4 here, since 0 means success and 1 is
        dnl indistinguishable from other errors.
        [AC_MSG_ERROR([*** Kokkos was not found, but --enable-kokkos-required was specified.], 4)])

  AC_SUBST([KOKKOS_CXX])
  AC_SUBST([KOKKOS_CPPFLAGS])
  AC_SUBST([KOKKOS_CXXFLAGS])
  AC_SUBST([KOKKOS_LDFLAGS])
  AC_SUBST([KOKKOS_LIBS])
  AC_SUBST([KOKKOS_MPI_CPPFLAGS])
  AM_CONDITIONAL(LIBMESH_ENABLE_KOKKOS, test x$enablekokkos = xyes)
])
