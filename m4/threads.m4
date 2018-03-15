# ------------------------------------------------------------------------------
# Choose the best thread option available, if any
# ------------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_THREAD],
[
  AC_ARG_WITH(thread-model,
              AS_HELP_STRING([--with-thread-model=tbb,pthread,openmp,auto,none],[Specify the thread model to use]),
              [case "${withval}" in
                tbb)     requested_thread_model=tbb     ;;
                pthread) requested_thread_model=pthread ;;
                pthreads) requested_thread_model=pthread ;;
                openmp)  requested_thread_model=openmp ;;
                auto)    requested_thread_model=auto    ;;
                none)    requested_thread_model=none    ;;
                *)       AC_MSG_ERROR(bad value ${withval} for --with-thread-model) ;;
              esac],
              requested_thread_model=auto)

  AC_MSG_RESULT([<<< User requested thread model: $requested_thread_model >>>])

  # We *always* detect the compiler's OpenMP support, because libMesh
  # is frequently compiled with other libraries that use OpenMP
  # directives (it can also use a few itself [0]) and if you don't set
  # OMP_NUM_THREADS in the environment or call omp_set_num_threads()
  # programmatically, you get the questionable default of "1 thread
  # per CPU" [1]. This can lead to drastic over-subscription of
  # hardware. For example, an MPI code run with "mpiexec -np 4" on a
  # machine with 4 physical cores will try and spawn up to 16
  # simultaneous threads in this configuration.
  #
  # [0]: The "pthread" threading model (see below) can optionally use
  # OpenMP pragmas for its parallel_for() and parallel_reduce()
  # implementations.
  #
  # [1]: https://gcc.gnu.org/onlinedocs/libgomp/OMP_005fNUM_005fTHREADS.html
  AX_OPENMP([],[enableopenmp=no])

  # The call to AX_OPENMP only sets OPENMP_CXXFLAGS, so if it worked,
  # we also set it for C, Fortran, and all the METHODs.
  if (test "x$OPENMP_CXXFLAGS" != x) ; then
     AC_MSG_RESULT(<<< Configuring library with OpenMP support >>>)
     OPENMP_CFLAGS=$OPENMP_CXXFLAGS
     OPENMP_FFLAGS=$OPENMP_CXXFLAGS
     CXXFLAGS_OPT="$CXXFLAGS_OPT $OPENMP_CXXFLAGS"
     CXXFLAGS_DBG="$CXXFLAGS_DBG $OPENMP_CXXFLAGS"
     CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL $OPENMP_CXXFLAGS"
     CXXFLAGS_PROF="$CXXFLAGS_PROF $OPENMP_CXXFLAGS"
     CXXFLAGS_OPROF="$CXXFLAGS_OPROF $OPENMP_CXXFLAGS"
     CFLAGS_OPT="$CFLAGS_OPT $OPENMP_CFLAGS"
     CFLAGS_DBG="$CFLAGS_DBG $OPENMP_CFLAGS"
     CFLAGS_DEVEL="$CFLAGS_DEVEL $OPENMP_CFLAGS"
     CFLAGS_PROF="$CFLAGS_PROF $OPENMP_CFLAGS"
     CFLAGS_OPROF="$CFLAGS_OPROF $OPENMP_CFLAGS"
     FFLAGS="$FFLAGS $OPENMP_FFLAGS"
     AC_SUBST(OPENMP_CXXFLAGS)
     AC_SUBST(OPENMP_CFLAGS)
     AC_SUBST(OPENMP_FFLAGS)
  fi

  # If the user explicitly requested the openmp threading model and
  # this compiler doesn't support OpenMP for some reason, we throw an
  # error instead of silently continuing.
  if (test "x$requested_thread_model" = "xopenmp"); then
    if (test "x$enableopenmp" = "xno"); then
      AC_MSG_ERROR([requested openmp threading model, but compiler does not support openmp.])
    fi
  fi

  # Set this variable when a threading model is found
  found_thread_model=none

  # First, try pthreads/openmp as long as the user requested it (or auto).
  if (test "x$requested_thread_model" = "xpthread" -o "x$requested_thread_model" = "xauto" -o "x$requested_thread_model" = "xopenmp"); then
    # Let the user explicitly specify --{enable,disable}-pthreads.
    AC_ARG_ENABLE(pthreads,
                  AS_HELP_STRING([--disable-pthreads],
                                 [build without POSIX threading (pthreads) support]),
                  [case "${enableval}" in
                    yes)  enablepthreads=yes ;;
                    no)  enablepthreads=no ;;
                    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-pthreads) ;;
                  esac],
                  [enablepthreads=$enableoptional])

    if (test "$enablepthreads" = yes) ; then
      AX_PTHREAD

      if (test x$ax_pthread_ok = xyes); then
        AC_DEFINE(USING_THREADS, 1,
                  [Flag indicating whether the library shall be compiled to use any particular thread API.])
        AC_MSG_RESULT(<<< Configuring library with pthread support >>>)
        libmesh_optional_INCLUDES="$PTHREAD_CFLAGS $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$PTHREAD_LIBS $libmesh_optional_LIBS"
        found_thread_model=pthread
      else
        enablepthreads=no
      fi
    fi

    # If pthreads were not enabled, but the user requested it, we treat that as an error:
    # we want to alert the user as soon as possible that their requested thread model
    # could not be configured correctly.
    if (test $enablepthreads = no -a "x$requested_thread_model" = "xpthread"); then
      AC_MSG_ERROR([requested threading model, pthreads, could not be found.])
    fi
    if (test $enablepthreads = no -a "x$requested_thread_model" = "xopenmp"); then
      AC_MSG_ERROR([openmp threading model requested, but required pthread support unavailable.])
    fi
  fi

  # Try to configure TBB if the user explicitly requested it, or if we
  # are doing auto detection and pthreads/openmp detection did not
  # succeed.
  if (test "$found_thread_model" = none) ; then
    if (test "x$requested_thread_model" = "xtbb" -o "x$requested_thread_model" = "xauto"); then
      CONFIGURE_TBB

      if (test $enabletbb = yes); then
        libmesh_optional_INCLUDES="$TBB_INCLUDE $libmesh_optional_INCLUDES"
        libmesh_optional_LIBS="$TBB_LIBRARY $libmesh_optional_LIBS"
        found_thread_model=tbb
      fi

      # If TBB was not enabled, but the user requested it, we treat that as an error:
      # we want to alert the user as soon as possible that their requested thread model
      # could not be configured correctly.
      if (test "x$enabletbb" = "xno" -a "x$requested_thread_model" = "xtbb"); then
        AC_MSG_ERROR([requested threading model, TBB, could not be found.])
      fi
    fi
  fi

  AC_MSG_RESULT([<<< Found thread model: $found_thread_model >>>])
])
