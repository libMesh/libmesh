# ------------------------------------------------------------------------------
# Choose the best thread option available, if any
# ------------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_THREAD],
[
  AC_ARG_WITH(thread-model,
              AS_HELP_STRING([--with-thread-model=tbb,pthread,auto,none],[Specify the thread model to use]),
              [case "${withval}" in
                tbb)     requested_thread_model=tbb     ;;
                pthread) requested_thread_model=pthread ;;
                pthreads) requested_thread_model=pthread ;;
                auto)    requested_thread_model=auto    ;;
                none)    requested_thread_model=none    ;;
                *)       AC_MSG_ERROR(bad value ${withval} for --with-thread-model) ;;
              esac],
              requested_thread_model=auto)

  AC_MSG_RESULT([<<< User requested thread model: $requested_thread_model >>>])

  # Set this variable when a threading model is found
  found_thread_model=none

  # Test different threading models, enabling only one.  For auto-detection, the
  # ordering of the tests is:
  # .) TBB
  # .) OpenMP
  # .) Pthread
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
    if (test $enabletbb = no -a "x$requested_thread_model" = "xtbb"); then
      AC_MSG_ERROR([requested threading model, TBB, could not be found.])
    fi
  fi

  # If TBB wasn't selected, try pthreads as long as the user requested it (or auto)
  if (test "$found_thread_model" = none) ; then
    if (test "x$requested_thread_model" = "xpthread" -o "x$requested_thread_model" = "xauto"); then

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
    fi
  fi

  # OpenMP support -- enabled unless the user has requested no
  # threading, or no valid threading models could be found.
  #
  # OpenMP can be enabled independently of whether we are using the
  # TBB or pthread threading models.  The pthread parallel_for() and
  # parallel_reduce() implementations use openmp pragmas directly,
  # while the TBB implementation does not.  We do the test here simply
  # because it makes sense to keep all the threading stuff together
  # rather than spreading it out across different m4 files.
  if (test "$found_thread_model" != none) ; then
    AC_ARG_ENABLE(openmp,
                  AS_HELP_STRING([--disable-openmp],
                                 [Build without OpenMP Support]),
                  enableopenmp=$enableval,
                  enableopenmp=yes)

    if (test "$enableopenmp" != no) ; then
      AX_OPENMP([],[enableopenmp=no])

      # The above call only sets the flag for C++
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
        break
      fi
    fi
  fi

  AC_MSG_RESULT([<<< Found thread model: $found_thread_model >>>])
])
