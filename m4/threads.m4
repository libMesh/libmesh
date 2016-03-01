# ------------------------------------------------------------------------------
# Check to see if the compiler can compile a test program using
# std::thread
# ------------------------------------------------------------------------------
AC_DEFUN([ACX_STD_THREAD],
[AC_CACHE_CHECK(whether the compiler supports std::thread,
ac_cv_cxx_thread,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([@%:@include <thread>],
[
  thread_local int i;
  std::thread t;
  t.join();
],
 ac_cv_cxx_thread=yes, ac_cv_cxx_thread=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_thread" = yes; then
  true
  [$1]
else
  false
  [$2]
fi
])


# ------------------------------------------------------------------------------
# Check to see if the compiler can invoke the TBB to compile a
# test program using tbb::tbb_thread
# ------------------------------------------------------------------------------
AC_DEFUN([ACX_TBB_STD_THREAD],
[AC_CACHE_CHECK(whether TBB supports tbb::tbb_thread,
ac_cv_tbb_cxx_thread,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 saveCXXFLAGS="$CXXFLAGS"
 CXXFLAGS="$TBB_INCLUDE"
 AC_TRY_COMPILE([@%:@include <tbb/tbb_thread.h>],
[
  tbb::tbb_thread t;
  t.join();
],
 ac_cv_tbb_cxx_thread=yes, ac_cv_tbb_cxx_thread=no)
 AC_LANG_RESTORE
]
CXXFLAGS="$saveCXXFLAGS"
)

# if we don't have functioning TBB avoid a false positive
if test "x$enabletbb" = "xno"; then
  ac_cv_tbb_cxx_thread=no
fi

if test "$ac_cv_tbb_cxx_thread" = yes; then
  true
  [$1]
else
  false
  [$2]
fi
])


# ------------------------------------------------------------------------------
# Choose the best thread option available, if any
# ------------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_THREAD],
[
  AC_ARG_WITH(thread-model,
              AS_HELP_STRING([--with-thread-model=tbb,openmp,pthread,auto],[Specify the thread model to use]),
              [case "${withval}" in
                tbb)     requested_thread_model=tbb     ;;
                openmp)  requested_thread_model=openmp  ;;
                pthread) requested_thread_model=pthread ;;
                auto)    requested_thread_model=auto    ;;
                *)       AC_MSG_ERROR(bad value ${withval} for --with-thread-model) ;;
              esac],
              requested_thread_model=auto)

  AC_MSG_RESULT([User requested thread model: $requested_thread_model])

  ACX_STD_THREAD([cppthreadflavor="std::thread"],
                  ACX_TBB_STD_THREAD([cppthreadflavor="tbb::tbb_thread"], [enablecppthreads=no]))
])
