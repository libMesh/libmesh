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
  AC_DEFINE(HAVE_STD_THREAD,1,
            [define if the compiler supports std::thread])
  [$1]
else
  false
  [$2]
fi
])

# ------------------------------------------------------------------------------
# Check to see if the compiler can compile a test program using
# pthreads
# ------------------------------------------------------------------------------
AC_DEFUN([ACX_PTHREAD],
[AC_CACHE_CHECK(whether the compiler supports pthreads,
ac_cv_pthread_thread,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([@%:@include <pthread.h>],
[
  pthread_t thread;
  pthread_join(thread, NULL);
],
 ac_cv_pthread_thread=yes, ac_cv_pthread_thread=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_pthread_thread" = yes; then
  AC_DEFINE(HAVE_PTHREAD_THREAD,1,
            [define if the compiler supports pthreads])
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
  AC_DEFINE(HAVE_TBB_CXX_THREAD,1,
            [define if the compiler supports tbb::tbb_thread])
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
ACX_STD_THREAD([cppthreadflavor="std::thread"],
ACX_TBB_STD_THREAD([cppthreadflavor="tbb::tbb_thread"],
[enablecppthreads=no]))
])
