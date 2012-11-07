# ------------------------------------------------------------------------------
# Check to see if the compiler can compile a test program using
# std::thread
# ------------------------------------------------------------------------------
AC_DEFUN([ACX_STD_THREAD],
[AC_CACHE_CHECK(whether the compiler supports std::thread,
ac_cv_cxx_thread,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <thread>],
[
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
# Choose the best unordered_map implementation available
# ------------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_THREAD],
[
ACX_STD_THREAD([])
# ACX_BOOST_THREAD([],[]))
])