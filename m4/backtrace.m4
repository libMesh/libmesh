dnl ----------------------------------------------------------------------------
dnl check for gcc backtrace functions
dnl ----------------------------------------------------------------------------
AC_DEFUN([AX_CXX_GLIBC_BACKTRACE],
[AC_CACHE_CHECK(whether the c++ compiler supports glibc backtrace,
ac_cv_cxx_glibc_backtrace,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([@%:@include <execinfo.h>],
[void *addresses[10];
int size = backtrace(addresses, 10);
char** strings = backtrace_symbols(addresses, size);
],
 ac_cv_cxx_glibc_backtrace=yes, ac_cv_cxx_glibc_backtrace=no)
 AC_LANG_RESTORE
])
AS_IF([test "x$ac_cv_cxx_glibc_backtrace" = "xyes"],
      [AC_DEFINE(HAVE_GLIBC_BACKTRACE,1, [define if the compiler supports glibc backtrace])])
])
