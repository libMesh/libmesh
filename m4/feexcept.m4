AC_DEFUN([AC_HAVE_FEEXCEPT],
[
dnl Try to compile a program with feenableexcept
AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([@%:@include <fenv.h>],
               [feenableexcept(FE_DIVBYZERO | FE_INVALID);])],
            [ac_cv_have_feenableexcept=yes],
            [ac_cv_have_feenableexcept=no])

AS_IF([test "x$ac_cv_have_feenableexcept" = "xyes"],
      [AC_DEFINE([HAVE_FEENABLEEXCEPT],[1],[define if the compiler supports feenableexcept])])

dnl Try to compile a program with fedisableexcept
AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([@%:@include <fenv.h>],
               [fedisableexcept(FE_DIVBYZERO | FE_INVALID);])],
            [ac_cv_have_fedisableexcept=yes],
            [ac_cv_have_fedisableexcept=no])

AS_IF([test "x$ac_cv_have_fedisableexcept" = "xyes"],
      [AC_DEFINE([HAVE_FEDISABLEEXCEPT],[1],[define if the compiler supports fedisableexcept])])
])
