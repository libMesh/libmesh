AC_DEFUN([AC_HAVE_FEEXCEPT],
[
# This test can apparently generate a false positive on Mac systems.
#AC_CHECK_LIB([m],[feenableexcept],[have_feenableexcept=yes AC_DEFINE([HAVE_FEENABLEEXCEPT],[1],[libm includes feenableexcept])])
#AC_CHECK_LIB([m],[fedisableexcept],[have_fedisableexcept=yes AC_DEFINE([HAVE_FEDISABLEEXCEPT],[1],[libm includes fedisableexcept])])

# So try to compile an actual program with these function calls.
AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([@%:@include <fenv.h>],
               [feenableexcept(FE_DIVBYZERO | FE_INVALID);])],
            [ac_cv_have_feenableexcept=yes],
            [ac_cv_have_feenableexcept=no])
if test "$ac_cv_have_feenableexcept" = yes; then
  AC_DEFINE([HAVE_FEENABLEEXCEPT],[1],[define if the compiler supports feenableexcept])
fi

# So try to compile an actual program with these function calls.
AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([@%:@include <fenv.h>],
               [fedisableexcept(FE_DIVBYZERO | FE_INVALID);])],
            [ac_cv_have_fedisableexcept=yes],
            [ac_cv_have_fedisableexcept=no])
if test "$ac_cv_have_fedisableexcept" = yes; then
  AC_DEFINE([HAVE_FEDISABLEEXCEPT],[1],[define if the compiler supports fedisableexcept])
fi

])
