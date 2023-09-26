dnl -------------------------------------------------------------
dnl AC_CXX_HAVE_LOCALE
dnl -------------------------------------------------------------
AC_DEFUN([AC_CXX_HAVE_LOCALE],
[AC_CACHE_CHECK(whether the compiler has locale,
ac_cv_cxx_have_locale,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG([C++])
 AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
 [
   @%:@include <locale>
   @%:@ifdef HAVE_NAMESPACES
   using namespace std;
   @%:@endif
 ],
 [
   locale loc;
   return 0;
 ])],
 ac_cv_cxx_have_locale=yes, ac_cv_cxx_have_locale=no)
 AC_LANG_RESTORE
])
AS_IF([test "x$ac_cv_cxx_have_locale" = "xyes"],
      [AC_DEFINE(HAVE_LOCALE,,[define if the compiler has locale])])
])
