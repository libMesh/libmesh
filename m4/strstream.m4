dnl -------------------------------------------------------------
dnl AC_CXX_HAVE_STRSTREAM
dnl -------------------------------------------------------------
AC_DEFUN([AC_CXX_HAVE_STRSTREAM],
[AC_CACHE_CHECK(whether the compiler has strstream,
ac_cv_cxx_have_strstream,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG([C++])
 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
   @%:@include <strstream>
   @%:@ifdef HAVE_NAMESPACES
   using namespace std;
   @%:@endif
 ],
 [
   strstream message;
   message << "Hello";
   return 0;
 ])],
 ac_cv_cxx_have_strstream=yes, ac_cv_cxx_have_strstream=no)
 AC_LANG_RESTORE
])
AS_IF([test "x$ac_cv_cxx_have_strstream" = "xyes"],
      [
        AC_DEFINE(HAVE_STRSTREAM,,[define if the compiler has the strstream header])
        AC_DEFINE(HAVE_STRINGSTREAM,,[define if the compiler has stringstream functionality])
      ])
])
