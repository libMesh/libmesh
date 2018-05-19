dnl -------------------------------------------------------------
dnl AC_CXX_HAVE_SSTREAM
dnl -------------------------------------------------------------
AC_DEFUN([AC_CXX_HAVE_SSTREAM],
[AC_CACHE_CHECK(whether the compiler has stringstream,
ac_cv_cxx_have_sstream,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE(
[
@%:@include <sstream>
@%:@ifdef HAVE_NAMESPACES
using namespace std;
@%:@endif
],
[
  stringstream message; 
  message << "Hello"; 
  return 0;
],
 ac_cv_cxx_have_sstream=yes, ac_cv_cxx_have_sstream=no)
 AC_LANG_RESTORE
])
AS_IF([test "x$ac_cv_cxx_have_sstream" = "xyes"],
      [
        AC_DEFINE(HAVE_SSTREAM,,[define if the compiler has the sstream header])
        AC_DEFINE(HAVE_STRINGSTREAM,,[define if the compiler has stringstream functionality])
      ],
      [
        dnl Some compilers implement strstream instead of sstream.  Check for it.
        AC_CXX_HAVE_STRSTREAM
      ])
])
