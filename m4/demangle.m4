dnl ----------------------------------------------------------------------------
dnl check for gcc name demangling function
dnl Copyright © 2004 Neil Ferguson <nferguso@eso.org>
dnl
dnl Copying and distribution of this file, with or without
dnl modification, are permitted in any medium without royalty
dnl provided the copyright notice and this notice are preserved.
dnl ----------------------------------------------------------------------------
AC_DEFUN([AX_CXX_GCC_ABI_DEMANGLE],
[AC_CACHE_CHECK(whether the compiler supports GCC C++ ABI name demangling,
ac_cv_cxx_gcc_abi_demangle,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE(
 [
   @%:@include <typeinfo>
   @%:@include <cxxabi.h>
   @%:@include <string>
   template<typename TYPE> class A {};
 ],
 [
   A<int> instance;
   int status = 0;
   char* c_name = 0;
   c_name = abi::__cxa_demangle(typeid(instance).name(), 0, 0, &status);
   std::string name(c_name);
   return name == "A<int>";
 ],
 ac_cv_cxx_gcc_abi_demangle=yes, ac_cv_cxx_gcc_abi_demangle=no)
 AC_LANG_RESTORE
])
AS_IF([test "x$ac_cv_cxx_gcc_abi_demangle" = "xyes"],
      [AC_DEFINE(HAVE_GCC_ABI_DEMANGLE,1,[define if the compiler supports GCC C++ ABI name demangling])])
])
