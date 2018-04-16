dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_STD_UNORDERED_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_map,
ac_cv_cxx_unordered_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([@%:@include <unordered_map>],
[
  std::unordered_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_unordered_map=yes, ac_cv_cxx_unordered_map=no)
 AC_LANG_RESTORE
])

AS_IF([test "$ac_cv_cxx_unordered_map" = "yes"],
      [
        AC_DEFINE(HAVE_STD_UNORDERED_MAP,1,[define if the compiler supports std::unordered_map])
        AC_DEFINE(BEST_UNORDERED_MAP,std::unordered_map,[definition of the final detected unordered_map type])
        AC_DEFINE(INCLUDE_UNORDERED_MAP,<unordered_map>,[header file for the final detected unordered_map type])
      ],
      [AC_MSG_ERROR([libMesh requires a working std::unordered_map implementation])])
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_multimap
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_STD_UNORDERED_MULTIMAP],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_multimap,
ac_cv_cxx_unordered_multimap,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([@%:@include <unordered_map>],
[
  std::unordered_multimap<int, int> m;
  m.insert(std::make_pair(1, 2));
  m.insert(std::make_pair(1, 3));
  if (m.size() != 2) return 1;
],
 ac_cv_cxx_unordered_multimap=yes, ac_cv_cxx_unordered_multimap=no)
 AC_LANG_RESTORE
])

AS_IF([test "$ac_cv_cxx_unordered_multimap" = "yes"],
      [
        AC_DEFINE(HAVE_STD_UNORDERED_MULTIMAP,1,[define if the compiler supports std::unordered_multimap])
        AC_DEFINE(BEST_UNORDERED_MULTIMAP,std::unordered_multimap,[definition of the final detected unordered_multimap type])
        AC_DEFINE(INCLUDE_UNORDERED_MULTIMAP,<unordered_map>,[header file for the final detected unordered_multimap type])
      ],
      [AC_MSG_ERROR([libMesh requires a working std::unordered_multimap implementation])])
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_multiset
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_STD_UNORDERED_MULTISET],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_multiset,
ac_cv_cxx_unordered_multiset,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([@%:@include <unordered_set>],
[
  std::unordered_multiset<int> s;
  s.insert(1);
  s.insert(1);
  if (s.size() != 2) return 1;
],
 ac_cv_cxx_unordered_multiset=yes,
 ac_cv_cxx_unordered_multiset=no)
 AC_LANG_RESTORE
])

AS_IF([test "$ac_cv_cxx_unordered_multiset" = "yes"],
      [
        AC_DEFINE(HAVE_STD_UNORDERED_MULTISET,1,[define if the compiler supports std::unordered_multiset])
        AC_DEFINE(BEST_UNORDERED_MULTISET,std::unordered_multiset,[definition of the final detected unordered_multiset type])
        AC_DEFINE(INCLUDE_UNORDERED_MULTISET,<unordered_set>,[header file for the final detected unordered_multiset type])
      ],
      [AC_MSG_ERROR([libMesh requires a working std::unordered_multiset implementation])])
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_STD_UNORDERED_SET],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_set,
ac_cv_cxx_unordered_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([@%:@include <unordered_set>],
[
  std::unordered_set<int> m;
  m.insert(1);
  m.insert(2);
],
 ac_cv_cxx_unordered_set=yes, ac_cv_cxx_unordered_set=no)
 AC_LANG_RESTORE
])

AS_IF([test "$ac_cv_cxx_unordered_set" = "yes"],
      [
        AC_DEFINE(HAVE_STD_UNORDERED_SET,1,[define if the compiler supports std::unordered_set])
        AC_DEFINE(BEST_UNORDERED_SET,std::unordered_set,[definition of the final detected unordered_set type])
        AC_DEFINE(INCLUDE_UNORDERED_SET,<unordered_set>,[header file for the final detected unordered_set type])
      ],
      [AC_MSG_ERROR([libMesh requires a working std::unordered_set implementation])])
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::hash
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_STD_HASH],
[AC_CACHE_CHECK(whether the compiler supports std::hash,
ac_cv_cxx_hash,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE(
[
  @%:@include <functional>
  @%:@include <string>
],
[
  std::hash<int> m;
  std::size_t hashed = m(12345);
  std::hash<std::string> m2;
  std::size_t hashed2 = m2(std::string("foo"));
],
 ac_cv_cxx_hash=yes, ac_cv_cxx_hash=no)
 AC_LANG_RESTORE
])

AS_IF([test "$ac_cv_cxx_hash" = "yes"],
      [
        AC_DEFINE(HAVE_STD_HASH,1,[define if the compiler supports std::hash])
        AC_DEFINE(BEST_HASH,std::hash,[definition of the final detected hash type])
        AC_DEFINE(INCLUDE_HASH,<functional>,[header file for the final detected hash type])
      ],
      [AC_MSG_ERROR([libMesh requires a working std::hash implementation])])
])




dnl ----------------------------------------------------------------------------
dnl Choose the best unordered_map implementation available
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_UNORDERED_MAP],
[
ACX_STD_UNORDERED_MAP
])



dnl ----------------------------------------------------------------------------
dnl Choose the best unordered_multimap implementation available
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_UNORDERED_MULTIMAP],
[
ACX_STD_UNORDERED_MULTIMAP
])



dnl ----------------------------------------------------------------------------
dnl Choose the best unordered_set implementation available
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_UNORDERED_SET],
[
ACX_STD_UNORDERED_SET
])



dnl ----------------------------------------------------------------------------
dnl Choose the best unordered_multiset implementation available
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_UNORDERED_MULTISET],
[
ACX_STD_UNORDERED_MULTISET
])



dnl ----------------------------------------------------------------------------
dnl Choose the best hash implementation available
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_HASH],
[
ACX_STD_HASH
])
