dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::tr1::unordered_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_TR1_UNORDERED_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::tr1::unordered_map,
ac_cv_cxx_tr1_unordered_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <tr1/unordered_map>],
[
  std::tr1::unordered_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_tr1_unordered_map=yes, ac_cv_cxx_tr1_unordered_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_tr1_unordered_map" = yes; then
  AC_DEFINE(HAVE_TR1_UNORDERED_MAP,1,
            [define if the compiler supports std::tr1::unordered_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::tr1::unordered_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_TR1_UNORDERED_SET],
[AC_CACHE_CHECK(whether the compiler supports std::tr1::unordered_set,
ac_cv_cxx_tr1_unordered_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <tr1/unordered_set>],
[
  std::tr1::unordered_set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_tr1_unordered_set=yes, ac_cv_cxx_tr1_unordered_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_tr1_unordered_set" = yes; then
  AC_DEFINE(HAVE_TR1_UNORDERED_SET,1,
            [define if the compiler supports std::tr1::unordered_set])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_UNORDERED_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_map,
ac_cv_cxx_unordered_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <unordered_map>],
[
  std::unordered_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_unordered_map=yes, ac_cv_cxx_unordered_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_unordered_map" = yes; then
  AC_DEFINE(HAVE_UNORDERED_MAP,1,
            [define if the compiler supports std::unordered_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_UNORDERED_SET],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_set,
ac_cv_cxx_unordered_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <unordered_set>],
[
  std::unordered_set<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_unordered_set=yes, ac_cv_cxx_unordered_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_unordered_set" = yes; then
  AC_DEFINE(HAVE_UNORDERED_SET,1,
            [define if the compiler supports std::unordered_set])
fi
])




dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl __gnu_cxx::hash_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_EXT_HASH_MAP],
[AC_CACHE_CHECK(whether the compiler supports __gnu_cxx::hash_map,
ac_cv_cxx_ext_hash_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <ext/hash_map>],
[
  __gnu_cxx::hash_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_ext_hash_map=yes, ac_cv_cxx_ext_hash_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_ext_hash_map" = yes; then
  AC_DEFINE(HAVE_EXT_HASH_MAP,1,
            [define if the compiler supports __gnu_cxx::hash_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl __gnu_cxx::hash_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_EXT_HASH_SET],
[AC_CACHE_CHECK(whether the compiler supports __gnu_cxx::hash_set,
ac_cv_cxx_ext_hash_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <ext/hash_set>],
[
  __gnu_cxx::hash_set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_ext_hash_set=yes, ac_cv_cxx_ext_hash_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_ext_hash_set" = yes; then
  AC_DEFINE(HAVE_EXT_HASH_SET,1,
            [define if the compiler supports __gnu_cxx::hash_set])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::hash_map (This is unlikely...)
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_HASH_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::hash_map,
ac_cv_cxx_hash_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <hash_map>],
[
  std::hash_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_hash_map=yes, ac_cv_cxx_hash_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_hash_map" = yes; then
  AC_DEFINE(HAVE_HASH_MAP,1,
            [define if the compiler supports std::hash_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::hash_set (This is unlikely...)
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_HASH_SET],
[AC_CACHE_CHECK(whether the compiler supports std::hash_set,
ac_cv_cxx_hash_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <hash_set>],
[
  std::hash_set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_hash_set=yes, ac_cv_cxx_hash_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_hash_set" = yes; then
  AC_DEFINE(HAVE_HASH_SET,1,
            [define if the compiler supports std::hash_set])
fi
])
