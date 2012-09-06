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
  AC_DEFINE(BEST_UNORDERED_MAP,std::tr1::unordered_map,
            [definition of the final detected unordered_map type])
  AC_DEFINE(INCLUDE_UNORDERED_MAP,<tr1/unordered_map>,
            [header file for the final detected unordered_map type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::tr1::unordered_multimap
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_TR1_UNORDERED_MULTIMAP],
[AC_CACHE_CHECK(whether the compiler supports std::tr1::unordered_multimap,
ac_cv_cxx_tr1_unordered_multimap,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <tr1/unordered_map>],
[
  std::tr1::unordered_multimap<int, int> m;
  m.insert(std::make_pair(1, 2));
  m.insert(std::make_pair(1, 3));
  if (m.size() != 2) return 1;
],
 ac_cv_cxx_tr1_unordered_multimap=yes, ac_cv_cxx_tr1_unordered_multimap=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_tr1_unordered_multimap" = yes; then
  AC_DEFINE(HAVE_TR1_UNORDERED_MULTIMAP,1,
            [define if the compiler supports std::tr1::unordered_multimap])
  AC_DEFINE(BEST_UNORDERED_MULTIMAP,std::tr1::unordered_multimap,
            [definition of the final detected unordered_multimap type])
  AC_DEFINE(INCLUDE_UNORDERED_MULTIMAP,<tr1/unordered_map>,
            [header file for the final detected unordered_multimap type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
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
  AC_DEFINE(BEST_UNORDERED_SET,std::tr1::unordered_set,
            [definition of the final detected unordered_set type])
  AC_DEFINE(INCLUDE_UNORDERED_SET,<tr1/unordered_set>,
            [header file for the final detected unordered_set type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::tr1::hash
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_TR1_HASH],
[AC_CACHE_CHECK(whether the compiler supports std::tr1::hash,
ac_cv_cxx_tr1_hash,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <tr1/functional>],
[
  std::tr1::hash<int> m;
  std::size_t hashed = m(12345);
],
 ac_cv_cxx_tr1_hash=yes, ac_cv_cxx_tr1_hash=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_tr1_hash" = yes; then
  AC_DEFINE(HAVE_TR1_HASH,1,
            [define if the compiler supports std::tr1::hash])
  AC_DEFINE(BEST_HASH,std::tr1::hash,
            [definition of the final detected hash type])
  AC_DEFINE(INCLUDE_HASH,<tr1/functional>,
            [header file for the final detected hash type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_STD_UNORDERED_MAP],
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
  AC_DEFINE(HAVE_STD_UNORDERED_MAP,1,
            [define if the compiler supports std::unordered_map])
  AC_DEFINE(BEST_UNORDERED_MAP,std::unordered_map,
            [definition of the final detected unordered_map type])
  AC_DEFINE(INCLUDE_UNORDERED_MAP,<unordered_map>,
            [header file for the final detected unordered_map type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
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
 AC_TRY_COMPILE([#include <unordered_map>],
[
  std::unordered_multimap<int, int> m;
  m.insert(std::make_pair(1, 2));
  m.insert(std::make_pair(1, 3));
  if (m.size() != 2) return 1;
],
 ac_cv_cxx_unordered_multimap=yes, ac_cv_cxx_unordered_multimap=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_unordered_multimap" = yes; then
  AC_DEFINE(HAVE_STD_UNORDERED_MULTIMAP,1,
            [define if the compiler supports std::unordered_multimap])
  AC_DEFINE(BEST_UNORDERED_MULTIMAP,std::unordered_multimap,
            [definition of the final detected unordered_multimap type])
  AC_DEFINE(INCLUDE_UNORDERED_MULTIMAP,<unordered_map>,
            [header file for the final detected unordered_multimap type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
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
 AC_TRY_COMPILE([#include <unordered_set>],
[
  std::unordered_set<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_unordered_set=yes, ac_cv_cxx_unordered_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_unordered_set" = yes; then
  AC_DEFINE(HAVE_STD_UNORDERED_SET,1,
            [define if the compiler supports std::unordered_set])
  AC_DEFINE(BEST_UNORDERED_SET,std::unordered_set,
            [definition of the final detected unordered_set type])
  AC_DEFINE(INCLUDE_UNORDERED_SET,<unordered_set>,
            [header file for the final detected unordered_set type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
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
 AC_TRY_COMPILE([#include <functional>],
[
  std::hash<int> m;
  std::size_t hashed = m(12345);
],
 ac_cv_cxx_hash=yes, ac_cv_cxx_hash=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_hash" = yes; then
  AC_DEFINE(HAVE_STD_HASH,1,
            [define if the compiler supports std::hash])
  AC_DEFINE(BEST_HASH,std::hash,
            [definition of the final detected hash type])
  AC_DEFINE(INCLUDE_HASH,<functional>,
            [header file for the final detected hash type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
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
  AC_DEFINE(BEST_UNORDERED_MAP,__gnu_cxx::hash_map,
            [definition of the final detected unordered_map type])
  AC_DEFINE(INCLUDE_UNORDERED_MAP,<ext/hash_map>,
            [header file for the final detected unordered_map type])
  AC_DEFINE(DEFINE_HASH_STRING,[namespace __gnu_cxx {template<> struct hash<std::string>{size_t operator()(const std::string& s)const{return hash<const char*>()(s.c_str());}};}],
            [workaround for potentially missing hash<string>])
  ac_cv_cxx_hash_string=yes
  [$1]
else
  false
  [$2]
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl __gnu_cxx::hash_multimap
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_EXT_HASH_MULTIMAP],
[AC_CACHE_CHECK(whether the compiler supports __gnu_cxx::hash_multimap,
ac_cv_cxx_ext_hash_multimap,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <ext/hash_map>],
[
  __gnu_cxx::hash_multimap<int, int> m;
  m.insert(std::make_pair(1, 2));
  m.insert(std::make_pair(1, 3));
  if (m.size() != 2) return 1;
],
 ac_cv_cxx_ext_hash_multimap=yes, ac_cv_cxx_ext_hash_multimap=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_ext_hash_multimap" = yes; then
  AC_DEFINE(HAVE_EXT_HASH_MULTIMAP,1,
            [define if the compiler supports __gnu_cxx::hash_multimap])
  AC_DEFINE(BEST_UNORDERED_MULTIMAP,__gnu_cxx::hash_multimap,
            [definition of the final detected unordered_multimap type])
  AC_DEFINE(INCLUDE_UNORDERED_MULTIMAP,<ext/hash_map>,
            [header file for the final detected unordered_multimap type])
  AC_DEFINE(DEFINE_HASH_STRING,[namespace __gnu_cxx {template<> struct hash<std::string>{size_t operator()(const std::string& s)const{return hash<const char*>()(s.c_str());}};}],
            [workaround for potentially missing hash<string>])
  ac_cv_cxx_hash_string=yes
  [$1]
else
  false
  [$2]
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
  AC_DEFINE(BEST_UNORDERED_SET,__gnu_cxx::hash_set,
            [definition of the final detected unordered_set type])
  AC_DEFINE(INCLUDE_UNORDERED_SET,<ext/hash_set>,
            [header file for the final detected unordered_set type])
  AC_DEFINE(DEFINE_HASH_STRING,[namespace __gnu_cxx {template<> struct hash<std::string>{size_t operator()(const std::string& s)const{return hash<const char*>()(s.c_str());}};}],
            [workaround for potentially missing hash<string>])
  ac_cv_cxx_hash_string=yes
  [$1]
else
  false
  [$2]
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl __gnu_cxx::hash
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_EXT_HASH],
[AC_CACHE_CHECK(whether the compiler supports __gnu_cxx::hash,
ac_cv_cxx_ext_hash,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <ext/hash_set>],
[
  __gnu_cxx::hash<int> m;
  std::size_t hashed = m(12345);
],
 ac_cv_cxx_ext_hash=yes, ac_cv_cxx_ext_hash=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_ext_hash" = yes; then
  AC_DEFINE(HAVE_EXT_HASH,1,
            [define if the compiler supports __gnu_cxx::hash])
  AC_DEFINE(BEST_HASH,__gnu_cxx::hash,
            [definition of the final detected hash type])
  AC_DEFINE(INCLUDE_HASH,<ext/hash_set>,
            [header file for the final detected hash type])
  AC_DEFINE(DEFINE_HASH_STRING,[namespace __gnu_cxx {template<> struct hash<std::string>{size_t operator()(const std::string& s)const{return hash<const char*>()(s.c_str());}};}],
            [workaround for potentially missing hash<string>])
  ac_cv_cxx_hash_string=yes
  [$1]
else
  false
  [$2]
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
  AC_DEFINE(BEST_UNORDERED_MAP,std::hash_map,
            [definition of the final detected unordered_map type])
  AC_DEFINE(INCLUDE_UNORDERED_MAP,<hash_map>,
            [header file for the final detected unordered_map type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::hash_multimap (This is unlikely...)
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_HASH_MULTIMAP],
[AC_CACHE_CHECK(whether the compiler supports std::hash_multimap,
ac_cv_cxx_hash_multimap,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <hash_map>],
[
  std::hash_multimap<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_hash_multimap=yes, ac_cv_cxx_hash_multimap=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_hash_multimap" = yes; then
  AC_DEFINE(HAVE_HASH_MULTIMAP,1,
            [define if the compiler supports std::hash_multimap])
  AC_DEFINE(BEST_UNORDERED_MULTIMAP,std::hash_multimap,
            [definition of the final detected unordered_multimap type])
  AC_DEFINE(INCLUDE_UNORDERED_MULTIMAP,<hash_map>,
            [header file for the final detected unordered_multimap type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
])



dnl ----------------------------------------------------------------------------
dnl Make sure the compiler can compile a test program using
dnl std::set as a "fall-back" unordered_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_STD_SET],
[AC_CACHE_CHECK(whether the compiler supports std::set,
ac_cv_cxx_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <set>],
[
  std::set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_set=yes, ac_cv_cxx_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_set" = yes; then
  AC_DEFINE(BEST_UNORDERED_SET,std::set,
            [definition of the final detected unordered_set type])
  AC_DEFINE(INCLUDE_UNORDERED_SET,<set>,
            [header file for the final detected unordered_set type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
])


dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::map as a "fall-back" unordered_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_STD_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::map,
ac_cv_cxx_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <map>],
[
  std::map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_map=yes, ac_cv_cxx_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_map" = yes; then
  AC_DEFINE(BEST_UNORDERED_MAP,std::map,
            [definition of the final detected unordered_map type])
  AC_DEFINE(INCLUDE_UNORDERED_MAP,<map>,
            [header file for the final detected unordered_map type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::multimap as a "fall-back" unordered_multimap
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_STD_MULTIMAP],
[AC_CACHE_CHECK(whether the compiler supports std::multimap,
ac_cv_cxx_multimap,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <map>],
[
  std::multimap<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_multimap=yes, ac_cv_cxx_multimap=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_multimap" = yes; then
  AC_DEFINE(BEST_UNORDERED_MULTIMAP,std::multimap,
            [definition of the final detected unordered_multimap type])
  AC_DEFINE(INCLUDE_UNORDERED_MULTIMAP,<map>,
            [header file for the final detected unordered_multimap type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
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
  AC_DEFINE(BEST_UNORDERED_SET,std::hash_set,
            [definition of the final detected unordered_set type])
  AC_DEFINE(INCLUDE_UNORDERED_SET,<hash_set>,
            [header file for the final detected unordered_set type])
  if test "$ac_cv_cxx_hash_string" != yes; then
    AC_DEFINE(DEFINE_HASH_STRING,,
              [workaround for potentially missing hash<string>])
  fi
  [$1]
else
  false
  [$2]
fi
])



dnl ----------------------------------------------------------------------------
dnl Choose the best unordered_map implementation available
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_UNORDERED_MAP],
[
ACX_STD_UNORDERED_MAP([],
ACX_TR1_UNORDERED_MAP([],
ACX_EXT_HASH_MAP([],
ACX_HASH_MAP([],
ACX_STD_MAP([],[])))))
])



dnl ----------------------------------------------------------------------------
dnl Choose the best unordered_multimap implementation available
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_UNORDERED_MULTIMAP],
[
ACX_STD_UNORDERED_MULTIMAP([],
ACX_TR1_UNORDERED_MULTIMAP([],
ACX_EXT_HASH_MULTIMAP([],
ACX_HASH_MULTIMAP([],
ACX_STD_MULTIMAP([],[])))))
])



dnl ----------------------------------------------------------------------------
dnl Choose the best unordered_set implementation available
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_UNORDERED_SET],
[
ACX_STD_UNORDERED_SET([],
ACX_TR1_UNORDERED_SET([],
ACX_EXT_HASH_SET([],
ACX_HASH_SET([],
ACX_STD_SET([],[])))))
])



dnl ----------------------------------------------------------------------------
dnl Choose the best hash implementation available
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_BEST_HASH],
[
ACX_STD_HASH([],
ACX_TR1_HASH([],
ACX_EXT_HASH([],[])))
])
