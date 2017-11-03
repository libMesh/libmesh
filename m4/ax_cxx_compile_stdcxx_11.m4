# ============================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx_11.html
# ============================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX_11([ext|noext],[mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the C++11
#   standard; if necessary, add switches to CXXFLAGS to enable support.
#
#   The first argument, if specified, indicates whether you insist on an
#   extended mode (e.g. -std=gnu++11) or a strict conformance mode (e.g.
#   -std=c++11).  If neither is specified, you get whatever works, with
#   preference for an extended mode.
#
#   The second argument, if specified 'mandatory' or if left unspecified,
#   indicates that baseline C++11 support is required and that the macro
#   should error out if no mode with that support is found.  If specified
#   'optional', then configuration proceeds regardless, after defining
#   HAVE_CXX11 if and only if a supporting mode is found.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2014 Alexey Sokolov <sokolov@google.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 4

m4_define([_AX_CXX_COMPILE_STDCXX_11_testbody], [[
  template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);

    auto d = a;
    auto l = [](){};

    // The compiler may evaluate something like:
    // const int val = multiply(10, 10);
    // at compile time.
    constexpr int multiply (int x, int y) { return x * y; }

    // A minimally-conforming C++11 compiler must support alias declarations
    template <typename T>
    using MyCheck = check<T>;
]])

AC_DEFUN([AX_CXX_COMPILE_STDCXX_11], [dnl
  m4_if([$1], [], [],
        [$1], [ext], [],
        [$1], [noext], [],
        [m4_fatal([invalid argument `$1' to AX_CXX_COMPILE_STDCXX_11])])dnl
  m4_if([$2], [], [ax_cxx_compile_cxx11_required=true],
        [$2], [mandatory], [ax_cxx_compile_cxx11_required=true],
        [$2], [optional], [ax_cxx_compile_cxx11_required=false],
        [m4_fatal([invalid second argument `$2' to AX_CXX_COMPILE_STDCXX_11])])
  AC_LANG_PUSH([C++])dnl
  ac_success=no
  AC_CACHE_CHECK(whether $CXX supports C++11 features by default,
  ax_cv_cxx_compile_cxx11,
  [AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_11_testbody])],
    [ax_cv_cxx_compile_cxx11=yes],
    [ax_cv_cxx_compile_cxx11=no])])
  if test x$ax_cv_cxx_compile_cxx11 = xyes; then
    ac_success=yes
  fi

  m4_if([$1], [noext], [], [dnl
  if test x$ac_success = xno; then
    for switch in -std=gnu++11 -std=gnu++0x; do
      cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx11_$switch])
      AC_CACHE_CHECK(whether $CXX supports C++11 features with $switch,
                     $cachevar,
        [ac_save_CXXFLAGS="$CXXFLAGS"
         CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_11_testbody])],
          [eval $cachevar=yes],
          [eval $cachevar=no])
         CXXFLAGS="$ac_save_CXXFLAGS"])
      if eval test x\$$cachevar = xyes; then
        dnl libmesh propagates CXXFLAGS_OPT, CXXFLAGS_DEVEL, and CXXFLAGS_DBG
        dnl to its config script, so make sure they know about the C++11 flag
        dnl we found.  Note that we still determine the proper switch even if
        dnl the user has not enabled C++11, but we only only propagate it to
        dnl the flags if the user has enabled C++11.
        if (test "x$enablecxx11" = "xyes"); then
          CXXFLAGS="$CXXFLAGS $switch"
          CXXFLAGS_OPT="$CXXFLAGS_OPT $switch"
          CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL $switch"
          CXXFLAGS_DBG="$CXXFLAGS_DBG $switch"
        fi
        ac_success=yes
        break
      fi
    done
  fi])

  dnl The only difference between this block of code and the one above appears
  dnl to be the noext vs. ext arguments to m4_if.  I think the intention was to
  dnl allow people to specify whether they want special GNU extensions of the
  dnl C++11 standard vs. "vanilla" C++11, but I'm not sure how important this
  dnl distinction is nowadays...
  m4_if([$1], [ext], [], [dnl
  if test x$ac_success = xno; then
    for switch in -std=c++11 -std=c++0x; do
      cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx11_$switch])
      AC_CACHE_CHECK(whether $CXX supports C++11 features with $switch,
                     $cachevar,
        [ac_save_CXXFLAGS="$CXXFLAGS"
         CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_11_testbody])],
          [eval $cachevar=yes],
          [eval $cachevar=no])
         CXXFLAGS="$ac_save_CXXFLAGS"])
      if eval test x\$$cachevar = xyes; then
        dnl libmesh propagates CXXFLAGS_OPT, CXXFLAGS_DEVEL, and CXXFLAGS_DBG
        dnl to its config script, so make sure they know about the C++11 flag
        dnl we found.  Note that we still determine the proper switch even if
        dnl the user has not enabled C++11, but we only only propagate it to
        dnl the flags if the user has enabled C++11.
        if (test "x$enablecxx11" = "xyes"); then
          CXXFLAGS="$CXXFLAGS $switch"
          CXXFLAGS_OPT="$CXXFLAGS_OPT $switch"
          CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL $switch"
          CXXFLAGS_DBG="$CXXFLAGS_DBG $switch"
        fi
        ac_success=yes
        break
      fi
    done
  fi])
  AC_LANG_POP([C++])

  if test x$ax_cxx_compile_cxx11_required = xtrue; then
    if test x$ac_success = xno; then
      dnl We return error code 2 here, since 0 means success and 1 is
      dnl indistinguishable from other errors.  Ideally, all of the
      dnl AC_MSG_ERROR calls in our m4 files would return a different
      dnl error code, but currently this is not implemented.
      AC_MSG_ERROR([*** A compiler with support for C++11 language features is required.], 2)
    else
      HAVE_CXX11=1
      AC_DEFINE(HAVE_CXX11, 1, [define if the compiler supports basic C++11 syntax])
    fi
  else
    if test x$ac_success = xno; then
      HAVE_CXX11=0
      AC_MSG_NOTICE([No compiler with C++11 support was found])
      dnl We couldn't find a working C++11 compiler, so we need to set
      dnl enablecxx11=no. This way, other C++11 features can potentially
      dnl still be detected, and listed as "have but disabled" by
      dnl configure.
      enablecxx11=no
    else
      dnl Only set LIBMESH_HAVE_CXX11 if the user has actually configured with --enable-cxx11.
      if (test "x$enablecxx11" = "xyes"); then
        HAVE_CXX11=1
        AC_DEFINE(HAVE_CXX11,1,
                  [define if the compiler supports basic C++11 syntax])
      fi
    fi

    AC_SUBST(HAVE_CXX11)
  fi
])
