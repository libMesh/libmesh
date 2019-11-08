dnl -------------------------------------------------------------
dnl Determine the C++ compiler in use. Return the name and possibly
dnl version of this compiler in GXX_VERSION.
dnl
dnl Usage: DETERMINE_CXX_BRAND
dnl
dnl -------------------------------------------------------------
AC_DEFUN([DETERMINE_CXX_BRAND],
[
  dnl Set this flag once the compiler "brand" has been detected to skip remaining tests.
  compiler_brand_detected=no

  dnl First check for gcc version, prevents intel's icc from
  dnl pretending to be gcc
  REAL_GXX=`($CXX -v 2>&1) | grep "gcc version"`

  dnl Do not allow Intel to masquerade as a "real" GCC.
  is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'Compiler'`"
  AS_IF([test "x$is_intel_icc" != "x"],
        [REAL_GXX=""])

  AS_IF([test "$GXX" = "yes" && test "x$REAL_GXX" != "x"],
        [
          dnl find out the right version
          GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "gcc version"`

          AS_CASE("$GXX_VERSION_STRING",
                  [*gcc\ version\ 7.*], [AC_MSG_RESULT(<<< C++ compiler is gcc-7.x >>>)
                                         GXX_VERSION=gcc7],
                  [*gcc\ version\ 6.*], [AC_MSG_RESULT(<<< C++ compiler is gcc-6.x >>>)
                                         GXX_VERSION=gcc6],
                  [*gcc\ version\ 5.*], [AC_MSG_RESULT(<<< C++ compiler is gcc-5.x >>>)
                                         GXX_VERSION=gcc5],
                  [*4.9.*], [AC_MSG_RESULT(<<< C++ compiler is gcc-4.9 >>>)
                             GXX_VERSION=gcc4.9],
                  [*4.8.*], [AC_MSG_RESULT(<<< C++ compiler is gcc-4.8 >>>)
                             GXX_VERSION=gcc4.8],
                  [*4.7.*], [AC_MSG_RESULT(<<< C++ compiler is gcc-4.7 >>>)
                             GXX_VERSION=gcc4.7],
                  [*4.6.*], [AC_MSG_RESULT(<<< C++ compiler is gcc-4.6 >>>)
                             GXX_VERSION=gcc4.6],
                  [AC_MSG_RESULT(<<< C++ compiler is unknown but accepted gcc version >>>)
                   GXX_VERSION=gcc-other])

          dnl Detection was successful, so set the flag.
          compiler_brand_detected=yes
        ])

  dnl Clang/LLVM C++?
  AS_IF([test "x$compiler_brand_detected" = "xno"],
        [
          clang_version="`($CXX --version 2>&1)`"
          is_clang="`echo $clang_version | grep 'clang'`"

          AS_IF([test "x$is_clang" != "x"],
                [
                  dnl Detect if clang is the version built by
                  dnl Apple, because then the version number means
                  dnl something different...
                  is_apple_clang="`echo $clang_version | grep 'Apple'`"
                  clang_vendor="clang"
                  AS_IF([test "x$is_apple_clang" != "x"],
                        [clang_vendor="Apple clang"])

                  dnl If we have perl, we can also pull out the clang version number using regexes.
                  dnl Note that @S|@ is a quadrigraph for the dollar sign.
                  clang_major_minor=unknown

                  AS_IF([test "x$PERL" != "x"],
                        [
                          clang_major_minor=`echo $clang_version | $PERL -ne 'print @S|@1 if /version\s(\d+\.\d+)/'`
                          AS_IF([test "x$clang_major_minor" = "x"],
                                [clang_major_minor=unknown])
                        ])

                  AC_MSG_RESULT([<<< C++ compiler is ${clang_vendor}, version ${clang_major_minor} >>>])
                  GXX_VERSION=clang
                  compiler_brand_detected=yes
                ])
        ])

  dnl Intel's ICC C++ compiler?
  AS_IF([test "x$compiler_brand_detected" = "xno"],
        [
          is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'Compiler'`"
          AS_IF([test "x$is_intel_icc" != "x"],
                [
                  GXX_VERSION_STRING="`($CXX -V 2>&1) | grep 'Version '`"
                  AS_CASE("$GXX_VERSION_STRING",
                  [*19.*], [AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 19 >>>)
                            GXX_VERSION=intel_icc_v19.x],
                  [*18.*], [AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 18 >>>)
                            GXX_VERSION=intel_icc_v18.x],
                  [*17.*], [AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 17 >>>)
                            GXX_VERSION=intel_icc_v17.x],
                  [*16.*], [AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 16 >>>)
                            GXX_VERSION=intel_icc_v16.x],
                  [*15.*], [AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 15 >>>)
                            GXX_VERSION=intel_icc_v15.x],
                  [*14.*], [AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 14 >>>)
                            GXX_VERSION=intel_icc_v14.x],
                  [*13.*], [AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 13 >>>)
                            GXX_VERSION=intel_icc_v13.x],
                  [AC_MSG_ERROR(Unsupported Intel compiler detected)])
                  compiler_brand_detected=yes
                ])
        ])

  dnl Check for IBM xlC. Depending on environment
  dnl variables, moon position, and other reasons unknown to me, this
  dnl compiler can display different names in the first line of output, so
  dnl check various possibilities.  Calling xlC with no arguments displays
  dnl the man page.  Grepping for case-sensitive xlc is not enough if the
  dnl user wants xlC, so we use case-insensitive grep instead.
  AS_IF([test "x$compiler_brand_detected" = "xno"],
        [
          is_ibm_xlc="`($CXX 2>&1) | egrep -i 'xlc'`"
          AS_IF([test "x$is_ibm_xlc" != "x"],
                [
                  AC_MSG_RESULT(<<< C++ compiler is IBM xlC >>>)
                  GXX_VERSION=ibm_xlc
                  compiler_brand_detected=yes
                ])
        ])

  dnl Cray C++?
  AS_IF([test "x$compiler_brand_detected" = "xno"],
        [
          is_cray_cc="`($CXX -V 2>&1) | grep 'Cray '`"
          AS_IF([test "x$is_cray_cc" != "x"],
                [
                  AC_MSG_RESULT(<<< C++ compiler is Cray C++ >>>)
                  GXX_VERSION=cray_cc
                  compiler_brand_detected=yes
                ])
        ])

  dnl Portland Group C++?
  AS_IF([test "x$compiler_brand_detected" = "xno"],
        [
          is_pgcc="`($CXX -V 2>&1) | grep 'Portland Group\|PGI'`"
          AS_IF([test "x$is_pgcc" != "x"],
          [
            AC_MSG_RESULT(<<< C++ compiler is Portland Group C++ >>>)
            GXX_VERSION=portland_group
            compiler_brand_detected=yes
          ])
        ])

  dnl Could not recognize the compiler. Warn the user and continue.
  AS_IF([test "x$compiler_brand_detected" = "xno"],
        [
          AC_MSG_RESULT( WARNING:)
          AC_MSG_RESULT( >>> Unrecognized compiler: "$CXX" <<<)
          AC_MSG_RESULT( You will likely need to modify)
          AC_MSG_RESULT( Make.common directly to specify)
          AC_MSG_RESULT( proper compiler flags)
          GXX_VERSION=unknown
        ])
])
