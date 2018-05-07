# -------------------------------------------------------------
# fparser
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_FPARSER],
[
  AC_ARG_ENABLE([fparser],
                AS_HELP_STRING([--disable-fparser],
                               [build without C++ function parser support]),
                [AS_CASE("${enableval}",
                         [yes], [enablefparser=yes],
                         [no],  [enablefparser=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-fparser)])],
                [enablefparser=$enableoptional])

  AC_ARG_WITH([fparser],
               AS_HELP_STRING([--with-fparser=<release|none|devel>],
                              [Determine which version of the C++ function parser to use]),
               [AS_CASE("${withval}",
                        [release], [enablefparserdevel=no],
                        [devel],   [enablefparserdevel=yes],
                        [none],    [enablefparser=no],
                        [AC_MSG_ERROR(bad value ${withval} for --with-fparser)])],
               [enablefparserdevel=no])


  AS_IF([test "x$enablefparser" = "xyes"],
        [
          AC_ARG_ENABLE(fparser-debugging,
                        AS_HELP_STRING([--enable-fparser-debugging],
                                       [Build fparser with bytecode debugging functions]),
                        [AS_CASE("${enableval}",
                                 [yes], [enablefparserdebugging=yes],
                                 [no],  [enablefparserdebugging=no],
                                 [AC_MSG_ERROR(bad value ${enableval} for --enable-fparser-debugging)])],
                        [enablefparserdebugging=no])

          dnl The FPARSER API is distributed with libmesh, so we don't have to guess
          dnl where it might be installed...
          AC_PROG_MKDIR_P
          AC_PROG_SED
          AC_PROG_YACC

          FPARSER_INCLUDE="-I\$(top_srcdir)/contrib/fparser"
          FPARSER_LIBRARY="\$(EXTERNAL_LIBDIR)/libfparser\$(libext)"
          AC_DEFINE(HAVE_FPARSER, 1, [Flag indicating whether the library will be compiled with FPARSER support])

          AS_IF([test "x$enablefparserdevel" = "xyes"],
                [
                  AC_DEFINE(HAVE_FPARSER_DEVEL, 1, [Flag indicating whether FPARSER will build the full development version])
                  AC_MSG_RESULT(<<< Configuring library with fparser support (development version) >>>)
                ],
                [
                  AC_DEFINE(HAVE_FPARSER_DEVEL, 0, [Flag indicating whether FPARSER will build the full development version])
                  AC_MSG_RESULT(<<< Configuring library with fparser support (release version) >>>)
                ])

          dnl According to the autoconf docs, "the third argument must have no
          dnl side effects except for setting the variable cache-id"
          AC_CACHE_CHECK([for dlopen support], [ac_cv_cxx_dlopen], AX_CXX_DLOPEN)

          dnl JIT requires dlopen, use the result of the AX_CXX_DLOPEN test.
          AS_IF([test "x$ac_cv_cxx_dlopen" = "xyes"],
                [
                  AC_DEFINE(HAVE_FPARSER_JIT, 1, [Flag indicating whether FPARSER will be built with JIT compilation enabled])
                  AC_MSG_RESULT(<<< Configuring library with fparser JIT compilation support >>>)
                  enablefparserjit=yes
                ],
                [
                  AC_MSG_RESULT(<<< dlopen() not found, configuring library without fparser JIT compilation support >>>)
                  enablefparserjit=no
                ])

          dnl This define in libmesh_config.h is used internally in fparser.hh and various source files
          AS_IF([test "x$enablefparserdebugging" = "xyes"],
                [
                  AC_DEFINE(FPARSER_SUPPORT_DEBUGGING, 1, [Enable fparser debugging functions])
                  AC_MSG_RESULT(<<< Configuring library with fparser debugging functions >>>)
                ])
        ],
        [
          FPARSER_INCLUDE=""
          FPARSER_LIBRARY=""
          enablefparser=no
        ])

  AC_SUBST(FPARSER_INCLUDE)
  AC_SUBST(FPARSER_LIBRARY)

  AM_CONDITIONAL(FPARSER_RELEASE,              test x$enablefparserdevel = xno)
  AM_CONDITIONAL(FPARSER_DEVEL,                test x$enablefparserdevel = xyes)
  AM_CONDITIONAL(FPARSER_SUPPORT_DEBUGGING,    test x$enablefparserdebugging = xyes)
  AM_CONDITIONAL(FPARSER_SUPPORT_JIT,    test x$enablefparserjit = xyes)
])
