dnl ----------------------------------------------------------------
dnl Tests for various C++17 features.
dnl ----------------------------------------------------------------
AC_DEFUN([LIBMESH_TEST_CXX17_FALLTHROUGH_ATTRIBUTE],
  [
    have_cxx17_fallthrough_attribute=no
    have_double_underscore_attribute_fallthrough=no

    AC_LANG_PUSH([C++])

    # We don't want to use [[fallthrough]] if it generates a compiler warning, so
    # turn on -Werror for the test. This way, if a warning is generated it will be
    # treated as an error and the test will fail. We also try to make sure as many
    # warnings as possible are turned on, so the test has every opportunity to fail.
    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS -Werror -Wall -Wextra -pedantic"

    AC_MSG_CHECKING(for C++17 fallthrough attribute support)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        ]], [[
          int i=1, j;
          switch(i) {
            case 0:
              j=0;
              [[fallthrough]];
            case 1:
              j=1;
              break;
            default:
              j=2;
          }
          (void)j; // Avoid unused variable warning
        ]]
        )],[
      if (test "x$enablecxx11" = "xyes"); then
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX17_FALLTHROUGH_ATTRIBUTE, 1, [Flag indicating whether compiler supports fallthrough attribute])
        have_cxx17_fallthrough_attribute=yes
      else
        AC_MSG_RESULT([yes, but disabled.])
        AC_DEFINE(HAVE_CXX17_FALLTHROUGH_ATTRIBUTE_BUT_DISABLED, 1, [Compiler supports fallthrough attribute, but it is disabled in libmesh])
      fi
    ],[
      AC_MSG_RESULT(no)
    ])

    # If standard [[fallthrough]]; didn't work, maybe __attribute__((fallthrough));
    # will. This definitely works for GCC7 but generates an error on
    # GCC6. Clang 3.9.0 warns about it, which we also don't want, so
    # we are using the -Werror flag when testing again.
    if (test "$have_cxx17_fallthrough_attribute" = no); then
      AC_MSG_CHECKING([for __attribute__ ((fallthrough)) support])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
          ]], [[
            int i=1, j;
            switch(i) {
              case 0:
                j=0;
                __attribute__((fallthrough));
              case 1:
                j=1;
                break;
              default:
                j=2;
            }
          (void)j; // Avoid unused variable warning
          ]]
          )],[
          AC_MSG_RESULT(yes)
          AC_DEFINE(HAVE_DOUBLE_UNDERSCORE_ATTRIBUTE_FALLTHROUGH, 1, [Flag indicating whether compiler supports __attribute__ ((fallthrough))])
          have_double_underscore_attribute_fallthrough=yes
      ],[
        AC_MSG_RESULT(no)
      ])
    fi


    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX17_FALLTHROUGH_ATTRIBUTE, test x$have_cxx17_fallthrough_attribute == xyes)
    AM_CONDITIONAL(HAVE_DOUBLE_UNDERSCORE_ATTRIBUTE_FALLTHROUGH, test x$have_double_underscore_attribute_fallthrough == xyes)
  ])
