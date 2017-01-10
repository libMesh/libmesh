# --------------------------------------------------------------
# Test for compiler which does not support namespace-wrapped inclusion
# of errno.h.  This is currently only known to affect GCC 5.2.0 on OSX
# Yosemite and above.  One can work around the issue by simply
# including errno.h separately, outside of any namespace.
# --------------------------------------------------------------
AC_DEFUN([CHECK_FOR_BROKEN_ERRNO_T],
[
  # Try to compile a test program that exhibits the error.
  AC_LANG_PUSH([C++])

  AC_MSG_CHECKING([if errno.h can be wrapped in namespace])

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
  namespace Foo {
  @%:@include <errno.h>
  }
  @%:@include <cstring>
  ]], [[
  ]])],[
    AC_MSG_RESULT(yes)
  ],[
    errno_t_works=no
    AC_MSG_RESULT(no)
  ])

  # If that test failed, verify that including the header outside of any namespace fixes the issue.
  if (test x$errno_t_works = xno); then
    AC_MSG_CHECKING([if workaround fixes the issue])

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <errno.h> // the fix
    namespace Foo {
    @%:@include <errno.h>
    }
    @%:@include <cstring>
    ]], [[
    ]])],[
      typedef_errno_t_fixes_issue=yes
      AC_MSG_RESULT(yes)
    ],[
      AC_MSG_RESULT(no)
    ])
  fi

  AC_LANG_POP([C++])

  # If the compiler has the issue and the workaround fixes it, set the #define
  if (test "x$errno_t_works" = "xno" -a "x$typedef_errno_t_fixes_issue" = "xyes"); then
    AC_DEFINE(COMPILER_HAS_BROKEN_ERRNO_T, 1, [define if errno.h cannot be included in a namespace])
  fi

  # On the other hand, if the compiler has the issue and the workaround *doesn't* fix it, error out.
  if (test "x$errno_t_works" = "xno" -a "x$typedef_errno_t_fixes_issue" = "x"); then
    AC_MSG_ERROR([Cannot work around errno.h inclusion issue.  This
                  compiler will most likely not be able to build libmesh correctly.])
  fi
])
