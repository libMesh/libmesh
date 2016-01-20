# --------------------------------------------------------------
# Test for compiler with broken errno_t.  This is currently only known
# to affect GCC 5.2.0 on OSX Yosemite and above.  It seems to be
# possible to workaround the issue by declaring an errno_t type of our
# own as necessary.
# --------------------------------------------------------------
AC_DEFUN([CHECK_FOR_BROKEN_ERRNO_T],
[
  # Try to compile a test program that exhibits the error.
  AC_LANG_PUSH([C++])

  AC_MSG_CHECKING([if compiler has working errno_t declaration])

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

  # If errno_t does not work, verify that "typedef int errno_t;" fixes it.
  if (test x$errno_t_works = xno); then
    AC_MSG_CHECKING([if errno_t typedef fixes the issue])

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    namespace Foo {
    @%:@include <errno.h>
    }
    typedef int errno_t; // the fix
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

  # If the compiler has the errno_t issue and the typedef fixes it, set the #define
  if (test "x$errno_t_works" = "xno" -a "x$typedef_errno_t_fixes_issue" = "xyes"); then
    AC_DEFINE(COMPILER_HAS_BROKEN_ERRNO_T, 1, [define if the compiler has a broken errno_t])
  fi

  # On the other hand, if the compiler has the errno_t issue and the typedef *doesn't* fix it, error out.
  if (test "x$errno_t_works" = "xno" -a "x$typedef_errno_t_fixes_issue" = "x"); then
    AC_MSG_ERROR([Compiler has the errno_t issue, but typedef does not fix it.  This
                  compiler will most likely not be able to build libmesh correctly.])
  fi
])
