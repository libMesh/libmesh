# --------------------------------------------------------------
# Make sure that the safe_bool<T> installation we're using prevents
# derived classes from being compared with operator==.
# --------------------------------------------------------------
AC_DEFUN([CONFIGURE_SAFE_BOOL],
[
  # Try to compile a test program that uses the unique_ptr.hpp header file
  AC_LANG_PUSH([C++])

  # Save any original value that CXXFLAGS had
  saveCXXFLAGS="$CXXFLAGS"

  # Add location of contrib header file to CXXFLAGS
  CXXFLAGS="$saveCXXFLAGS -I$top_srcdir/include/utils";

  AC_MSG_CHECKING([if safe_bool<T> allows boolean comparisons])

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
  @%:@include "safe_bool.h"
  using namespace libMesh;
  class Foo : public safe_bool<Foo>
  {
  public:
    bool boolean_test() const
    {
      // This dummy class just returns true for operator bool()
      return true;
    }
  };
  ]], [[
  {
    Foo foo1;
    if (foo1)
      return 0;
    return 1;
  }
  ]])],[
    safe_bool_comparison_works=yes
    AC_MSG_RESULT(yes)
  ],[
    AC_MSG_RESULT(no)
  ])

  # This test should not compile, so a failure indicates success :-)
  AC_MSG_CHECKING([if safe_bool<T> prevents implicit conversions])

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
  @%:@include "safe_bool.h"
  using namespace libMesh;
  class Foo : public safe_bool<Foo>
  {
  public:
    bool boolean_test() const
    {
      // This dummy class just returns true for operator bool()
      return true;
    }
  };
  ]], [[
  {
    Foo foo1, foo2;
    if (foo1 == foo2)
      return 1;
    return 0;
  }
  ]])],[
    AC_MSG_RESULT(no)
  ],[
    safe_bool_equality_fails=yes
    AC_MSG_RESULT(yes)
  ])

  # Restore the original flags, whatever they were.
  CXXFLAGS="$saveCXXFLAGS"

  AC_LANG_POP([C++])

  # safe_bool<T> must provide at least basic functionality.
  if (test x$safe_bool_comparison_works = xyes); then
    AC_MSG_RESULT([<<< Compiling libmesh with safe_bool<T> support >>>])
  else
    AC_MSG_ERROR([libMesh requires a working safe_bool<T> implementation])
  fi

  # It would be *nice* if safe_bool prevented incorrect code from
  # compiling, but we don't require that.
  if (test x$safe_bool_equality_fails = xno); then
    AC_MSG_WARN([safe_bool<T> does not prevent incorrect equality comparisons])
  fi
])
