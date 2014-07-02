dnl ----------------------------------------------------------------
dnl Tests for various C++11 features.  These will probably only work
dnl if they are run after the autoconf test that sets -std=c++11.
dnl ----------------------------------------------------------------

AC_DEFUN([LIBMESH_TEST_CXX11_MOVE],
  [
    AC_MSG_CHECKING(for C++11 std::move support)

    AC_LANG_PUSH([C++])

    have_cxx11_move=no

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <utility>
    template <class T>
    void move_swap(T& a, T& b)
    {
      T tmp(std::move(a));
      a = std::move(b);
      b = std::move(tmp);
    }
    ]], [[
        int one = 1, two = 2;
        move_swap(one,two);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_MOVE, 1, [Flag indicating whether compiler supports std::move])
        have_cxx11_move=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    AC_LANG_POP([C++])
    AM_CONDITIONAL(HAVE_CXX11_MOVE, test x$have_cxx11_move == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_SHARED_PTR],
  [
    AC_MSG_CHECKING(for C++11 std::shared_ptr support)

    AC_LANG_PUSH([C++])

    have_cxx11_shared_ptr=no

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <memory>
    ]], [[
        std::shared_ptr<int> p1;
        std::shared_ptr<int> p2 (new int);
        std::shared_ptr<int> p3 (p2);
        p3.reset(new int);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_SHARED_PTR, 1, [Flag indicating whether compiler supports std::shared_ptr])
        have_cxx11_shared_ptr=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    AC_LANG_POP([C++])
    AM_CONDITIONAL(HAVE_CXX11_SHARED_PTR, test x$have_cxx11_shared_ptr == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_UNIQUE_PTR],
  [
    AC_MSG_CHECKING(for C++11 std::unique_ptr support)

    AC_LANG_PUSH([C++])

    have_cxx11_unique_ptr=no

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <iostream>
    @%:@include <memory>
struct Foo
{
  Foo()      { std::cout << "Foo::Foo\n";  }
  ~Foo()     { std::cout << "Foo::~Foo\n"; }
};
    ]], [[
{
  // up now owns a Foo
  std::unique_ptr<Foo> up(new Foo);
} // Foo deleted when up goes out of scope
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_UNIQUE_PTR, 1, [Flag indicating whether compiler supports std::unique_ptr])
        have_cxx11_unique_ptr=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    AC_LANG_POP([C++])
    AM_CONDITIONAL(HAVE_CXX11_UNIQUE_PTR, test x$have_cxx11_unique_ptr == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_REGEX],
  [
    AC_MSG_CHECKING(for C++11 std::regex support)

    AC_LANG_PUSH([C++])

    have_cxx11_regex=no

    dnl We actually have to try and *run* the test program, since
    dnl GCC up to 4.7 will compile this but then is not able to run it.
    AC_RUN_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <regex>
    ]], [[
std::regex integer_regex("(\\+|-)?[[:digit:]]+");
regex_match("abc", integer_regex);
regex_match("123", integer_regex);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_REGEX, 1, [Flag indicating whether compiler supports std::regex])
        have_cxx11_regex=yes
    ],[
        AC_MSG_RESULT(no)
    ],[
        dnl The test program is not run when cross-compiling, so you are supposed to
        dnl provide a "pessimistic" action here.  We'll just assume the compiler does
        dnl not support C++11 regexes in this case.
        AC_MSG_RESULT(no)
    ])

    AC_LANG_POP([C++])
    AM_CONDITIONAL(HAVE_CXX11_REGEX, test x$have_cxx11_regex == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_OVERRIDE],
  [
    AC_MSG_CHECKING(for C++11 override keyword support)

    AC_LANG_PUSH([C++])

    have_cxx11_override=no

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    ]], [[
    struct Base {
    virtual void f() {}
    };
    struct Child : public Base {
    virtual void f() override {}
    };
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_OVERRIDE, 1, [Flag indicating whether compiler supports the override keyword])
        have_cxx11_override=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    AC_LANG_POP([C++])
    AM_CONDITIONAL(HAVE_CXX11_OVERRIDE, test x$have_cxx11_override == xyes)
  ])
