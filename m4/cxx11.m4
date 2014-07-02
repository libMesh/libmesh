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


