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
    AC_LANG_PUSH([C++])

    have_cxx11_shared_ptr=init

    # Save any original value that CXXFLAGS had
    saveCXXFLAGS="$CXXFLAGS"

    # Try compiling the test code in all methods requested by the user
    for method in ${METHODS}; do
        case "${method}" in
            optimized|opt)
              CXXFLAGS="$saveCXXFLAGS $CXXFLAGS_OPT $CPPFLAGS_OPT";;

            debug|dbg)
              CXXFLAGS="$saveCXXFLAGS $CXXFLAGS_DBG $CPPFLAGS_DBG";;

            devel)
              CXXFLAGS="$saveCXXFLAGS $CXXFLAGS_DEVEL $CPPFLAGS_DEVEL";;

            profiling|pro|prof)
              CXXFLAGS="$saveCXXFLAGS $CXXFLAGS_PROF $CPPFLAGS_PROF";;

            oprofile|oprof)
              CXXFLAGS="$saveCXXFLAGS $CXXFLAGS_OPROF $CPPFLAGS_OPROF";;

            *)
            AC_MSG_ERROR([bad value ${method} for --with-methods])
            ;;
        esac

        # If compilation fails for *any* of the methods, we'll disable
        # shared_ptr support for *all* methods.
        if test "x$have_cxx11_shared_ptr" != xno; then

          AC_MSG_CHECKING([for C++11 std::shared_ptr support with ${method} flags])

          AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
          @%:@include <memory>
          ]], [[
              std::shared_ptr<int> p1;
              std::shared_ptr<int> p2 (new int);
              std::shared_ptr<int> p3 (p2);
              p3.reset(new int);
          ]])],[
              have_cxx11_shared_ptr=yes
              AC_MSG_RESULT(yes)
          ],[
              have_cxx11_shared_ptr=no
              AC_MSG_RESULT(no)
          ])

        fi
    done

    # Only set the header file variable if our flag was set to 'yes'.
    if test "x$have_cxx11_shared_ptr" = xyes; then
      AC_DEFINE(HAVE_CXX11_SHARED_PTR, 1, [Flag indicating whether compiler supports std::shared_ptr])
    fi

    # Restore the original flags, whatever they were.
    CXXFLAGS="$saveCXXFLAGS"

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


AC_DEFUN([LIBMESH_TEST_CXX11_INITIALIZER_LIST],
  [
    AC_MSG_CHECKING(for C++11 initializer list support)

    AC_LANG_PUSH([C++])

    have_cxx11_initializer_list=no

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
@%:@include <vector>
@%:@include <string>
    ]], [[
std::vector<std::string> v = { "xyzzy", "plugh", "abracadabra" };
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_INITIALIZER_LIST, 1, [Flag indicating whether compiler supports initializer lists])
        have_cxx11_initializer_list=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    AC_LANG_POP([C++])
    AM_CONDITIONAL(HAVE_CXX11_INITIALIZER_LIST, test x$have_cxx11_initializer_list == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_VARIADIC_TEMPLATES],
  [
    AC_MSG_CHECKING(for C++11 variadic template support)

    AC_LANG_PUSH([C++])

    have_cxx11_variadic_templates=no

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
// Base case
template <typename T>
T sum(T t) { return t; }

// Compute sum of arbitary number of passed parameters.
template <typename T, typename ...P>
T sum(T t, P ...p)
{
  t += sum(p...);
  return t;
}
    ]], [[
sum(1, 2, 3, 4, 5);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_VARIADIC_TEMPLATES, 1, [Flag indicating whether compiler supports variadic templates])
        have_cxx11_variadic_templates=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    AC_LANG_POP([C++])
    AM_CONDITIONAL(HAVE_CXX11_VARIADIC_TEMPLATES, test x$have_cxx11_variadic_templates == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_THREAD],
  [
    AC_MSG_CHECKING(for C++11 <thread> support)

    AC_LANG_PUSH([C++])

    have_cxx11_thread=no

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
@%:@include <thread>
void my_thread_func() {}
    ]], [[
std::thread t(my_thread_func);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_THREAD, 1, [Flag indicating whether compiler supports std::thread])
        have_cxx11_thread=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    AC_LANG_POP([C++])
    AM_CONDITIONAL(HAVE_CXX11_THREAD, test x$have_cxx11_thread == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_TYPE_TRAITS],
  [
    AC_MSG_CHECKING(for C++11 <type_traits> support)

    AC_LANG_PUSH([C++])

    have_cxx11_type_traits=no

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
@%:@include <type_traits>
    ]], [[
bool a = std::is_void<char>::value;
bool b = std::is_integral<char>::value;
bool c = std::is_floating_point<char>::value;
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_TYPE_TRAITS, 1, [Flag indicating whether compiler supports std::thread])
        have_cxx11_type_traits=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    AC_LANG_POP([C++])
    AM_CONDITIONAL(HAVE_CXX11_TYPE_TRAITS, test x$have_cxx11_type_traits == xyes)
  ])
