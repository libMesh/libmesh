dnl ----------------------------------------------------------------
dnl Tests for various C++11 features.  These will probably only work
dnl if they are run after the autoconf test that sets -std=c++11.
dnl ----------------------------------------------------------------

dnl Test C++11 std::isnan, std::isinf
AC_DEFUN([LIBMESH_TEST_CXX11_ISNAN_ISINF],
  [
    have_cxx11_isnan=no
    have_cxx11_isinf=no

    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_MSG_CHECKING(for C++11 std::isnan)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <cmath>
    ]], [[
    if (std::isnan(0.0))
      return 1;
    ]])],[
        AC_MSG_RESULT(yes)
        have_cxx11_isnan=yes
        dnl AC_DEFINE(HAVE_CXX11_ISNAN, 1, [Flag indicating whether compiler supports std::isnan])
    ],[
        AC_MSG_RESULT(no)
    ])

    AC_MSG_CHECKING(for C++11 std::isinf)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <cmath>
    ]], [[
    if (std::isinf(0.0))
      return 1;
    ]])],[
        AC_MSG_RESULT(yes)
        have_cxx11_isinf=yes
        dnl AC_DEFINE(HAVE_CXX11_ISINF, 1, [Flag indicating whether compiler supports std::isinf])
    ],[
        AC_MSG_RESULT(no)
    ])

    dnl Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])
  ])

dnl Test C++11 std::array
AC_DEFUN([LIBMESH_TEST_CXX11_ARRAY],
  [
    have_cxx11_array=no

    AC_MSG_CHECKING(for C++11 std::array)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <array>
    ]], [[
    std::array<double, 4> a;
    a[0] = 42.0;
    double * begin = a.data();
    ]])],[
        AC_MSG_RESULT(yes)
        have_cxx11_array=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    dnl Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])
  ])

dnl Test C++11 std::vector::data()
AC_DEFUN([LIBMESH_TEST_CXX11_VECTOR_DATA],
  [
    have_cxx11_vector_data=no

    AC_MSG_CHECKING(for C++11 std::vector::data() API)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <vector>
    ]], [[
    std::vector<int> v(10);
    int * begin = v.data();
    ]])],[
        AC_MSG_RESULT(yes)
        have_cxx11_vector_data=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    dnl Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])
  ])

dnl Test C++11 std::iota
AC_DEFUN([LIBMESH_TEST_CXX11_IOTA],
  [
    have_cxx11_iota=no

    AC_MSG_CHECKING(for C++11 std::iota algorithm)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <vector>
    @%:@include <numeric>
    ]], [[
    std::vector<int> v(10);
    std::iota(v.begin(), v.end(), 0);
    ]])],[
        AC_MSG_RESULT(yes)
        have_cxx11_iota=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    dnl Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])
  ])

dnl Test C++11 std::map,set,multimap,multiset iterator-returning erase() APIs.
AC_DEFUN([LIBMESH_TEST_CXX11_CONTAINER_ERASE],
  [
    have_cxx11_container_erase=no

    AC_MSG_CHECKING(for C++11 std container erase() functions returning iterators)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <map>
    @%:@include <set>
    ]], [[
    {
      std::map<int, int> m;
      m.insert(std::make_pair(1,2));
      std::map<int, int>::iterator it = m.erase(m.begin());
    }
    {
      std::set<int> s;
      s.insert(1);
      std::set<int>::iterator it = s.erase(s.begin());
    }
    {
      std::multimap<int, int> m;
      m.insert(std::make_pair(1,2));
      std::multimap<int, int>::iterator it = m.erase(m.begin());
    }
    {
      std::multiset<int> s;
      s.insert(1);
      std::multiset<int>::iterator it = s.erase(s.begin());
    }
    ]])],[
        AC_MSG_RESULT(yes)
        have_cxx11_container_erase=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    dnl Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])
  ])

dnl Test C++11 std::tuple and several related helper functions.
AC_DEFUN([LIBMESH_TEST_CXX11_BEGIN_END],
  [
    have_cxx11_begin_end=no

    AC_MSG_CHECKING(for C++11 std::begin/end support for arrays)
    AC_LANG_PUSH([C++])

    dnl For this and all of the C++ standards tests: Save the original
    dnl CXXFLAGS (if any) before appending the $switch determined by
    dnl AX_CXX_COMPILE_STDCXX_11, and any compiler flags specified by
    dnl the user in the libmesh_CXXFLAGS environment variable, letting
    dnl that override everything else.
    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <iterator>
    @%:@include <algorithm>
    ]], [[
    int array[5] = {3, 1, 5, 2, 4};
    std::sort(std::begin(array), std::end(array));
    ]])],[
        AC_MSG_RESULT(yes)
        have_cxx11_begin_end=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    dnl Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])
  ])

dnl Test C++11 std::tuple and several related helper functions.
AC_DEFUN([LIBMESH_TEST_CXX11_TUPLE],
  [
    have_cxx11_tuple=no

    AC_MSG_CHECKING(for C++11 std::tuple support)
    AC_LANG_PUSH([C++])

    dnl For this and all of the C++ standards tests: Save the original
    dnl CXXFLAGS (if any) before appending the $switch determined by
    dnl AX_CXX_COMPILE_STDCXX_11, and any compiler flags specified by
    dnl the user in the libmesh_CXXFLAGS environment variable, letting
    dnl that override everything else.
    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <tuple>
    @%:@include <map>
    @%:@include <string>
    ]], [[
    // Test std::make_tuple
    std::map<int, std::tuple<double, char, std::string>> students;
    students[0] = std::make_tuple(3.8, 'A', "Lisa Simpson");
    students[1] = std::make_tuple(2.9, 'C', "Milhouse Van Houten");

    // Test templated std::get() method.
    std::get<0>(students[0]);
    std::get<1>(students[0]);
    std::get<2>(students[0]);

    // Test std::tie and std::ignore which are declared in the <tuple> header.
    double gpa1;
    std::string name1;
    std::tie(gpa1, std::ignore, name1) = students[1];

    // Test std::tuple_cat(), which combines together one or more
    // objects, not necessarily tuples, into a single tuple.
    auto t = std::tuple_cat(students[0], students[1], std::make_pair("foo", 42));
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_TUPLE, 1, [Flag indicating whether compiler supports std::move])
        have_cxx11_tuple=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    dnl Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_TUPLE, test x$have_cxx11_tuple == xyes)
  ])

dnl Test C++11 fixed type enumerations.
AC_DEFUN([LIBMESH_TEST_CXX11_FIXED_TYPE_ENUM],
  [
    have_cxx11_fixed_type_enum=no

    AC_MSG_CHECKING(for C++11 fixed type enumeration support)
    AC_LANG_PUSH([C++])

    dnl For this and all of the C++ standards tests: Save the original
    dnl CXXFLAGS (if any) before appending the $switch determined by
    dnl AX_CXX_COMPILE_STDCXX_11, and any compiler flags specified by
    dnl the user in the libmesh_CXXFLAGS environment variable, letting
    dnl that override everything else.
    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    enum Fruit : int {APPLE=0, ORANGE=1};
    ]], [[
    Fruit f = APPLE;
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_FIXED_TYPE_ENUM, 1, [Flag indicating whether compiler supports fixed type enumerations])
        have_cxx11_fixed_type_enum=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    dnl Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_FIXED_TYPE_ENUM, test x$have_cxx11_fixed_type_enum == xyes)
  ])

dnl Test C++11 fixed type enumeration forward declarations.
AC_DEFUN([LIBMESH_TEST_CXX11_FIXED_TYPE_ENUM_FWD],
  [
    have_cxx11_fixed_type_enum_fwd=no

    AC_MSG_CHECKING(for C++11 fixed type enumeration forward declaration support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <iostream>
    enum Foo : int;
    void func(Foo f) { std::cout << f << std::endl; }
    enum Foo : int { FOO0 = 0, FOO1 = 1, FOO2 = 2 };
    ]], [[
    func(FOO0);
    func(FOO1);
    func(FOO2);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_FIXED_TYPE_ENUM_FWD, 1, [Flag indicating whether compiler supports fixed type enumerations])
        have_cxx11_fixed_type_enum_fwd=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    dnl Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_FIXED_TYPE_ENUM_FWD, test x$have_cxx11_fixed_type_enum_fwd == xyes)
  ])

dnl Properly implemented move constructors require rvalue references,
dnl std::move, and noexcept, so this tests for all of those features.
AC_DEFUN([LIBMESH_TEST_CXX11_MOVE_CONSTRUCTORS],
  [
    have_cxx11_move_constructors=no

    AC_MSG_CHECKING(for C++11 move constructor support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <utility>
    class move_constructable_base
    {
    public:
      move_constructable_base() {}
      move_constructable_base(move_constructable_base && other) noexcept {}
    };
    class move_constructable : public move_constructable_base
    {
    public:
      move_constructable() {}
      move_constructable(move_constructable && other) noexcept : move_constructable_base(std::move(other)) {}
    };
    ]], [[
        move_constructable m1;
        move_constructable m2(std::move(m1));
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_MOVE_CONSTRUCTORS, 1, [Flag indicating whether compiler supports move construction])
        have_cxx11_move_constructors=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_MOVE_CONSTRUCTORS, test x$have_cxx11_move_constructors == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_RANGEFOR],
  [
    have_cxx11_rangefor=no

    AC_MSG_CHECKING(for C++11 range-based for loop support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <iostream>
        @%:@include <vector>
        void print(const std::vector<int> & v)
        {
          for (const int & x : v)
            std::cout << x << ' ';
          std::cout << std::endl;
        }
    ]], [[
        std::vector<int> v(3);
        print(v);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_RANGEFOR, 1, [Flag indicating whether compiler supports range-based for loops])
        have_cxx11_rangefor=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_RANGEFOR, test x$have_cxx11_rangefor == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_DECLTYPE],
  [
    have_cxx11_decltype=no

    AC_MSG_CHECKING(for C++11 decltype support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <vector>
    ]], [[
        int a;
        decltype(a) b;
        std::vector<int> vec(10);
        for (auto i = decltype(vec.size())(0); i < vec.size(); ++i)
          vec[i] += i;
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_DECLTYPE, 1, [Flag indicating whether compiler supports decltype])
        have_cxx11_decltype=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_DECLTYPE, test x$have_cxx11_decltype == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_RVALUE_REFERENCES],
  [
    have_cxx11_rvalue_references=no

    AC_MSG_CHECKING(for C++11 rvalue references support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      int foo(int && x) { return x; }
      int bar() { return 4; }
    ]], [[
      // Call function that takes an rvalue reference.
      foo (bar());
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_RVALUE_REFERENCES, 1, [Flag indicating whether compiler supports rvalue references])
        have_cxx11_rvalue_references=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_RVALUE_REFERENCES, test x$have_cxx11_rvalue_references == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_AUTO],
  [
    have_cxx11_auto=no

    AC_MSG_CHECKING(for C++11 auto keyword support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    ]], [[
      int x = 5;
      auto y = x;
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_AUTO, 1, [Flag indicating whether compiler supports the auto keyword])
        have_cxx11_auto=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_AUTO, test x$have_cxx11_auto == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_LAMBDA],
  [
    have_cxx11_lambda=no

    AC_MSG_CHECKING(for C++11 lambda support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      // typedef for a function pointer that takes int and returns bool.
      typedef bool (*FunctionPointer) (int);

      // A function that takes a pointer to a function that takes an int,
      // calls it with the number 4, and returns the result.
      bool f(FunctionPointer g) { return g(4); }
    ]], [[
      // Call f, passing it a lambda constructed on the fly instead
      // of a standard function pointer.  The result should be true.
      f ( [](int x) { return x > 3; } );
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_LAMBDA, 1, [Flag indicating whether compiler supports lambdas])
        have_cxx11_lambda=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_LAMBDA, test x$have_cxx11_lambda == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_CONSTEXPR],
  [
    have_cxx11_constexpr=no

    AC_MSG_CHECKING(for C++11 constexpr support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      constexpr int multiply (int x, int y) { return x * y; }
    ]], [[
      // The compiler should compute "val" at compile time.
      const int val = multiply(10, 10);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_CONSTEXPR, 1, [Flag indicating whether compiler supports constexpr])
        have_cxx11_constexpr=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_CONSTEXPR, test x$have_cxx11_constexpr == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_ALIAS_DECLARATIONS],
  [
    have_cxx11_alias_declarations=no

    AC_MSG_CHECKING(for C++11 alias declarations support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      template <typename T>
      struct check
      {
        T t;
      };

      // An alias declaration is like a templated typedef
      template <typename T>
      using MyCheck = check<T>;
    ]], [[
      MyCheck<int> mc;
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_ALIAS_DECLARATIONS, 1, [Flag indicating whether compiler supports alias declarations])
        have_cxx11_alias_declarations=yes
    ],[
        AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_ALIAS_DECLARATIONS, test x$have_cxx11_alias_declarations == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_SHARED_PTR],
  [
    have_cxx11_shared_ptr=no

    AC_MSG_CHECKING(for C++11 std::shared_ptr support)
    AC_LANG_PUSH([C++])

    # Save any original value that CXXFLAGS had
    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

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

    dnl Only set the header file variable if our flag was set to 'yes'.
    AS_IF([test "x$have_cxx11_shared_ptr" = "xyes"],
          [AC_DEFINE(HAVE_CXX11_SHARED_PTR, 1, [Flag indicating whether compiler supports std::shared_ptr])])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_SHARED_PTR, test x$have_cxx11_shared_ptr == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_UNIQUE_PTR],
  [
    have_cxx11_unique_ptr=no

    AC_MSG_CHECKING(for C++11 std::unique_ptr support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

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

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_UNIQUE_PTR, test x$have_cxx11_unique_ptr == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX14_MAKE_UNIQUE],
  [
    have_cxx14_make_unique=no

    # std::make_unique is actually part of the C++14 standard, but some
    # compilers might (?) support it via the -std=c++11 flag, or eventually
    # with no flag at all.
    AC_MSG_CHECKING(for C++14 std::make_unique support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <memory>
        ]], [[
    {
      // Normally, you would use "auto" on the LHS here to avoid
      // repeating the type name, but we are not testing auto here.
      std::unique_ptr<int> up = std::make_unique<int>(42);
    } // Foo deleted when up goes out of scope
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX14_MAKE_UNIQUE, 1, [Flag indicating whether compiler supports std::make_unique])
        have_cxx14_make_unique=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX14_MAKE_UNIQUE, test x$have_cxx14_make_unique == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_MAKE_UNIQUE_WORKAROUND],
  [
    have_cxx11_make_unique_workaround=no

    # This is a simple workaround for no std::make_unique in C++11:
    # http://stackoverflow.com/questions/7038357/make-unique-and-perfect-forwarding
    # Requires working rvalue references, std::forward, variadic
    # templates, and std::unique_ptr from C++11.
    AC_MSG_CHECKING(for C++11 std::make_unique workaround support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <memory>
    namespace local
    {
      template<typename T, typename... Args>
      std::unique_ptr<T> make_unique(Args&&... args)
      {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
      }
    }
        ]], [[
    {
      // Normally, you would use "auto" on the LHS here to avoid
      // repeating the type name, but we are not testing auto here.
      std::unique_ptr<int> up = local::make_unique<int>(42);
    } // Foo deleted when up goes out of scope
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_MAKE_UNIQUE_WORKAROUND, 1, [Flag indicating whether compiler supports C++11 std::make_unique workaround])
        have_cxx11_make_unique_workaround=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_MAKE_UNIQUE_WORKAROUND, test x$have_cxx11_make_unique_workaround == xyes)
  ])




AC_DEFUN([LIBMESH_TEST_CXX11_REGEX],
  [
    have_cxx11_regex=no

    AC_MSG_CHECKING(for C++11 std::regex support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    dnl We actually have to try and *run* the test program, since
    dnl GCC up to 4.8 will compile this but then is not able to run it.
    dnl GCC 4.9.1 and Clang 3.5 are actually able to run this test code.
    dnl
    dnl Note the quadruple backslash below -- this is needed so it
    dnl expands to 2 backslashes in the test program generated by
    dnl Autoconf...
    AC_RUN_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <regex>
    ]], [[
      std::regex integer_regex("(\\\\+|-)?[[:digit:]]+");
      std::regex_match("abc", integer_regex);
      std::regex_match("123", integer_regex);
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

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_REGEX, test x$have_cxx11_regex == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_OVERRIDE],
  [
    have_cxx11_override=no

    AC_MSG_CHECKING(for C++11 override keyword support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

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

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_OVERRIDE, test x$have_cxx11_override == xyes)
  ])



AC_DEFUN([LIBMESH_TEST_CXX11_INITIALIZER_LIST],
  [
    have_cxx11_initializer_list=no

    AC_MSG_CHECKING(for C++11 initializer list support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

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

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_INITIALIZER_LIST, test x$have_cxx11_initializer_list == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_VARIADIC_TEMPLATES],
  [
    have_cxx11_variadic_templates=no

    AC_MSG_CHECKING(for C++11 variadic template support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      // Base case
      template <typename T>
      T sum(T t) { return t; }

      // Compute sum of arbitrary number of passed parameters.
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

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_VARIADIC_TEMPLATES, test x$have_cxx11_variadic_templates == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_THREAD],
  [
    have_cxx11_thread=no

    AC_MSG_CHECKING(for C++11 <thread> support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <thread>
      @%:@include <atomic>
      @%:@include <mutex>
      void my_thread_func() {}
    ]], [[
      thread_local int i;
      std::thread t(my_thread_func);
      t.join();

      std::atomic<bool> ab1, ab2;
      ab1.store(true, std::memory_order_relaxed);
      ab2.store(false, std::memory_order_relaxed);
      ab1.exchange(ab2);

      std::mutex m;
      std::lock_guard<std::mutex> lock(m);

      std::atomic_thread_fence(std::memory_order_acquire);
      std::atomic_thread_fence(std::memory_order_release);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_THREAD, 1, [Flag indicating whether compiler supports std::thread])
        have_cxx11_thread=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_THREAD, test x$have_cxx11_thread == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_CONDITION_VARIABLE],
  [
    have_cxx11_condition_variable=no

    AC_MSG_CHECKING(for C++11 <condition_variable> support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    # Test code is from the accepted answer on:
    # http://stackoverflow.com/questions/16350473/why-do-i-need-stdcondition-variable
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <iostream>
      @%:@include <condition_variable>
      @%:@include <mutex>
      @%:@include <thread>
      @%:@include <chrono>

      bool is_ready = false;
      std::mutex m;
      std::condition_variable cv;

      void
      test()
      {
        std::this_thread::sleep_for(std::chrono::seconds(30));
        std::unique_lock<std::mutex> ulock(m);
        is_ready = true;
        cv.notify_one();
      }
    ]], [[
      std::thread t(test);
      std::unique_lock<std::mutex> ulock(m);
      while (!is_ready)
        {
          cv.wait(ulock);
          if (!is_ready)
            std::cout << "Spurious wake up!\n";
        }
      t.join();
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_CONDITION_VARIABLE, 1, [Flag indicating whether compiler supports std::condition_variable])
        have_cxx11_condition_variable=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_CONDITION_VARIABLE, test x$have_cxx11_condition_variable == xyes)
  ])



AC_DEFUN([LIBMESH_TEST_CXX11_TYPE_TRAITS],
  [
    # This is designed to be an exhaustive test of the capabilities of
    # the <type_traits> header, as defined by the C++11 standard.  We
    # don't use all of these in libmesh, but we want application code
    # to be able to rely on the value of LIBMESH_HAVE_CXX11_TYPE_TRAITS.
    #
    # Not all compilers fully implement the header.  For example, GCC
    # 4.9.1 does not provide the is_trivially_copyable() function, but
    # it does provide support for other <type_traits> functions.
    #
    # See also:
    # http://en.cppreference.com/w/cpp/header/type_traits
    have_cxx11_type_traits=no

    AC_MSG_CHECKING(for C++11 <type_traits> support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <iostream>
      @%:@include <type_traits>

      // std::enable_if - the return type of the function is only defined if
      // T is an integral type.  Therefore, it's a *compile* error if you
      // try to call it with a non-integral type.
      template <class T>
      typename std::enable_if<std::is_integral<T>::value, bool>::type
      is_odd (T i)
      {
        return static_cast<bool>(i%2);
      }

      // std::underlying_type - names the underlying type of an enum
      enum e1 {};

      // typedef (named fn_ptr) for a function that takes nothing and returns char.
      typedef char (*fn_ptr)();
    ]], [[
      std::cout << std::is_void<char>::value
        // << std::is_null_pointer<char>::value                  // C++14
                << std::is_integral<char>::value
                << std::is_floating_point<char>::value
                << std::is_array<char>::value
                << std::is_enum<char>::value
                << std::is_union<char>::value
                << std::is_class<char>::value
                << std::is_function<char>::value
                << std::is_pointer<char>::value
                << std::is_lvalue_reference<char>::value
                << std::is_rvalue_reference<char>::value
                << std::is_member_object_pointer<char>::value
                << std::is_fundamental<char>::value
                << std::is_arithmetic<char>::value
                << std::is_scalar<char>::value
                << std::is_object<char>::value
                << std::is_compound<char>::value
                << std::is_reference<char>::value
                << std::is_member_pointer<char>::value
                << std::is_const<char>::value
                << std::is_volatile<char>::value
                << std::is_trivial<char>::value
                << std::is_trivially_copyable<char>::value // Not supported by GCC 4.6.3 with -std=c++0x
                << std::is_standard_layout<char>::value
                << std::is_pod<char>::value
                << std::is_literal_type<char>::value
                << std::is_empty<char>::value
                << std::is_polymorphic<char>::value
                << std::is_abstract<char>::value
                << std::is_signed<char>::value
                << std::is_unsigned<char>::value
                << std::is_constructible<char>::value
                << std::is_trivially_constructible<char>::value
                << std::is_nothrow_constructible<char>::value
                << std::is_default_constructible<char>::value
                << std::is_trivially_default_constructible<char>::value
                << std::is_nothrow_default_constructible<char>::value
                << std::is_copy_constructible<char>::value
                << std::is_trivially_copy_constructible<char>::value
                << std::is_nothrow_copy_constructible<char>::value
                << std::is_move_constructible<char>::value
                << std::is_trivially_move_constructible<char>::value
                << std::is_nothrow_move_constructible<char>::value
                << std::is_assignable<char, char>::value
                << std::is_trivially_assignable<char, char>::value
                << std::is_nothrow_assignable<char, char>::value
                << std::is_copy_assignable<char>::value
                << std::is_trivially_copy_assignable<char>::value
                << std::is_nothrow_copy_assignable<char>::value
                << std::is_move_assignable<char>::value
                << std::is_trivially_move_assignable<char>::value
                << std::is_nothrow_move_assignable<char>::value
                << std::is_destructible<char>::value
                << std::is_trivially_destructible<char>::value
                << std::is_nothrow_destructible<char>::value
                << std::has_virtual_destructor<char>::value
                << std::alignment_of<char>::value
                << std::rank<char>::value
                << std::extent<char>::value
                << std::is_same<char, char>::value
                << std::is_base_of<char, char>::value
                << std::is_convertible<char, char>::value
                << std::is_same<char, std::remove_cv<const char>::type>::value // std::remove_cv
                << std::is_same<char, std::remove_const<const char>::type>::value // std::remove_const
                << std::is_same<char, std::remove_volatile<volatile char>::type>::value // std::remove_volatile
                << std::is_same<const volatile char, std::add_cv<char>::type>::value // std::add_cv
                << std::is_same<const char, std::add_const<char>::type>::value // std::add_const
                << std::is_same<volatile char, std::add_volatile<char>::type>::value // std::add_volatile
                << std::is_same<char, std::remove_reference<char &>::type>::value // std::remove_reference
                << std::is_same<char &, std::add_lvalue_reference<char>::type>::value // std::add_lvalue_reference
                << std::is_same<char &&, std::add_rvalue_reference<char>::type>::value // std::add_rvalue_reference
                << std::is_same<char, std::remove_pointer<char *>::type>::value // std::remove_pointer
                << std::is_same<char *, std::add_pointer<char>::type>::value // std::add_pointer
                << std::is_same<char, std::make_signed<unsigned char>::type>::value // std::make_signed
                << std::is_same<unsigned char, std::make_unsigned<char>::type>::value // std::make_unsigned
                << std::is_same<char, std::remove_extent<char>::type>::value // std::remove_extent
                << std::is_same<char, std::remove_all_extents<char>::type>::value // std::remove_all_extents
                << std::is_same<char, std::decay<const char &>::type>::value // std::decay
                << is_odd(13) // std::enable_if
                << std::is_same<char, std::conditional<true, /*type if true*/char, /*type if false*/int>::type>::value // std::conditional
                << std::is_same<long, std::common_type<char, short, int, long>::type>::value // std::common_type
                << std::is_same<int, std::underlying_type<e1>::type>::value // std::underlying_type
                << std::is_same<char, std::result_of<fn_ptr()>::type>::value // std::result_of
                << std::false_type::value
                << std::true_type::value
                << std::endl;

      // std::aligned_storage
      typedef std::aligned_storage</*store objects of length=*/1>::type aligned_t;

      // std::aligned_union
      typedef std::aligned_union</*size of at least*/32, int, char, double>::type union_t;
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_TYPE_TRAITS, 1, [Flag indicating whether compiler supports <type_traits>])
        have_cxx11_type_traits=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_TYPE_TRAITS, test x$have_cxx11_type_traits == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_MATH_FUNCS],
  [
    have_cxx11_inverse_hyperbolic_sine=no
    have_cxx11_inverse_hyperbolic_cosine=no
    have_cxx11_inverse_hyperbolic_tangent=no

    have_cxx11_inverse_hyperbolic_sine_complex=no
    have_cxx11_inverse_hyperbolic_cosine_complex=no
    have_cxx11_inverse_hyperbolic_tangent_complex=no

    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    # Test for asinh
    AC_MSG_CHECKING(for C++11 std::asinh support in <cmath>)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <cmath>
    ]], [[
      double x = std::asinh(1.);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_INVERSE_HYPERBOLIC_SINE, 1, [Flag indicating whether compiler supports std::asinh])
        have_cxx11_inverse_hyperbolic_sine=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Test for acosh
    AC_MSG_CHECKING(for C++11 std::acosh support in <cmath>)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <cmath>
    ]], [[
      double x = std::acosh(1.);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_INVERSE_HYPERBOLIC_COSINE, 1, [Flag indicating whether compiler supports std::acosh])
        have_cxx11_inverse_hyperbolic_cosine=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Test for atanh
    AC_MSG_CHECKING(for C++11 std::atanh support in <cmath>)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <cmath>
    ]], [[
      double x = std::atanh(0.);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_INVERSE_HYPERBOLIC_TANGENT, 1, [Flag indicating whether compiler supports std::atanh])
        have_cxx11_inverse_hyperbolic_tangent=yes
    ],[
      AC_MSG_RESULT(no)
    ])


    # Test for asinh(complex)
    AC_MSG_CHECKING(for C++11 std::asinh(complex) support in <complex>)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <complex>
    ]], [[
      std::complex<double> z(0, -2);
      std::complex<double> x = std::asinh(z);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_INVERSE_HYPERBOLIC_SINE_COMPLEX, 1, [Flag indicating whether compiler supports std::asinh(complex)])
        have_cxx11_inverse_hyperbolic_sine_complex=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Test for acosh(complex)
    AC_MSG_CHECKING(for C++11 std::acosh(complex) support in <complex>)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <complex>
    ]], [[
      std::complex<double> z(0.5, 0);
      std::complex<double> x = std::acosh(z);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_INVERSE_HYPERBOLIC_COSINE_COMPLEX, 1, [Flag indicating whether compiler supports std::acosh(complex)])
        have_cxx11_inverse_hyperbolic_cosine_complex=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Test for atanh(complex)
    AC_MSG_CHECKING(for C++11 std::atanh(complex) support in <complex>)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <complex>
    ]], [[
      std::complex<double> z(2, 0);
      std::complex<double> x = std::atanh(z);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_INVERSE_HYPERBOLIC_TANGENT_COMPLEX, 1, [Flag indicating whether compiler supports std::atanh(complex)])
        have_cxx11_inverse_hyperbolic_tangent_complex=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Test for erf
    AC_MSG_CHECKING(for C++11 std::erf support in <cmath>)
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <cmath>
    ]], [[
      double val = std::erf(1.);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_ERF, 1, [Flag indicating whether compiler supports std::erf])
        have_cxx11_erf=yes
    ],[
      AC_MSG_RESULT(no)
    ])


    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_INVERSE_HYPERBOLIC_SINE, test x$have_cxx11_inverse_hyperbolic_sine == xyes)
    AM_CONDITIONAL(HAVE_CXX11_INVERSE_HYPERBOLIC_COSINE, test x$have_cxx11_inverse_hyperbolic_cosine == xyes)
    AM_CONDITIONAL(HAVE_CXX11_INVERSE_HYPERBOLIC_TANGENT, test x$have_cxx11_inverse_hyperbolic_tangent == xyes)

    AM_CONDITIONAL(HAVE_CXX11_INVERSE_HYPERBOLIC_SINE_COMPLEX, test x$have_cxx11_inverse_hyperbolic_sine_complex == xyes)
    AM_CONDITIONAL(HAVE_CXX11_INVERSE_HYPERBOLIC_COSINE_COMPLEX, test x$have_cxx11_inverse_hyperbolic_cosine_complex == xyes)
    AM_CONDITIONAL(HAVE_CXX11_INVERSE_HYPERBOLIC_TANGENT_COMPLEX, test x$have_cxx11_inverse_hyperbolic_tangent_complex == xyes)

    AM_CONDITIONAL(HAVE_CXX11_ERF, test x$have_cxx11_erf == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_DELETED_FUNCTIONS],
  [
    have_cxx11_deleted_functions=no

    AC_MSG_CHECKING(for C++11 deleted functions support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    class Foo
    {
      Foo(const Foo &) = delete;
    };
    ]], [[
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_DELETED_FUNCTIONS, 1, [Flag indicating whether compiler supports f() = delete;])
        have_cxx11_deleted_functions=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_DELETED_FUNCTIONS, test x$have_cxx11_deleted_functions == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_DEFAULTED_FUNCTIONS],
  [
    have_cxx11_defaulted_functions=no

    AC_MSG_CHECKING(for C++11 defaulted functions support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    class Foo
    {
      Foo(const Foo &) = default;
      ~Foo();
    };
    Foo::~Foo() = default;
    ]], [[
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_DEFAULTED_FUNCTIONS, 1, [Flag indicating whether compiler supports defaulted functions])
        have_cxx11_defaulted_functions=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_DEFAULTED_FUNCTIONS, test x$have_cxx11_defaulted_functions == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_FINAL],
  [
    have_cxx11_final=no

    AC_MSG_CHECKING(for C++11 'final' keyword support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    // Test that a function can be declared final.
    struct A
    {
      virtual void foo() final;
    };

    // Test that a struct can be declared final.
    struct B final : A
    {
    };
    ]], [[
    ]])],[
      have_cxx11_final=yes
    ],[
      have_cxx11_final=no
    ])

    # Confirm that you cannot declare a non-virtual function 'final'.
    AS_IF([test "x$have_cxx11_final" != "xno"],
          [
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
            struct A
            {
              // Error: non-virtual function cannot be final
              void bar() final;
            };
            ]], [[
            ]])],[
              # If this code compiles, 'final' is not working correctly.
              have_cxx11_final=no
            ],[
              have_cxx11_final=yes
            ])
          ])

    dnl Confirm that you cannot override a final function.
    AS_IF([test "x$have_cxx11_final" != "xno"],
          [
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
            struct A
            {
              virtual void foo() final;
            };
            struct B : A
            {
              // Error: foo cannot be overridden as it's final in A
              void foo();
            };
            ]], [[
            ]])],[
              # If this code compiles, 'final' is not working correctly.
              have_cxx11_final=no
            ],[
              have_cxx11_final=yes
            ])
          ])

    dnl Confirm that you cannot inherit from a 'final' class.
    AS_IF([test "x$have_cxx11_final" != "xno"],
          [
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
            struct A
            {
            };

            // struct B is final
            struct B final : A
            {
            };

            // Error: B is final
            struct C : B
            {
            };
            ]], [[
            ]])],[
              # If this code compiles, 'final' is not working correctly.
              have_cxx11_final=no
            ],[
              have_cxx11_final=yes
            ])
          ])

    # If the flag is still 'yes' after all the tests, set the #define.
    AS_IF([test "x$have_cxx11_final" = "xyes"],
          [
            AC_MSG_RESULT(yes)
            AC_DEFINE(HAVE_CXX11_FINAL, 1, [Flag indicating whether compiler supports f() final;])
          ],
          [AC_MSG_RESULT(no)])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_FINAL, test x$have_cxx11_final == xyes)
  ])


AC_DEFUN([LIBMESH_TEST_CXX11_NULLPTR],
  [
    have_cxx11_nullptr=no

    AC_MSG_CHECKING(for C++11 nullptr support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    @%:@include <cstddef>
    void f(int * pi) {}
    void f(double * pd) {}
    void f(std::nullptr_t nullp) {}
    ]], [[
    // would be ambiguous without void f(nullptr_t)
    f(nullptr);
    ]])],[
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_CXX11_NULLPTR, 1, [Flag indicating whether compiler supports nullptr])
      have_cxx11_nullptr=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_NULLPTR, test x$have_cxx11_nullptr == xyes)
  ])



AC_DEFUN([LIBMESH_TEST_CXX11_TO_STRING],
  [
    have_cxx11_to_string=no

    AC_MSG_CHECKING(for C++11 std::to_string() support)
    AC_LANG_PUSH([C++])

    old_CXXFLAGS="$CXXFLAGS"
    CXXFLAGS="$CXXFLAGS $switch $libmesh_CXXFLAGS"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      @%:@include <string>
    ]], [[
      // tiny="0.000000".  Note: std::to_string(double) is required to produce
      // a std::string with the same contents as std::sprintf(buf, "%f", value)
      // would produce, given a sufficiently large buf.  This is *different* from
      // what you get from a std::stringstream using default formatting and
      // precision flags, i.e.
      //   std::ostringstream oss;
      //   oss << 1.e-40;
      //   std::string tiny = oss.str();
      // will produce the string "1e-40".
      std::string tiny = std::to_string(1.e-40);
    ]])],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_CXX11_TO_STRING, 1, [Flag indicating whether compiler supports std::to_string()])
        have_cxx11_to_string=yes
    ],[
      AC_MSG_RESULT(no)
    ])

    # Reset the flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])

    AM_CONDITIONAL(HAVE_CXX11_TO_STRING, test x$have_cxx11_to_string == xyes)
  ])
