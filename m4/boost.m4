# --------------------------------------------------------------
# Look for a user or system-provided boost.  If there is not
# one available then use the minimal ./contrib/boost provided.
# --------------------------------------------------------------
AC_DEFUN([CONFIGURE_BOOST],
[
  AC_ARG_ENABLE(boost,
                AS_HELP_STRING([--disable-boost],
                               [build without either external or built-in BOOST support]),
                [AS_CASE("${enableval}",
                         [yes], [enableboost=yes],
                         [no],  [enableboost=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-boost)])],
                enableboost=$enableoptional)

  install_internal_boost=no
  BOOST_INCLUDE=""
  AS_IF([test "x$enableboost" != "xno"],
        [
          dnl --------------------------------------------------------------
          dnl Look for a user or system-provided boost.  If there is not
          dnl one available then use the minimal ./contrib/boost provided.
          dnl --------------------------------------------------------------
          external_boost_found=yes
          AX_BOOST_BASE([1.57.0],
                        [AC_MSG_RESULT(<<< Using external boost installation >>>)],
                        [external_boost_found=no],
                        [])

          dnl Set an extra define if an external boost installation (as opposed to
          dnl libmesh's internal boost subset) was found. This will allow users to
          dnl determine at compile time when they are safe to use boost features
          dnl that are not available in the libmesh subset. Note: we currently don't
          dnl distinguish which features are available in the external boost (we assume
          dnl it is a full installation) but that could be added later.
          AS_IF([test "x$external_boost_found" = "xyes"],
                [AC_DEFINE(HAVE_EXTERNAL_BOOST, 1, [Flag indicating whether an external Boost was found])])

          dnl If that did not work, try using our builtin boost.
          AS_IF([test "x$external_boost_found" = "xno"],
                [
                  AC_MSG_RESULT(<<< External boost installation *not* found... will try to configure for libmesh's internal boost >>>)

                  dnl Note: 4th argument is libmesh's builtin boost.
                  internal_boost_found=yes
                  AX_BOOST_BASE([1.57.0],
                                [AC_MSG_RESULT(<<< Using libmesh-provided boost in ./contrib >>>)],
                                [internal_boost_found=no],
                                [$top_srcdir/contrib/boost])

                  AS_IF([test "x$internal_boost_found" = "xno"],
                        [
                          AC_MSG_RESULT(<<< Libmesh boost installation *not* found >>>)
                          enableboost=no
                        ],
                        [
                          install_internal_boost=yes
                          BOOST_INCLUDE="-I\$(top_srcdir)/contrib/boost/include"
                        ])
                ])
        ])


  dnl if we are installing our own Boost, add it to the contrib search path
  dnl which is not exported during install.  If we are using an external boost,
  dnl add its (absolute) path, as determined by AX_BOOST_BASE to the
  dnl libmesh_optional_INCLUDES variable.
  AS_IF([test "x$enableboost" = "xyes"],
        [
          AS_IF([test "x$install_internal_boost" = "xyes"],
                [libmesh_contrib_INCLUDES="$BOOST_INCLUDE $libmesh_contrib_INCLUDES"],
                [libmesh_optional_INCLUDES="$BOOST_CPPFLAGS $libmesh_optional_INCLUDES"])
        ])

  AM_CONDITIONAL(LIBMESH_INSTALL_INTERNAL_BOOST, test x$install_internal_boost = xyes)
])

AC_DEFUN([LIBMESH_TEST_BOOST_MOVELIB_UNIQUE_PTR],
[
  have_boost_unique_ptr=no
  AS_IF([test "x$enableboost" = "xyes"],
        [
          AC_MSG_CHECKING(for boost::movelib::unique_ptr support)
          AC_LANG_PUSH([C++])

          dnl Store the old value of CXXFLAGS, then append BOOST_CPPFLAGS which was
          dnl determined by AX_BOOST_BASE.
          old_CXXFLAGS="$CXXFLAGS"
          CXXFLAGS="$CXXFLAGS $BOOST_CPPFLAGS"

          AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
          @%:@include "boost/move/unique_ptr.hpp"
          ]], [[
              boost::movelib::unique_ptr<double> x(new double);
          ]])],[
              have_boost_unique_ptr=yes
              AC_MSG_RESULT(yes)
              AC_DEFINE(HAVE_BOOST_MOVELIB_UNIQUE_PTR, 1, [Flag indicating whether the library will use Boost Move's unique_ptr implementation])
          ],[
              AC_MSG_RESULT(no)
          ])

          dnl Reset old flags
          CXXFLAGS="$old_CXXFLAGS"
          AC_LANG_POP([C++])
        ])
])
