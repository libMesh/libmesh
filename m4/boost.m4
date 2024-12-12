# --------------------------------------------------------------
# Look for a user or system-provided boost.  If there is not
# one available then use the minimal ./contrib/boost provided.
# --------------------------------------------------------------
AC_DEFUN([CONFIGURE_BOOST],
[
  AC_ARG_ENABLE(boost,
                AS_HELP_STRING([--disable-boost],
                               [build without external BOOST support]),
                [AS_CASE("${enableval}",
                         [yes], [enableboost=yes],
                         [no],  [enableboost=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-boost)])],
                enableboost=$enableoptional)

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
                [AC_DEFINE(HAVE_EXTERNAL_BOOST, 1, [Flag indicating whether an external Boost was found])],
                [enableboost=no])
        ])


  dnl if we are installing our own Boost, add it to the contrib search path
  dnl which is not exported during install.  If we are using an external boost,
  dnl add its (absolute) path, as determined by AX_BOOST_BASE to the
  dnl libmesh_optional_INCLUDES variable.
  dnl In either case, we also add Boost's absolute path to timpi_CPPFLAGS so
  dnl that, if we are using quadruple precision, TIMPI can use Boost's headers.
  AS_IF([test "x$enableboost" = "xyes"],
        [
         libmesh_optional_INCLUDES="$BOOST_CPPFLAGS $libmesh_optional_INCLUDES"
         export timpi_CPPFLAGS="$timpi_CPPFLAGS $BOOST_CPPFLAGS"
        ])
])
