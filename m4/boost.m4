# --------------------------------------------------------------
# Look for a user or system-provided boost.  If there is not
# one available then use the minimal ./contrib/boost provided.
# --------------------------------------------------------------
AC_DEFUN([CONFIGURE_BOOST],
[
  AC_ARG_ENABLE(boost,
                AS_HELP_STRING([--disable-boost],
                               [build without either external or built-in BOOST support]),
                [case "${enableval}" in
                  yes)  enableboost=yes ;;
                  no)  enableboost=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-boost) ;;
                esac],
                enableboost=$enableoptional)

  install_internal_boost=no
  BOOST_INCLUDE=""
  if (test "$enableboost" !=  no) ; then
    # --------------------------------------------------------------
    # Look for a user or system-provided boost.  If there is not
    # one available then use the minimal ./contrib/boost provided.
    # --------------------------------------------------------------
    external_boost_found=yes
    AX_BOOST_BASE([1.57.0],
                  [AC_MSG_RESULT(<<< Using external boost installation >>>)],
                  [external_boost_found=no],
                  [])

    # Set an extra define if an external boost installation (as opposed to
    # libmesh's internal boost subset) was found. This will allow users to
    # determine at compile time when they are safe to use boost features
    # that are not available in the libmesh subset. Note: we currently don't
    # distinguish which features are available in the external boost (we assume
    # it is a full installation) but that could be added later.
    if test "$external_boost_found" = "yes" ; then
      AC_DEFINE(HAVE_EXTERNAL_BOOST, 1, [Flag indicating whether an external Boost was found])
    fi

    # If that did not work, try using our builtin boost.
    if test "$external_boost_found" = "no" ; then
      AC_MSG_RESULT(<<< External boost installation *not* found... will try to configure for libmesh's internal boost >>>)

      # Note: 4th argument is libmesh's builtin boost.
      internal_boost_found=yes
      AX_BOOST_BASE([1.57.0],
                    [AC_MSG_RESULT(<<< Using libmesh-provided boost in ./contrib >>>)],
                    [internal_boost_found=no],
                    [$top_srcdir/contrib/boost])

      if test "$internal_boost_found" = "no" ; then
        AC_MSG_RESULT(<<< Libmesh boost installation *not* found >>>)
        enableboost=no
      else
        install_internal_boost=yes
        BOOST_INCLUDE="-I\$(top_srcdir)/contrib/boost/include"
      fi
    fi
  fi


  # if we are installing our own Boost, add it to the contrib search path
  # which is not exported during install.  If we are using an external boost,
  # add its (absolute) path, as determined by AX_BOOST_BASE to the
  # libmesh_optional_INCLUDES variable.
  if (test x$enableboost = xyes); then
    if (test x$install_internal_boost = xyes); then
      libmesh_contrib_INCLUDES="$BOOST_INCLUDE $libmesh_contrib_INCLUDES"
    else
      libmesh_optional_INCLUDES="$BOOST_CPPFLAGS $libmesh_optional_INCLUDES"
    fi
  fi

  AM_CONDITIONAL(LIBMESH_INSTALL_INTERNAL_BOOST, test x$install_internal_boost = xyes)
])

AC_DEFUN([LIBMESH_TEST_BOOST_MOVELIB_UNIQUE_PTR],
[
  have_boost_unique_ptr=no
  if (test x$enableboost = xyes); then
    AC_MSG_CHECKING(for boost::movelib::unique_ptr support)
    AC_LANG_PUSH([C++])

    # Store the old value of CXXFLAGS, then append BOOST_CPPFLAGS which was
    # determined by AX_BOOST_BASE.
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

    # Reset old flags
    CXXFLAGS="$old_CXXFLAGS"
    AC_LANG_POP([C++])
  fi
])
