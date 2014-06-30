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
  if (test "$enableboost" !=  no) ; then
    # --------------------------------------------------------------
    # Look for a user or system-provided boost.  If there is not
    # one available then use the minimal ./contrib/boost provided.
    # --------------------------------------------------------------
    external_boost_found=yes
    AX_BOOST_BASE([1.20.0], [AC_MSG_RESULT(<<< Using external boost installation >>>)], [external_boost_found=no], [])

    # If that did not work, try using our builtin boost.
    if test "$external_boost_found" = "no" ; then
      AC_MSG_RESULT(<<< External boost installation *not* found... will try to configure for libmesh's internal boost >>>)

      # Note: 4th argument is libmesh's builtin boost.
      internal_boost_found=yes
      AX_BOOST_BASE([1.20.0], [AC_MSG_RESULT(<<< Using libmesh-provided boost in ./contrib >>>)], [internal_boost_found=no], [$srcdir/contrib/boost])

      if test "$internal_boost_found" = "no" ; then
        AC_MSG_RESULT(<<< Libmesh boost installation *not* found >>>)
	enableboost=no
      else
        install_internal_boost=yes
      fi
    fi

    # currently only using boost headers
    libmesh_optional_INCLUDES="$BOOST_CPPFLAGS $libmesh_optional_INCLUDES"
    #libmesh_optional_LIBS="$BOOST_LDFLAGS $libmesh_optional_LIBS"
  fi
  AM_CONDITIONAL(LIBMESH_INSTALL_INTERNAL_BOOST, test x$install_internal_boost = xyes)
  # --------------------------------------------------------------
])
