# ----------------------------------------------------------------
# Locate header files for the C++ linear algebra library Eigen.
# Eigen is a header-only template library. By default we check for the
# Eigen files in the --with-eigen-include=xxx argument provided to
# configure, or if those don't exist in the $EIGEN_INC/Eigen directory,
# or in /usr/include.
#
# If an argument to this macro is specified, it should give the minimum useful
# version number.
#
# If a second argument to this macro is specified as 'yes', then this is a
# required dependency and configure will exit with an error if it is not
# satisfied.
#
# Note: Eigen is installed (by default) at the location
# /path/to/eigen/Eigen, i.e. with path ending in capital 'Eigen'.
# You should specify --with-eigen-include=/path/to/eigen
# during configure, or set your $EIGEN_INC environment variable
# to /path/to/eigen.
# ----------------------------------------------------------------

AC_DEFUN([CONFIGURE_EIGEN],
[
  AC_ARG_ENABLE(eigen,
                AS_HELP_STRING([--disable-eigen],
                               [build without Eigen linear algebra support]),
		[case "${enableval}" in
		  yes)  enableeigen=yes ;;
		   no)  enableeigen=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-eigen) ;;
		 esac],
		 [enableeigen=$enableoptional])

  # package requirement; if not specified, the default is to assume that
  # the package is optional

  is_package_required=ifelse([$2], ,no, $2 )

  install_internal_eigen=no
  if (test x$enableeigen = xyes); then

    # User-specific include path
    AC_ARG_WITH(eigen-include,
                AS_HELP_STRING([--with-eigen-include=PATH],[Specify the path for EIGEN header files]),
                witheigeninc=$withval,
                witheigeninc=no)

    # Fall back on default paths to Eigen's include files
    if (test $witheigeninc != no); then
	EIGEN_INC="$witheigeninc"

    elif (test "x$EIGEN_INC" != x -a -f $EIGEN_INC/Eigen/Eigen); then
	echo "Environment EIGEN_INC=$EIGEN_INC"

    elif (test "x$EIGEN3_INCLUDE" != x -a -f $EIGEN3_INCLUDE/Eigen/Eigen); then
	EIGEN_INC=$EIGEN3_INCLUDE
	echo "Environment EIGEN_INC=$EIGEN_INC"

    elif (test "x$EIGEN_INCLUDE" != x -a -f $EIGEN_INCLUDE/Eigen/Eigen); then
	EIGEN_INC=$EIGEN_INCLUDE
	echo "Environment EIGEN_INC=$EIGEN_INC"

    elif (test -f /usr/include/eigen3/Eigen/Eigen); then
	EIGEN_INC="/usr/include/eigen3"

    else
	EIGEN_INC="/usr/include"
    fi

    # Initialize Makefile/config.h substitution variables
    EIGEN_INCLUDE=""

    # Check for existence of a header file in the specified location.  Note: here
    # we are checking for the header file "Eigen" in the Eigen directory.
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS

    externaleigenincFound=no;
    AC_CHECK_HEADERS($EIGEN_INC/Eigen/Eigen, externaleigenincFound=yes)

    # Check to make sure the external header files are sufficiently up
    # to date - this fixes our Eigen detection on Scientific Linux 6
    if (test x$externaleigenincFound = xyes); then
        ac_eigen_save_CPPFLAGS="$CPPFLAGS"
	CPPFLAGS="-I${EIGEN_INC} ${CPPFLAGS}"

	AC_CHECK_HEADERS([Eigen/Dense],[],[enableeigenincFound=no])

        if (test x$enableeigensparse = xyes); then
	    AC_CHECK_HEADERS([Eigen/Sparse],[],[enableeigenincFound=no])
        fi

	CPPFLAGS="${ac_eigen_save_CPPFLAGS}"
    fi


    if (test x$externaleigenincFound = xyes); then
        EIGEN_INCLUDE="-I$EIGEN_INC"
    elif (test -d $top_srcdir/contrib/eigen/eigen); then
        AC_MSG_RESULT([<<< external Eigen header files not found, using Eigen from ./contrib >>>])
	EIGEN_INC=$top_srcdir/contrib/eigen/eigen
	EIGEN_INCLUDE="-I\$(top_srcdir)/contrib/eigen/eigen"
	install_internal_eigen=yes
    else
	enableeigen=no
    fi


    # OK, we have a usable eigen path, make sure the headers we want are good.
    if (test x$enableeigen = xyes); then

        ac_eigen_save_CPPFLAGS="$CPPFLAGS"
	CPPFLAGS="-I${EIGEN_INC} ${CPPFLAGS}"

	AC_CHECK_HEADERS([Eigen/Dense],[],[enableeigen=no])

        if (test x$enableeigensparse = xyes); then
	    AC_CHECK_HEADERS([Eigen/Sparse],[],[enableeigen=no])
        fi

        #-----------------------
        # Minimum version check
        #----------------------

        min_eigen_version=ifelse([$1], ,0.0.0, $1)

        # looking for major.minor.micro (which Eigen calls world.major.minor) style versioning

        MAJOR_VER=`echo $min_eigen_version | sed 's/^\([[0-9]]*\).*/\1/'`
        if test "x${MAJOR_VER}" = "x" ; then
           MAJOR_VER=0
        fi

        MINOR_VER=`echo $min_eigen_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
        if test "x${MINOR_VER}" = "x" ; then
           MINOR_VER=0
        fi

        MICRO_VER=`echo $min_eigen_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
        if test "x${MICRO_VER}" = "x" ; then
           MICRO_VER=0
        fi


        AC_MSG_CHECKING(for eigen - version >= $min_eigen_version)

        AC_LANG_PUSH([C++])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include "Eigen/Core"
            ]], [[
            #if EIGEN_WORLD_VERSION > $MAJOR_VER
            #elif (EIGEN_WORLD_VERSION >= $MAJOR_VER) && (EIGEN_MAJOR_VERSION > $MINOR_VER)
            #elif (EIGEN_WORLD_VERSION >= $MAJOR_VER) && (EIGEN_MAJOR_VERSION >= $MINOR_VER) && (EIGEN_MINOR_VERSION >= $MICRO_VER)
            #else
            #  error version is too old
            #endif
        ]])],[
            AC_MSG_RESULT(yes)
        ],[
            AC_MSG_RESULT(no)
            enableeigen=no
        ])
        AC_LANG_POP([C++])

	CPPFLAGS="${ac_eigen_save_CPPFLAGS}"

	# if we survived, we really have Eigen
 	if (test x$enableeigen = xyes); then
            HAVE_EIGEN=1
	    AC_DEFINE(HAVE_EIGEN, 1, [Flag indicating whether the library will be compiled with Eigen support])
	    AC_MSG_RESULT(<<< Configuring library with Eigen support >>>)
        elif test "$is_package_required" = yes; then
            AC_MSG_ERROR([

   Your EIGEN version ($EIGEN_INC) does not meet the minimum versioning
   requirements ($min_eigen_version).  Please use --with-eigen-include to
   specify the location of an updated installation.

                ])
	fi
    fi

    AC_LANG_RESTORE
  fi

  # Substitute the substitution variables
  AC_SUBST(EIGEN_INCLUDE)
])
