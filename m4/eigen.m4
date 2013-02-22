# ----------------------------------------------------------------
# Locate header files for the C++ linear algebra library Eigen.
# Eigen is a header-only template library. By default we check for the
# Eigen files in the --with-eigen-include=xxx argument provided to
# configure, or if those don't exist in the $EIGEN_INC/Eigen directory,
# or in /usr/include.  
dnl
# Note: Eigen is installed (by default) at the location
# /path/to/eigen/Eigen, i.e. with path ending in capital 'Eigen'.
# You should specify --with-eigen-include=/path/to/eigen
# during configure, or set your $EIGEN_INC environment variable
# to /path/to/eigen.
# ----------------------------------------------------------------

AC_DEFUN([CONFIGURE_EIGEN], 
[
  AC_ARG_ENABLE(eigen,
                AC_HELP_STRING([--enable-eigen],
                               [build with Eigen linear algebra support]),
		[case "${enableval}" in
		  yes)  enableeigen=yes ;;
		   no)  enableeigen=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-eigen) ;;
		 esac],
		 [enableeigen=$enableoptional])


  install_internal_eigen=no
  if (test $enableeigen = yes); then
  
    # User-specific include path
    AC_ARG_WITH(eigen-include,
                AC_HELP_STRING([--with-eigen-include=PATH],[Specify the path for EIGEN header files]),
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
	
    elif (test -f /usr/include/eigen3/Eigen/Eigen); then
	EIGEN_INC="/usr/include/eigen3"
	
    else
	EIGEN_INC="/usr/include"
    fi
  
    # Initialize Makefile/config.h substitution variables
    EIGEN_INCLUDE=""
  
    # Check for existence of a header file in the specified location.  Note: here
    # we are checking for the header file "Eigen" in the Eigen directory.
    externaleigenincFound=no;
    AC_CHECK_HEADERS($EIGEN_INC/Eigen/Eigen, externaleigenincFound=yes)
  
    if (test x$externaleigenincFound = xyes); then
        EIGEN_INCLUDE="-I$EIGEN_INC"
    elif (test -d $top_srcdir/contrib/eigen/eigen); then
        AC_MSG_RESULT([<<< external Eigen header files not found, using libmesh-provided Eigen in ./contrib >>>])
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

	 AC_CHECK_HEADERS([Eigen/Dense Eigen/Sparse],[],[enableeigen=no])
			  
	 CPPFLAGS="${ac_eigen_save_CPPFLAGS}"

	 # if we survived, we really have Eigen
 	 if (test x$enableeigen = xyes); then				    
	   AC_DEFINE(HAVE_EIGEN, 1, [Flag indicating whether the library will be compiled with Eigen support])
	   AC_MSG_RESULT(<<< Configuring library with Eigen support >>>)
	 fi
    fi
  fi

  # Substitute the substitution variables
  AC_SUBST(EIGEN_INCLUDE)
])
