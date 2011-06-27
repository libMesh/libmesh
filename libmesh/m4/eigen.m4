dnl -------------------------------------------------------------
dnl $Id: eigen.m4 4255 2011-03-12 21:13:18Z benkirk $
dnl -------------------------------------------------------------

dnl ----------------------------------------------------------------
dnl Locate header files for the C++ linear algebra library Eigen.
dnl Eigen is a header-only template library. By default we check for the
dnl Eigen files in the --with-eigen-include=xxx argument provided to
dnl configure, or if those don't exist in the $EIGEN_INC directory,
dnl or in /usr/include.
dnl ----------------------------------------------------------------

AC_DEFUN([CONFIGURE_EIGEN], 
[
  dnl User-specific include path
  AC_ARG_WITH(eigen-include,
              AC_HELP_STRING([--with-eigen-include=PATH],[Specify the path for EIGEN header files]),
              witheigeninc=$withval,
              witheigeninc=no)

  dnl Fall back on default paths to Eigen's include files
  if (test $witheigeninc != no); then
    EIGEN_INC="$witheigeninc"
  elif test "x$EIGEN_INC" != x -a -f $EIGEN_INC/Eigen; then
    echo "Environment EIGEN_INC=$EIGEN_INC"
  else
    EIGEN_INC="/usr/include"
  fi

  dnl Initialize Makefile/config.h substitution variables
  EIGEN_INCLUDE=""

  dnl Properly let the substitution variables
  if (test $enableeigen = yes); then
  
     dnl Check for existence of a header file in the specified location
     dnl AC_CHECK_FILE([$EIGEN_INC/Eigen], [eigenincFound="OK"], [eigenincFound="FAIL"])
     eigenincFound=no;
     AC_CHECK_HEADERS($EIGEN_INC/Eigen, eigenincFound=yes)

     if (test $eigenincFound = no); then
       AC_MSG_RESULT(Eigen header files not found!)
       enableeigen=no;
     fi

     dnl Print out a helpful message if we successfully located the headers
     if (test $enableeigen = yes); then
       AC_MSG_RESULT(<<< Configuring library with Eigen support >>>)
     fi
     
     dnl If the Eigen headers were found, continue.
     if (test $enableeigen = yes); then
       EIGEN_INCLUDE="-I$EIGEN_INC"
       AC_DEFINE(HAVE_EIGEN, 1, [Flag indicating whether the library will be compiled with Eigen support])
       AC_MSG_RESULT(<<< Configuring library with Eigen support >>>)
     fi
  fi

  dnl Substitute the substitution variables
  AC_SUBST(EIGEN_INCLUDE)	
  AC_SUBST(enableeigen)
])
