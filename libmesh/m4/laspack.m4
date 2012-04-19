dnl -------------------------------------------------------------
dnl LASPACK Iterative Solvers
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_LASPACK], 
[
  AC_ARG_ENABLE(laspack,
                AC_HELP_STRING([--enable-laspack],
                               [build with LASPACK iterative solver suppport]),
		[case "${enableval}" in
		  yes)  enablelaspack=yes ;;
		   no)  enablelaspack=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-laspack) ;;
		 esac],
		 [enablelaspack=$enableoptional])



  dnl The LASPACK API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablelaspack = yes); then
     LASPACK_INCLUDE="-I\$(top_srcdir)/contrib/laspack"
     #LASPACK_LIBRARY="\$(EXTERNAL_LIBDIR)/liblaspack\$(libext)"
     AC_DEFINE(HAVE_LASPACK, 1, [Flag indicating whether the library will be compiled with LASPACK support])
     AC_MSG_RESULT(<<< Configuring library with Laspack support >>>)
     libmesh_contrib_INCLUDES="$LASPACK_INCLUDE $libmesh_contrib_INCLUDES"
  else
     LASPACK_INCLUDE=""
     #LASPACK_LIBRARY=""
     enablelaspack=no
  fi

  AC_SUBST(LASPACK_INCLUDE)
  #AC_SUBST(LASPACK_LIBRARY)	
  AC_SUBST(enablelaspack)

  AM_CONDITIONAL(LIBMESH_ENABLE_LASPACK, test x$enablelaspack = xyes)		 
])


# AC_DEFUN([CONFIGURE_LASPACK], 
# [

# dnl Initialize variables
# LASPACK_INCLUDE=""
# LASPACK_LIB=""

# dnl Sanity check: make sure the user really has the contrib directory
# if (test $enablelaspack = yes); then
#   AC_CHECK_FILE(./contrib/laspack/lastypes.h, [enablelaspack=yes], [enablelaspack=no])
# fi


# if (test $enablelaspack = yes); then

#   LASPACK_INCLUDE="-I\$(top_srcdir)/contrib/laspack"
#   LASPACK_LIB="\$(EXTERNAL_LIBDIR)/liblaspack\$(libext)"
#   AC_DEFINE(HAVE_LASPACK, 1, [Flag indicating whether or not LASPACK iterative solvers are available])
#   laspack_version=`grep "define LASPACK_VERSION " $PWD/contrib/laspack/version.h | sed -e "s/[[^0-9.]]*//g"`
#   AC_MSG_RESULT(<<< Configuring library with LASPACK version $laspack_version support >>>)

# fi

# AC_SUBST(LASPACK_INCLUDE)
# AC_SUBST(LASPACK_LIB)	
# AC_SUBST(enablelaspack)

# ])
