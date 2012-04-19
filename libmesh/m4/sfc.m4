dnl -------------------------------------------------------------
dnl Space Filling Curves
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_SFC], 
[
  AC_ARG_ENABLE(sfc,
                AC_HELP_STRING([--enable-sfc],
                               [build with space-filling curves suppport]),
		[case "${enableval}" in
		  yes)  enablesfc=yes ;;
		   no)  enablesfc=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-sfc) ;;
		 esac],
		 [enablesfc=$enableoptional])



  dnl The SFC API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablesfc = yes); then
     SFC_INCLUDE="-I\$(top_srcdir)/contrib/sfc"
     #SFC_LIBRARY="\$(EXTERNAL_LIBDIR)/libsfc\$(libext)"
     AC_DEFINE(HAVE_SFCURVES, 1, [Flag indicating whether the library will be compiled with SFC support])
     AC_MSG_RESULT(<<< Configuring library with SFC support >>>)
     libmesh_contrib_INCLUDES="$SFC_INCLUDE $libmesh_contrib_INCLUDES"
  else
     SFC_INCLUDE=""
     #SFC_LIBRARY=""
     enablesfc=no
  fi

  AC_SUBST(SFC_INCLUDE)
  #AC_SUBST(SFC_LIBRARY)	
  AC_SUBST(enablesfc)

  AM_CONDITIONAL(LIBMESH_ENABLE_SFC, test x$enablesfc = xyes)		 
])


# [
#   dnl Initialize variables
#   SFC_INCLUDE=""
#   SFC_LIB=""

#   dnl Sanity check: make sure the user really has the contrib directory
#   if (test $enablesfc = yes); then
#     AC_CHECK_FILE(./contrib/sfcurves/sfcurves.h, [enablesfc=yes], [enablesfc=no])
#   fi


#   if (test $enablesfc = yes); then
#      SFC_INCLUDE="-I\$(top_srcdir)/contrib/sfcurves"
#      SFC_LIB="\$(EXTERNAL_LIBDIR)/libsfcurves\$(libext)"
#      AC_DEFINE(HAVE_SFCURVES, 1, [Flag indicating whether or not Space filling curves are available])
#      AC_MSG_RESULT(<<< Configuring library with SFC support >>>)
#   fi

#   AC_SUBST(SFC_INCLUDE)
#   AC_SUBST(SFC_LIB)	
#   AC_SUBST(enablesfc)

# ])
