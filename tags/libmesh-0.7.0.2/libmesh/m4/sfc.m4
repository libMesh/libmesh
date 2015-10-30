dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl Space Filling Curves
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_SFC], 
[
  dnl Initialize variables
  SFC_INCLUDE=""
  SFC_LIB=""

  dnl Sanity check: make sure the user really has the contrib directory
  if (test $enablesfc = yes); then
    AC_CHECK_FILE(./contrib/sfcurves/sfcurves.h, [enablesfc=yes], [enablesfc=no])
  fi


  if (test $enablesfc = yes); then
     SFC_INCLUDE="-I$PWD/contrib/sfcurves"
     SFC_LIB="\$(EXTERNAL_LIBDIR)/libsfcurves\$(libext)"
     AC_DEFINE(HAVE_SFCURVES, 1, [Flag indicating whether or not Space filling curves are available])
     AC_MSG_RESULT(<<< Configuring library with SFC support >>>)
  fi

  AC_SUBST(SFC_INCLUDE)
  AC_SUBST(SFC_LIB)	
  AC_SUBST(enablesfc)

])
