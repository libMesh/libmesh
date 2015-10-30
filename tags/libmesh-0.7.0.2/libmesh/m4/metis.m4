dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl Metis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_METIS], 
[
  AC_CHECK_FILE(./contrib/metis/Lib/metis.h,
	        [ 
	          METIS_INCLUDE_PATH=$PWD/contrib/metis/Lib
                  METIS_INCLUDE=-I$METIS_INCLUDE_PATH
                  METIS_LIB="\$(EXTERNAL_LIBDIR)/libmetis\$(libext)"
		  AC_SUBST(METIS_INCLUDE)
                  AC_SUBST(METIS_LIB)
                  AC_DEFINE(HAVE_METIS, 1,
	                     [Flag indicating whether or not Metis is available])
                  AC_MSG_RESULT(<<< Configuring library with Metis support >>>)
	          enablemetis=yes
                ],
                [enablemetis=no])
])
