dnl -------------------------------------------------------------
dnl Parmetis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PARMETIS], 
[
  AC_REQUIRE([ACX_MPI])
	
  dnl We require a valid MPI installation for Parmetis
  if (test "x$MPI_IMPL" != x) ; then

    dnl need Metis for Parmetis
    AC_REQUIRE([CONFIGURE_METIS])

    if (test $enablemetis = yes) ; then
      AC_CHECK_FILE(./contrib/parmetis/Lib/parmetis.h,
      	        [
      		  
      	          PARMETIS_INCLUDE_PATH=$PWD/contrib/parmetis/Lib
                      PARMETIS_INCLUDE=-I$PARMETIS_INCLUDE_PATH
                      PARMETIS_LIB="\$(EXTERNAL_LIBDIR)/libparmetis\$(libext)"
      		  AC_SUBST(PARMETIS_INCLUDE)
                      AC_SUBST(PARMETIS_LIB)
                      AC_DEFINE(HAVE_PARMETIS, 1,
      	                     [Flag indicating whether or not ParMetis is available])
                      AC_MSG_RESULT(<<< Configuring library with ParMetis support >>>)
      	          enableparmetis=yes
                    ],
                    [enableparmetis=no])
    else
      enableparmetis=no
    fi  
  else
    enableparmetis=no
  fi
])
