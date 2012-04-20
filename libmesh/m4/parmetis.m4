dnl -------------------------------------------------------------
dnl Parmetis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PARMETIS], 
[
  AC_ARG_ENABLE(parmetis,
                AC_HELP_STRING([--enable-parmetis],
                               [build with Parmetis parallel graph partitioning suppport]),
		[case "${enableval}" in
		  yes)  enableparmetis=yes ;;
		   no)  enableparmetis=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-parmetis) ;;
		 esac],
		 [enableparmetis=$enableoptional])

  dnl Trump --enable-parmetis with --disable-mpi
  if (test "x$enablempi" = xno); then
    enableparmetis=no
  fi	


  dnl The PARMETIS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enableparmetis = yes); then
     PARMETIS_INCLUDE="-I\$(top_srcdir)/contrib/parmetis/Lib"
     PARMETIS_LIB="\$(EXTERNAL_LIBDIR)/libparmetis\$(libext)"
     AC_DEFINE(HAVE_PARMETIS, 1, [Flag indicating whether the library will be compiled with Parmetis support])
     AC_MSG_RESULT(<<< Configuring library with Parmetis support >>>)
     libmesh_contrib_INCLUDES="$PARMETIS_INCLUDE $libmesh_contrib_INCLUDES"
  else
     PARMETIS_INCLUDE=""
     PARMETIS_LIB=""
     enableparmetis=no
  fi

  AC_SUBST(PARMETIS_INCLUDE)
  AC_SUBST(PARMETIS_LIB)	
  AC_SUBST(enableparmetis)

  AM_CONDITIONAL(LIBMESH_ENABLE_PARMETIS, test x$enableparmetis = xyes)		 
])


# AC_DEFUN([CONFIGURE_PARMETIS], 
# [
#   AC_REQUIRE([ACX_MPI])
	
#   dnl We require a valid MPI installation for Parmetis
#   if (test "x$MPI_IMPL" != x) ; then

#     dnl need Metis for Parmetis
#     AC_REQUIRE([CONFIGURE_METIS])

#     if (test $enablemetis = yes) ; then
#       AC_CHECK_FILE(./contrib/parmetis/Lib/parmetis.h,
#       	        [
      		  
#       	          PARMETIS_INCLUDE_PATH=$PWD/contrib/parmetis/Lib
#                       PARMETIS_INCLUDE=-I$PARMETIS_INCLUDE_PATH
#                       PARMETIS_LIB="\$(EXTERNAL_LIBDIR)/libparmetis\$(libext)"
#       		  AC_SUBST(PARMETIS_INCLUDE)
#                       AC_SUBST(PARMETIS_LIB)
#                       AC_DEFINE(HAVE_PARMETIS, 1,
#       	                     [Flag indicating whether or not ParMetis is available])
#                       AC_MSG_RESULT(<<< Configuring library with ParMetis support >>>)
#       	          enableparmetis=yes
#                     ],
#                     [enableparmetis=no])
#     else
#       enableparmetis=no
#     fi  
#   else
#     enableparmetis=no
#   fi
# ])
