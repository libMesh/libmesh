dnl -------------------------------------------------------------
dnl Metis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_METIS],
[
  AC_ARG_ENABLE(metis,
                AC_HELP_STRING([--enable-metis],
                               [build with Metis graph partitioning suppport]),
		[case "${enableval}" in
		  yes)  enablemetis=yes ;;
		   no)  enablemetis=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-metis) ;;
		 esac],
		 [enablemetis=$enableoptional])

 AC_ARG_WITH(metis,
             AC_HELP_STRING([--with-metis=<internal,PETSc>],
                            [metis to use.
			      interal: build from contrib.
			      PETSc: rely on PETSc]),
             [case "${withval}" in
                  internal)   build_metis=yes ;;
		  PETSc)      build_metis=no ;;
                      *) AC_MSG_ERROR(bad value ${withval} for --with-metis) ;;
                  esac],
                 [build_metis=yes])

  dnl The METIS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablemetis = yes); then
     METIS_INCLUDE="-I\$(top_srcdir)/contrib/metis/include"
     METIS_LIB="\$(EXTERNAL_LIBDIR)/libmetis\$(libext) \$(EXTERNAL_LIBDIR)/libGK\$(libext)"
     AC_DEFINE(HAVE_METIS, 1, [Flag indicating whether the library will be compiled with Metis support])
     AC_MSG_RESULT(<<< Configuring library with Metis support >>>)

     dnl look for thread-local storage
     AX_TLS
 else
     METIS_INCLUDE=""
     METIS_LIB=""
     enablemetis=no
  fi
  AM_CONDITIONAL(BUILD_METIS, test x$build_metis = xyes)

  AC_SUBST(METIS_INCLUDE)
  AC_SUBST(METIS_LIB)
])
