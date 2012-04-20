dnl -------------------------------------------------------------
dnl nemesis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NEMESIS], 
[
  AC_ARG_ENABLE(nemesis,
                AC_HELP_STRING([--enable-nemesis],
                               [build with NemesisII API support]),
		[case "${enableval}" in
		  yes)  enablenemesis=yes ;;
		   no)  enablenemesis=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-nemesis) ;;
		 esac],
		 [enablenemesis=$enableexodus]) # if unspecified, depend on exodus


  dnl Trump --enable-nemesis with --disable-mpi
  if (test "x$enablempi" = xno); then
    enablenemesis=no
  fi	
		 		
  dnl The NEMESIS API is distributed with libmesh, so we don't have to guess
  dnl where it might be installed...
  if (test $enablenemesis = yes); then
     NEMESIS_INCLUDE="-I\$(top_srcdir)/contrib/nemesis/Lib"
     NEMESIS_LIBRARY="\$(EXTERNAL_LIBDIR)/libnemesis\$(libext)"
     AC_DEFINE(HAVE_NEMESIS_API, 1, [Flag indicating whether the library will be compiled with Nemesis support])
     AC_MSG_RESULT(<<< Configuring library with Nemesis support >>>)
     libmesh_contrib_INCLUDES="$NEMESIS_INCLUDE $libmesh_contrib_INCLUDES"
     have_nemesis=yes
  else
     NEMESIS_INCLUDE=""
     NEMESIS_LIBRARY=""
     enablenemesis=no
     have_nemesis=no
  fi

  AC_SUBST(NEMESIS_INCLUDE)
  AC_SUBST(NEMESIS_LIBRARY)
  AC_SUBST(enablenemesis)

  AM_CONDITIONAL(LIBMESH_ENABLE_NEMESIS, test x$enablenemesis = xyes)
])
