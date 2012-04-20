# -------------------------------------------------------------
# fparser
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_FPARSER], 
[
  AC_ARG_ENABLE(fparser,
                AC_HELP_STRING([--enable-fparser],
                               [build with fparser, from Juha Nieminen, Joel Yliluoma]),
		[case "${enableval}" in
		  yes)  enablefparser=yes ;;
		   no)  enablefparser=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-fparser) ;;
		 esac],
		 [enablefparser=$enableoptional])

  AC_ARG_ENABLE(fparser-optimizer,
                AC_HELP_STRING([--enable-fparser-optimizer],
                               [use fparser optimization where possible]),
		[case "${enableval}" in
		  yes)  enablefparseroptimizer=yes ;;
		   no)  enablefparseroptimizer=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-fparser-optimizer) ;;
		 esac],
		 [enablefparseroptimizer=yes])

  # The FPARSER API is distributed with libmesh, so we don't have to guess
  # where it might be installed...
  if (test $enablefparser = yes); then
     AC_PROG_MKDIR_P
     AC_PROG_SED
     AC_PROG_YACC

    # note - fparser optimization currently fails on OSX, to disable it regardless
    case "${host_os}" in
      *darwin*) 
        enablefparseroptimizer=no 
        AC_MSG_RESULT(>>> Disabling fparser optimization on ${host_os} <<<)
        ;;
        *) ;;
    esac	

     FPARSER_INCLUDE="-I\$(top_srcdir)/contrib/fparser"
     FPARSER_LIBRARY="\$(EXTERNAL_LIBDIR)/libfparser\$(libext)"
     AC_DEFINE(HAVE_FPARSER, 1, [Flag indicating whether the library will be compiled with FPARSER support])
     libmesh_contrib_INCLUDES="$FPARSER_INCLUDE $libmesh_contrib_INCLUDES"
     if (test $enablefparseroptimizer = yes); then
       AC_MSG_RESULT(<<< Configuring library with fparser support (with optimizer) >>>)
     else
       AC_MSG_RESULT(<<< Configuring library with fparser support (without optimizer) >>>)
     fi	
  else
     FPARSER_INCLUDE=""
     FPARSER_LIBRARY=""
     enablefparser=no
  fi

  AC_SUBST(FPARSER_INCLUDE)
  AC_SUBST(FPARSER_LIBRARY)	
  AC_SUBST(enablefparser)
		 		 
  AM_CONDITIONAL(FPARSER_NO_SUPPORT_OPTIMIZER, test x$enablefparseroptimizer = xno)
  AM_CONDITIONAL(FPARSER_SUPPORT_OPTIMIZER,    test x$enablefparseroptimizer = xyes)
  AM_CONDITIONAL(LIBMESH_ENABLE_FPARSER, test x$enablefparser = xyes)
])
