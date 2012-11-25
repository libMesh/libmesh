# -------------------------------------------------------------
# fparser
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_FPARSER],
[
  AC_ARG_ENABLE([fparser],
                AC_HELP_STRING([--enable-fparser],
                               [build with C++ function parser support]),
		[case "${enableval}" in
		  yes)  enablefparser=yes ;;
		   no)  enablefparser=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-fparser) ;;
		 esac],
		 [enablefparser=$enableoptional])

  AC_ARG_WITH([fparser],
	       AC_HELP_STRING([--with-fparser=<release|none|devel>],
                              [Determine which version of the C++ function parser to use]),
	      [case "${withval}" in
		  release)  enablefparserdevel=no  ;;
		    devel)  enablefparserdevel=yes ;;
		     none)  enablefparser=no ;;
 		    *)  AC_MSG_ERROR(bad value ${withval} for --with-fparser) ;;
		 esac],
	      [enablefparserdevel=no])

	      
  if (test x$enablefparser = xyes); then
  
    
    AC_ARG_ENABLE(fparser-optimizer,
                  AC_HELP_STRING([--enable-fparser-optimizer],
                                 [use fparser optimization where possible]),
  		[case "${enableval}" in
  		  yes)  enablefparseroptimizer=yes ;;
  		   no)  enablefparseroptimizer=no ;;
   		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-fparser-optimizer) ;;
  		 esac],
  		 [enablefparseroptimizer=yes])
  
    # note - fparser optimization currently fails on OSX, to disable it regardless	
    case "${host_os}" in	
      *darwin* | *cygwin*) 	
        enablefparseroptimizer=no 	
        AC_MSG_RESULT(>>> Disabling fparser optimization on ${host_os} <<<)	
        ;;	
        *) ;;	
      esac
  			 
    # The FPARSER API is distributed with libmesh, so we don't have to guess
    # where it might be installed...

     AC_PROG_MKDIR_P
     AC_PROG_SED
     AC_PROG_YACC

     FPARSER_INCLUDE="-I\$(top_srcdir)/contrib/fparser"
     FPARSER_LIBRARY="\$(EXTERNAL_LIBDIR)/libfparser\$(libext)"
     AC_DEFINE(HAVE_FPARSER, 1, [Flag indicating whether the library will be compiled with FPARSER support])

      if (test $enablefparserdevel = yes); then
        AC_DEFINE(HAVE_FPARSER_DEVEL, 1, [Flag indicating whether FPARSER will build the full development version])
        AC_MSG_RESULT(<<< Configuring library with fparser support (development version) >>>)
      else
        AC_DEFINE(HAVE_FPARSER_DEVEL, 0, [Flag indicating whether FPARSER will build the full development version])
        AC_MSG_RESULT(<<< Configuring library with fparser support (release version) >>>)
      fi

  else
     FPARSER_INCLUDE=""
     FPARSER_LIBRARY=""
     enablefparser=no
  fi

  AC_SUBST(FPARSER_INCLUDE)
  AC_SUBST(FPARSER_LIBRARY)
  AC_SUBST(enablefparser)
  AC_SUBST(enablefparserdevel)

  AM_CONDITIONAL(FPARSER_RELEASE,              test x$enablefparserdevel = xno) 
  AM_CONDITIONAL(FPARSER_DEVEL,                test x$enablefparserdevel = xyes)
  AM_CONDITIONAL(FPARSER_NO_SUPPORT_OPTIMIZER, test x$enablefparseroptimizer = xno)
  AM_CONDITIONAL(FPARSER_SUPPORT_OPTIMIZER,    test x$enablefparseroptimizer = xyes)
])
