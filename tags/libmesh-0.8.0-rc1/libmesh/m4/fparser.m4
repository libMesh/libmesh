# -------------------------------------------------------------
# fparser
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_FPARSER],
[
  AC_ARG_WITH([fparser],
	       AC_HELP_STRING([--with-fparser=<release|none|devel>],
                              [Determine which version of the C++ function parser to use]),
	      [fparser_value="$withval"],
	      [fparser_value=release])

  if test "$fparser_value" == none; then
    enablefparser="no"
    enablefparserdevel="no"
  elif test "$fparser_value" == release; then
    enablefparser="yes"
    enablefparserdevel="no"
  elif test "$fparser_value" == devel; then
    enablefparser="yes"
    enablefparserdevel="yes"
  else
    enablefparser="yes"
    enablefparserdevel="no"
  fi

  # The FPARSER API is distributed with libmesh, so we don't have to guess
  # where it might be installed...
  if (test $enablefparser = yes); then
     AC_PROG_MKDIR_P
     AC_PROG_SED
     AC_PROG_YACC

    # note - currently skipping fparser optimization on all platforms
    # note - fparser optimization currently fails on OSX, to disable it regardless
#    case "${host_os}" in
#      *darwin* | *cygwin*)
#        enablefparseroptimizer=no
#        AC_MSG_RESULT(<<< Disabling fparser optimization on ${host_os} >>>)
#        ;;
#        *) ;;
#    esac

     FPARSER_INCLUDE="-I\$(top_srcdir)/contrib/fparser"
     FPARSER_LIBRARY="\$(EXTERNAL_LIBDIR)/libfparser\$(libext)"
     AC_DEFINE(HAVE_FPARSER, 1, [Flag indicating whether the library will be compiled with FPARSER support])

    # note - currently skipping fparser optimization on all platforms
#     if (test $enablefparseroptimizer = yes); then
#       AC_MSG_RESULT(<<< Configuring library with fparser support (with optimizer) >>>)
#     else
#       AC_MSG_RESULT(<<< Configuring library with fparser support (without optimizer) >>>)
#     fi

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
     AC_MSG_RESULT(<<< Disabling fparser support >>>)
  fi

  AC_SUBST(FPARSER_INCLUDE)
  AC_SUBST(FPARSER_LIBRARY)
  AC_SUBST(enablefparser)
  AC_SUBST(enablefparserdevel)

#  AM_CONDITIONAL(FPARSER_NO_SUPPORT_OPTIMIZER, test x$enablefparseroptimizer = xno)
#  AM_CONDITIONAL(FPARSER_SUPPORT_OPTIMIZER,    test x$enablefparseroptimizer = xyes)
])
