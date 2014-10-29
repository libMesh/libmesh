# -------------------------------------------------------------
# fparser
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_FPARSER],
[
  AC_ARG_ENABLE([fparser],
                AS_HELP_STRING([--disable-fparser],
                               [build without C++ function parser support]),
                [case "${enableval}" in
                  yes)  enablefparser=yes ;;
                   no)  enablefparser=no ;;
                    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-fparser) ;;
                 esac],
                [enablefparser=$enableoptional])

  AC_ARG_WITH([fparser],
               AS_HELP_STRING([--with-fparser=<release|none|devel>],
                              [Determine which version of the C++ function parser to use]),
               [case "${withval}" in
                   release)  enablefparserdevel=no  ;;
                     devel)  enablefparserdevel=yes ;;
                      none)  enablefparser=no ;;
                         *)  AC_MSG_ERROR(bad value ${withval} for --with-fparser) ;;
                esac],
               [enablefparserdevel=no])


  if (test x$enablefparser = xyes); then


    AC_ARG_ENABLE(fparser-debugging,
                  AS_HELP_STRING([--enable-fparser-debugging],
                                 [Build fparser with bytecode debugging functions]),
                  [case "${enableval}" in
                    yes)  enablefparserdebugging=yes ;;
                     no)  enablefparserdebugging=no ;;
                      *)  AC_MSG_ERROR(bad value ${enableval} for --enable-fparser-debugging) ;;
                   esac],
                  [enablefparserdebugging=no])

    AC_ARG_ENABLE(fparser-optimizer,
                  AS_HELP_STRING([--disable-fparser-optimizer],
                                 [do not optimize parsed functions]),
                  [case "${enableval}" in
                    yes)  enablefparseroptimizer=yes ;;
                     no)  enablefparseroptimizer=no ;;
                      *)  AC_MSG_ERROR(bad value ${enableval} for --enable-fparser-optimizer) ;;
                   esac],
                  [enablefparseroptimizer=yes])

    # note - fparser optimization may fail on cygwin (please test), so disable it regardless
    case "${host_os}" in
      *cygwin*)
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

      AC_SEARCH_LIBS([dlopen], [dl dld], [
        AC_DEFINE(HAVE_FPARSER_JIT, 1, [Flag indicating whether FPARSER will be built with JIT compilation enabled])
        AC_MSG_RESULT(<<< Configuring library with fparser JIT compilation support >>>)
        enablefparserjit=yes
      ], [
        AC_DEFINE(HAVE_FPARSER_JIT, 0, [Flag indicating whether FPARSER will be built with JIT compilation enabled])
        AC_MSG_RESULT(<<< dlopen() not found, configuring library without fparser JIT compilation support >>>)
        enablefparserjit=no
      ])

      # This define in libmesh_config.h is used internally in fparser.hh and various source files
      if (test $enablefparserdebugging = yes); then
        AC_DEFINE(FPARSER_SUPPORT_DEBUGGING, 1, [Enable fparser debugging functions])
        AC_MSG_RESULT(<<< Configuring library with fparser debugging functions >>>)
      fi

  else
     FPARSER_INCLUDE=""
     FPARSER_LIBRARY=""
     enablefparser=no
  fi

  AC_SUBST(FPARSER_INCLUDE)
  AC_SUBST(FPARSER_LIBRARY)

  AM_CONDITIONAL(FPARSER_RELEASE,              test x$enablefparserdevel = xno)
  AM_CONDITIONAL(FPARSER_DEVEL,                test x$enablefparserdevel = xyes)
  AM_CONDITIONAL(FPARSER_NO_SUPPORT_OPTIMIZER, test x$enablefparseroptimizer = xno)
  AM_CONDITIONAL(FPARSER_SUPPORT_OPTIMIZER,    test x$enablefparseroptimizer = xyes)
  AM_CONDITIONAL(FPARSER_SUPPORT_DEBUGGING,    test x$enablefparserdebugging = xyes)
  AM_CONDITIONAL(FPARSER_SUPPORT_JIT,    test x$enablefparserjit = xyes)
])
