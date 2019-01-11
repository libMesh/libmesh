AC_DEFUN([AM_PATH_CPPUNIT],
[
  dnl Defaults that might work if cppunit headers are in /usr/include
  dnl and libraries are in /usr/lib, i.e. standard installation locations.
  CPPUNIT_CFLAGS=
  CPPUNIT_LIBS=-lcppunit

  dnl Check for the cppunit-config program, and if it exists, use it
  dnl to set compiler and linker flags.
  AC_CHECK_PROG(CPPUNIT_CONFIG, cppunit-config, cppunit-config, none, $PATH)
  AS_IF([test "x$CPPUNIT_CONFIG" = "xcppunit-config"],
        [
          CPPUNIT_CFLAGS=`$CPPUNIT_CONFIG --cflags`
          CPPUNIT_LIBS=`$CPPUNIT_CONFIG --libs`
        ])

  dnl User can override cppunit-config values with explicit values for:
  dnl --with-cppunit-include
  dnl --with-cppunit-lib
  AC_ARG_WITH(cppunit-include,
              AS_HELP_STRING([--with-cppunit-include=PATH],
                             [Specify a path for cppunit header files]),
              CPPUNIT_CFLAGS="-I$withval")

  AC_ARG_WITH(cppunit-lib,
              AS_HELP_STRING([--with-cppunit-lib=PATH],
                             [Specify a path for cppunit libs]),
              CPPUNIT_LIBS="-L$withval -lcppunit")

  AC_MSG_CHECKING(whether we can build a trivial CppUnit program)
  AC_LANG_PUSH([C++])
  saveCXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$saveCXXFLAGS $CPPUNIT_CFLAGS"
  saveLIBS="$LIBS"
  LIBS="$CPPUNIT_LIBS $saveLIBS"

  AC_LINK_IFELSE([AC_LANG_SOURCE([[
  @%:@include <cppunit/ui/text/TestRunner.h>
  int main(int argc, char **argv)
  {
    CppUnit::TextUi::TestRunner runner;

    if (runner.run())
      return 0;

    return 1;
  }
  ]])],[
    AC_MSG_RESULT(yes)
  ],[
    AC_MSG_RESULT(no)
    CPPUNIT_CFLAGS=
    CPPUNIT_LIBS=
    enablecppunit=no
  ])

  AC_LANG_POP
  LIBS="$saveLIBS"
  CXXFLAGS="$saveCXXFLAGS"

  AC_SUBST(CPPUNIT_CFLAGS)
  AC_SUBST(CPPUNIT_LIBS)
])
