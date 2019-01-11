AC_DEFUN([AM_PATH_CPPUNIT],
[
  dnl Defaults that might work if cppunit headers are in /usr/include
  dnl and libraries are in /usr/lib, i.e. standard installation locations.
  CPPUNIT_CFLAGS=
  CPPUNIT_LIBS=-lcppunit

  dnl User can specify --with-cppunit-include to specify path to cppunit headers.
  AC_ARG_WITH(cppunit-include,
              AS_HELP_STRING([--with-cppunit-include=PATH],
                             [Specify a path for cppunit header files]),
              CPPUNIT_CFLAGS="-I$withval",
              CPPUNIT_CFLAGS="")

  dnl User can specify --with-cppunit-lib to specify path to cppunit libs.
  AC_ARG_WITH(cppunit-lib,
              AS_HELP_STRING([--with-cppunit-lib=PATH],
                             [Specify a path for cppunit libs]),
              CPPUNIT_LIBS="-L$withval -lcppunit",
              CPPUNIT_LIBS="-lcppunit")

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
