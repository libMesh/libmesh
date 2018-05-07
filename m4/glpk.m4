# ----------------------------------------------------------------
# Certain parts of rbOOmit require GLPK, the GNU Linear Programming
# Kit.  By default we check for the GLPK installation files in
# --with-glpk-include=xxx and --with-glpk-lib=yyy arguments provided to
# configure, or if those don't exist in $GLPK_INC and $GPLK_LIB
# directories, or in $GLPK_DIR/include and $GLPK_DIR/lib directories, or
# in /usr/local, or in /usr.
# ----------------------------------------------------------------

AC_DEFUN([CONFIGURE_GLPK],
[
  AC_ARG_ENABLE(glpk,
                AS_HELP_STRING([--disable-glpk],
                               [build without GLPK support]),
                [AS_CASE("${enableval}",
                         [yes], [enableglpk=yes],
                         [no],  [enableglpk=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-glpk)])],
                [enableglpk=$enableoptional])

  AS_IF([test "x$enableglpk" = "xyes"],
        [
          dnl User-specific include path
          AC_ARG_WITH(glpk-include,
                      AS_HELP_STRING([--with-glpk-include=PATH],[Specify the path for GLPK header files]),
                      withglpkinc=$withval,
                      withglpkinc=no)

          dnl User-specific library path
          AC_ARG_WITH(glpk-lib,
                      AS_HELP_STRING([--with-glpk-lib=PATH],[Specify the path for GLPK libs]),
                      withglpklib=$withval,
                      withglpklib=no)

          dnl Fall back on default paths to GLPK's include and lib files
          AS_IF([test "x$withglpkinc" != "xno"],                                [GLPK_INC="$withglpkinc"],
                [test "x$GLPK_INC" != "x" && test -f $GLPK_INC/glpk.h],         [AS_ECHO(["Environment GLPK_INC=$GLPK_INC"])],
                [test "x$GLPK_DIR" != "x" && test -f $GLPK_DIR/include/glpk.h], [GLPK_INC="$GLPK_DIR/include"],
                [test -f /usr/include/glpk/glpk.h],                             [GLPK_INC="/usr/include/glpk"],
                [test -f /usr/local/include/glpk.h],                            [GLPK_INC="/usr/local/include"],
                [GLPK_INC="/usr/include"])

          AS_IF([test "x$withglpklib" != "xno"],     [GLPK_LIB="$withglpklib"],
                [test "x$GLPK_LIB" != "x"],          [AS_ECHO(["Environment GLPK_LIB=$GLPK_INC"])],
                [test "x$GLPK_DIR" != "x"],          [GLPK_LIB="$GLPK_DIR/lib"],
                [test -f /usr/include/glpk/glpk.h],  [AS_IF([test -f /usr/lib64/libglpk.so || test -f /usr/lib64/libglpk.a], [GLPK_LIB="/usr/lib64"],
                                                            [test -f /usr/lib/libglpk.so || test -f /usr/lib/libglpk.a], [GLPK_LIB="/usr/lib"])],
                [test -f /usr/local/include/glpk.h], [GLPK_LIB="/usr/local/lib"],
                [GLPK_LIB="/usr/lib"])

          dnl Initialize Makefile/config.h substitution variables
          GLPK_INCLUDE=""
          GLPK_LIBRARY=""

          dnl Properly let the substitution variables
          AS_IF([test "x$enableglpk" = "xyes"],
                [
                  dnl Check for existence of a header file in the specified location
                  glpkincFound=no;
                  AC_CHECK_HEADERS($GLPK_INC/glpk.h, glpkincFound=yes)

                  AS_IF([test "x$glpkincFound" = "xno"],
                        [
                          AC_MSG_RESULT(GLPK header files not found!)
                          enableglpk=no
                        ])

                  dnl Discover the major and minor version numbers of GLPK by looking in
                  dnl glpk.h.  This may eventually be useful for compiling against different
                  dnl GLPK APIs...
                  AS_IF([test "x$enableglpk" = "xyes"],
                        [
                          glpkmajor=`grep "define GLP_MAJOR_VERSION" $GLPK_INC/glpk.h | sed -e "s/#define GLP_MAJOR_VERSION[ ]*//g"`
                          glpkminor=`grep "define GLP_MINOR_VERSION" $GLPK_INC/glpk.h | sed -e "s/#define GLP_MINOR_VERSION[ ]*//g"`
                          glpkversion=$glpkmajor.$glpkminor
                          AC_MSG_RESULT(<<< Configuring library with GLPK version $glpkversion support >>>)
                        ])

                  AS_IF([test "x$enableglpk" = "xyes"],
                        [
                          dnl Also Check for existence of required libraries.
                          dnl Save original value of LIBS, then append $GLPK_LIB
                          old_LIBS="$LIBS"
                          LIBS="$old_LIBS -L$GLPK_LIB"

                          dnl AC_CHECK_LIB tries to link a test code that calls a certain function.
                          dnl AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
                          AC_CHECK_LIB([glpk], main, [enableglpk=yes], [enableglpk=no])

                          dnl Reset $LIBS
                          LIBS="$old_LIBS"
                       ])

                  dnl If both the header file and the required libs were found, continue.
                  AS_IF([test "x$enableglpk" = "xyes"],
                        [
                          GLPK_INCLUDE="-I$GLPK_INC"
                          GLPK_LIBRARY="-L$GLPK_LIB -lglpk"
                          dnl add the GLPK_LIB to the linker run path, if it is a directory
                          AS_IF([test "x$RPATHFLAG" != "x" && test -d $GLPK_LIB],
                                [
                                  AS_IF([test "$GLPK_LIB" != "/usr/lib" && test "$GLPK_LIB" != "/usr/lib64"],
                                        [GLPK_LIBRARY="${RPATHFLAG}${GLPK_LIB} $GLPK_LIBRARY"])
                                ])
                          AC_DEFINE(HAVE_GLPK, 1, [Flag indicating whether the library will be compiled with GLPK support])
                          AC_MSG_RESULT(<<< Configuring library with GLPK support >>>)
                        ])
                ])
        ])

  dnl Substitute the substitution variables
  AC_SUBST(GLPK_INCLUDE)
  AC_SUBST(GLPK_LIBRARY)
])
