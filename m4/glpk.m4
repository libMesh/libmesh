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
                [case "${enableval}" in
                  yes)  enableglpk=yes ;;
                  no)  enableglpk=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-glpk) ;;
                esac],
                [enableglpk=$enableoptional])

  if (test $enableglpk = yes); then

    # User-specific include path
    AC_ARG_WITH(glpk-include,
                AS_HELP_STRING([--with-glpk-include=PATH],[Specify the path for GLPK header files]),
                withglpkinc=$withval,
                withglpkinc=no)

    # User-specific library path
    AC_ARG_WITH(glpk-lib,
                AS_HELP_STRING([--with-glpk-lib=PATH],[Specify the path for GLPK libs]),
                withglpklib=$withval,
                withglpklib=no)

    # Fall back on default paths to GLPK's include and lib files
    if (test $withglpkinc != no); then
      GLPK_INC="$withglpkinc"
    elif test "x$GLPK_INC" != x -a -f $GLPK_INC/glpk.h; then
      echo "Environment GLPK_INC=$GLPK_INC"
    elif test "x$GLPK_DIR" != x -a -f $GLPK_DIR/include/glpk.h; then
      GLPK_INC="$GLPK_DIR/include"
    elif test -f /usr/include/glpk/glpk.h; then # RHEL6 puts glpk here
      GLPK_INC="/usr/include/glpk"
    elif test -f /usr/local/include/glpk.h; then
      GLPK_INC="/usr/local/include"
    else
      GLPK_INC="/usr/include"
    fi

    if (test $withglpklib != no); then
      GLPK_LIB="$withglpklib"
    elif test "x$GLPK_LIB" != x; then
      echo "Environment GLPK_LIB=$GLPK_INC"
    elif test "x$GLPK_DIR" != x; then
      GLPK_LIB="$GLPK_DIR/lib"
    elif test -f /usr/include/glpk/glpk.h; then # RHEL6 puts glpk here
      if test -f /usr/lib64/libglpk.so -o -f /usr/lib64/libglpk.a; then
        GLPK_LIB="/usr/lib64"
      elif test -f /usr/lib/libglpk.so -o -f /usr/lib/libglpk.a; then
        GLPK_LIB="/usr/lib"
      fi
    elif test -f /usr/local/include/glpk.h; then
      GLPK_LIB="/usr/local/lib"
    else
      GLPK_LIB="/usr/lib"
    fi

    # Initialize Makefile/config.h substitution variables
    GLPK_INCLUDE=""
    GLPK_LIBRARY=""

    # Properly let the substitution variables
    if (test $enableglpk = yes); then

       # Check for existence of a header file in the specified location
       glpkincFound=no;
       AC_CHECK_HEADERS($GLPK_INC/glpk.h, glpkincFound=yes)

       if (test $glpkincFound = no); then
         AC_MSG_RESULT(GLPK header files not found!)
         enableglpk=no;
       fi

       # Discover the major and minor version numbers of GLPK by looking in
       # glpk.h.  This may eventually be useful for compiling against different
       # GLPK APIs...
       if (test $enableglpk = yes); then
         glpkmajor=`grep "define GLP_MAJOR_VERSION" $GLPK_INC/glpk.h | sed -e "s/#define GLP_MAJOR_VERSION[ ]*//g"`
         glpkminor=`grep "define GLP_MINOR_VERSION" $GLPK_INC/glpk.h | sed -e "s/#define GLP_MINOR_VERSION[ ]*//g"`
         glpkversion=$glpkmajor.$glpkminor
         AC_MSG_RESULT(<<< Configuring library with GLPK version $glpkversion support >>>)
       fi

       if (test $enableglpk = yes); then
         # Also Check for existence of required libraries.

         # Save original value of LIBS, then append $GLPK_LIB
         old_LIBS="$LIBS"
         LIBS="$old_LIBS -L$GLPK_LIB"

         # AC_CHECK_LIB tries to link a test code that calls a certain function.
         # AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
         AC_CHECK_LIB([glpk], main, [enableglpk=yes], [enableglpk=no])

         # Reset $LIBS
         LIBS="$old_LIBS"
       fi

       # If both the header file and the required libs were found, continue.
       if (test x$enableglpk = xyes); then
         GLPK_INCLUDE="-I$GLPK_INC"
         GLPK_LIBRARY="-L$GLPK_LIB -lglpk"
         if (test "x$RPATHFLAG" != "x" -a -d $GLPK_LIB); then # add the GLPK_LIB to the linker run path, if it is a directory
           if (test "$GLPK_LIB" != "/usr/lib" -a "$GLPK_LIB" != "/usr/lib64"); then
             GLPK_LIBRARY="${RPATHFLAG}${GLPK_LIB} $GLPK_LIBRARY"
           fi
         fi
         AC_DEFINE(HAVE_GLPK, 1, [Flag indicating whether the library will be compiled with GLPK support])
         AC_MSG_RESULT(<<< Configuring library with GLPK support >>>)
       fi
    fi
  fi


  # Substitute the substitution variables
  AC_SUBST(GLPK_INCLUDE)
  AC_SUBST(GLPK_LIBRARY)
])
