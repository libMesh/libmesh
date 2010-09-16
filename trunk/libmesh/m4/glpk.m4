dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl ----------------------------------------------------------------
dnl Certain parts of rbOOmit require GLPK, the GNU Linear Programming
dnl Kit.  By default we check for the GLPK installation files in 
dnl --with-glpk-include=xxx and --with-glpk-lib=yyy arguments provided to
dnl configure, or if those don't exist in $GLPK_INC and $GPLK_LIB
dnl directories, or in $GLPK_DIR/include and $GLPK_DIR/lib directories, or
dnl in /usr/local, or in /usr.
dnl ----------------------------------------------------------------

AC_DEFUN([CONFIGURE_GLPK], 
[
  dnl User-specific include path
  AC_ARG_WITH(glpk-include,
              AC_HELP_STRING([--with-glpk-include=PATH],[Specify the path for GLPK header files]),
              withglpkinc=$withval,
              withglpkinc=no)
	      
  dnl User-specific library path
  AC_ARG_WITH(glpk-lib,
              AC_HELP_STRING([--with-glpk-lib=PATH],[Specify the path for GLPK libs]),
              withglpklib=$withval,
              withglpklib=no)

  dnl Fall back on default paths to GLPK's include and lib files
  if (test $withglpkinc != no); then
    GLPK_INC="$withglpkinc"
  elif test "x$GLPK_INC" != x -a -f $GLPK_INC/glpk.h; then
    echo "Environment GLPK_INC=$GLPK_INC"
  elif test "x$GLPK_DIR" != x -a -f $GLPK_DIR/include/glpk.h; then
    GLPK_INC="$GLPK_DIR/include"
  elif [ -f /usr/local/include/glpk.h ]; then
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
  elif [ -f /usr/local/include/glpk.h ]; then
    GLPK_LIB="/usr/local/lib"
  else
    GLPK_LIB="/usr/lib"
  fi

  dnl Initialize Makefile/config.h substitution variables
  GLPK_INCLUDE=""
  GLPK_LIBRARY=""

  dnl Properly let the substitution variables
  if (test $enableglpk = yes); then
  
     dnl Check for existence of a header file in the specified location
     dnl AC_CHECK_FILE([$GLPK_INC/glpk.h], [glpkincFound="OK"], [glpkincFound="FAIL"])
     glpkincFound=no;
     AC_CHECK_HEADERS($GLPK_INC/glpk.h, glpkincFound=yes)

     if (test $glpkincFound = no); then
       AC_MSG_RESULT(GLPK header files not found!)
       enableglpk=no;
     fi

     dnl Discover the major and minor version numbers of GLPK by looking in
     dnl glpk.h.  This may eventually be useful for compiling against different
     dnl GLPK APIs...
     if (test $enableglpk = yes); then
       glpkmajor=`grep "define GLP_MAJOR_VERSION" $GLPK_INC/glpk.h | sed -e "s/#define GLP_MAJOR_VERSION[ ]*//g"`
       glpkminor=`grep "define GLP_MINOR_VERSION" $GLPK_INC/glpk.h | sed -e "s/#define GLP_MINOR_VERSION[ ]*//g"`
       glpkversion=$glpkmajor.$glpkminor
       AC_MSG_RESULT(<<< Configuring library with GLPK version $glpkversion support >>>)
     fi

     if (test $enableglpk = yes); then
       dnl Also Check for existence of required libraries.  This is not really the
       dnl right way to do it -- it's not portable to Macs, where .so's are called
       dnl .dylib's instead.
       dnl AC_CHECK_FILE($GLPK_LIB/glpk.so, [enableglpk=yes], [enableglpk=no])
       
       dnl AC_HAVE_LIBRARY (library, [action-if-found], [action-if-not-found], [other-libraries])
       dnl Note: Basically tries to compile a function which calls main().  

       dnl Save original value of LIBS, then append $GLPK_LIB
       old_LIBS="$LIBS"
       LIBS="$old_LIBS -L$GLPK_LIB"

       dnl Try to compile test prog to check for existence of GLPK libraries
       dnl AC_HAVE_LIBRARY uses the LIBS variable.
       AC_HAVE_LIBRARY([glpk], [enableglpk=yes], [enableglpk=no])
       
       dnl Reset $LIBS
       LIBS="$old_LIBS"
     fi
     
     dnl If both the header file and the required libs were found, continue.
     if (test $enableglpk = yes); then
       GLPK_INCLUDE="-I$GLPK_INC"
       GLPK_LIBRARY="\$(libmesh_RPATHFLAG)$GLPK_LIB -L$GLPK_LIB -lglpk"
       AC_DEFINE(HAVE_GLPK, 1, [Flag indicating whether the library will be compiled with GLPK support])
       AC_MSG_RESULT(<<< Configuring library with GLPK support >>>)
     fi
  fi

  dnl Substitute the substitution variables
  AC_SUBST(GLPK_INCLUDE)
  AC_SUBST(GLPK_LIBRARY)	
  AC_SUBST(enableglpk)
])
