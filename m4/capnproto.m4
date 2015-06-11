# ----------------------------------------------------------------
# The Reduced Basis code uses Cap'n Proto to write training data to disk.
# By default we check for the Capn'n Proto installation files in
# --with-capnproto-include=xxx and --with-capnproto-lib=yyy arguments
# provided to configure, or if those don't exist in $CAPNPROTO_INC and
# $CAPNPROTO_LIB directories, or in $CAPNPROTO_DIR/include and
# $CAPNPROTO_DIR/lib directories, or in /usr/local, or in /usr.
# ----------------------------------------------------------------

AC_DEFUN([CONFIGURE_CAPNPROTO],
[
  AC_ARG_ENABLE(capnproto,
                AS_HELP_STRING([--disable-capnproto],
                               [build without Cap'n Proto support]),
		[case "${enableval}" in
		  yes)  enablecapnproto=yes ;;
		   no)  enablecapnproto=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-capnproto) ;;
		 esac],
		 [enablecapnproto=$enableoptional])


  if (test $enablecapnproto = yes); then

    # User-specific include path
    AC_ARG_WITH(capnproto-include,
                AS_HELP_STRING([--with-capnproto-include=PATH],[Specify the path for CAPNPROTO header files]),
                withcapnprotoinc=$withval,
                withcapnprotoinc=no)

    # User-specific library path
    AC_ARG_WITH(capnproto-lib,
                AS_HELP_STRING([--with-capnproto-lib=PATH],[Specify the path for CAPNPROTO libs]),
                withcapnprotolib=$withval,
                withcapnprotolib=no)

    # Use CAPNPROTO_DIR/include if it exists.
    if (test $withcapnprotoinc != no); then
      CAPNPROTO_INC="$withcapnprotoinc"
    elif test "x$CAPNPROTO_DIR" != x -a -f $CAPNPROTO_DIR/include/capnproto.h; then
      CAPNPROTO_INC="$CAPNPROTO_DIR/include"
    else
      CAPNPROTO_INC=""
    fi

    # Use CAPNPROTO_DIR/lib if it exists.
    if (test $withcapnprotolib != no); then
      CAPNPROTO_LIB="$withcapnprotolib"
    elif test "x$CAPNPROTO_DIR" != x; then
      CAPNPROTO_LIB="$CAPNPROTO_DIR/lib"
    else
      CAPNPROTO_LIB=""
    fi

    # Initialize Makefile/config.h substitution variables
    CAPNPROTO_INCLUDE=""
    CAPNPROTO_LIBRARY=""

    # Properly let the substitution variables
    if (test $enablecapnproto = yes); then

       # Check for existence of a header file in the specified location
       capnprotoincFound=no;
       AC_CHECK_HEADERS($CAPNPROTO_INC/capnproto.h, capnprotoincFound=yes)

       if (test $capnprotoincFound = no); then
         AC_MSG_RESULT(CAPNPROTO header files not found!)
         enablecapnproto=no;
       fi

       # If both the header file and the required libs were found, continue.
       if (test x$enablecapnproto = xyes); then
         CAPNPROTO_INCLUDE="-I$CAPNPROTO_INC"
         CAPNPROTO_LIBRARY="-L$CAPNPROTO_LIB -lcapnproto -lkj"

          # add the CAPNPROTO_LIB dir to the linker run path if it is a directory
          if (test "x$RPATHFLAG" != "x" -a -d $CAPNPROTO_LIB); then
            CAPNPROTO_LIBRARY="${RPATHFLAG}${CAPNPROTO_LIB} $CAPNPROTO_LIBRARY"
          fi

         AC_DEFINE(HAVE_CAPNPROTO, 1, [Flag indicating whether the library will be compiled with CAPNPROTO support])
         AC_MSG_RESULT(<<< Configuring library with CAPNPROTO support >>>)
       fi
    fi
  fi

  # Substitute the substitution variables
  AC_SUBST(CAPNPROTO_INCLUDE)
  AC_SUBST(CAPNPROTO_LIBRARY)
])
