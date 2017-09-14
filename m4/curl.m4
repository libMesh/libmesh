# The CURL API allows one to interact with URLs programmatically.
AC_DEFUN([CONFIGURE_CURL],
[
  AC_ARG_ENABLE(curl,
                AS_HELP_STRING([--enable-curl],
                               [link against libcurl, for using the cURL API]),
                [case "${enableval}" in
                  yes)  enablecurl=yes ;;
                  no)  enablecurl=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-curl) ;;
                esac],
                [enablecurl=no])

  if (test $enablecurl = yes); then

    # User-specific include path
    AC_ARG_WITH(curl-include,
                AS_HELP_STRING([--with-curl-include=PATH],[Specify the path for CURL header files]),
                withcurlinc=$withval,
                withcurlinc=no)

    # User-specific library path
    AC_ARG_WITH(curl-lib,
                AS_HELP_STRING([--with-curl-lib=PATH],[Specify the path for CURL libs]),
                withcurllib=$withval,
                withcurllib=no)

    # Start with reasonable defaults and let user-specified values
    # override them.  We don't support allowing user environment vars
    # to set these, as it just makes things too complicated and
    # error-prone.
    CURL_INC="/usr/include/curl"
    CURL_LIB="/usr/lib"

    # If the user specified --with-curl-include=foo, use that...
    if (test $withcurlinc != no); then
      CURL_INC="$withcurlinc"
    fi

    # If the user specified --with-curl-lib=foo, use that...
    if (test $withcurllib != no); then
      CURL_LIB="$withcurllib"
    fi

    # Initialize eventual Makefile/config.h substitution variables
    CURL_INCLUDE=""
    CURL_LIBRARY=""

    # Actually test for the existence of headers and libs.
    if (test $enablecurl = yes); then

       # Check for existence of a header file in the specified location
       curlincFound=no;
       AC_CHECK_HEADERS($CURL_INC/curl.h, curlincFound=yes)

       if (test $curlincFound = no); then
         AC_MSG_RESULT([CURL header files not found!])
         enablecurl=no;
       fi

       # Curl versions are defined in curlver.h (not curl.h!) and look like this:
       #define LIBCURL_VERSION_MAJOR 7
       #define LIBCURL_VERSION_MINOR 37
       #define LIBCURL_VERSION_PATCH 1
       if (test $enablecurl = yes); then
         curlmajor=`grep "define LIBCURL_VERSION_MAJOR" $CURL_INC/curlver.h | sed -e "s/#define LIBCURL_VERSION_MAJOR[ ]*//g"`
         curlminor=`grep "define LIBCURL_VERSION_MINOR" $CURL_INC/curlver.h | sed -e "s/#define LIBCURL_VERSION_MINOR[ ]*//g"`
         curlversion=$curlmajor.$curlminor
         AC_MSG_RESULT([<<< Configuring library with CURL version $curlversion support >>>])
       fi

       if (test $enablecurl = yes); then
         # Also Check for existence of required libraries.

         # Save original value of LIBS, then append $CURL_LIB
         old_LIBS="$LIBS"
         LIBS="$old_LIBS -L$CURL_LIB"

         # AC_CHECK_LIB tries to link a test code that calls a certain function.
         # AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
         AC_CHECK_LIB([curl], main, [enablecurl=yes], [enablecurl=no])

         # Reset $LIBS
         LIBS="$old_LIBS"
       fi

       # If both the header file and the required libs were found, continue.
       if (test x$enablecurl = xyes); then
         CURL_INCLUDE="-I$CURL_INC"
         CURL_LIBRARY="-L$CURL_LIB -lcurl"
         if (test "x$RPATHFLAG" != "x" -a -d $CURL_LIB); then # add the CURL_LIB to the linker run path, if it is a directory
           if (test "$CURL_LIB" != "/usr/lib" -a "$CURL_LIB" != "/usr/lib64"); then
             CURL_LIBRARY="${RPATHFLAG}${CURL_LIB} $CURL_LIBRARY"
           fi
         fi
         AC_DEFINE(HAVE_CURL, 1, [Flag indicating whether the library will be compiled with CURL support])
         AC_MSG_RESULT([<<< Configuring library with CURL support >>>])
       fi
    fi
  fi


  # Substitute the substitution variables
  AC_SUBST(CURL_INCLUDE)
  AC_SUBST(CURL_LIBRARY)
])
