# The CURL API allows one to interact with URLs programmatically.
AC_DEFUN([CONFIGURE_CURL],
[
  AC_ARG_ENABLE(curl,
                AS_HELP_STRING([--enable-curl],
                               [link against libcurl, for using the cURL API]),
                [AS_CASE("${enableval}",
                         [yes], [enablecurl=yes],
                         [no],  [enablecurl=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-curl)])],
                [enablecurl=no])

  AS_IF([test "x$enablecurl" = "xyes"],
        [
          dnl User-specific include path
          AC_ARG_WITH(curl-include,
                      AS_HELP_STRING([--with-curl-include=PATH],[Specify the path for CURL header files]),
                      withcurlinc=$withval,
                      withcurlinc=no)

          dnl User-specific library path
          AC_ARG_WITH(curl-lib,
                      AS_HELP_STRING([--with-curl-lib=PATH],[Specify the path for CURL libs]),
                      withcurllib=$withval,
                      withcurllib=no)

          dnl Start with reasonable defaults and let user-specified values
          dnl override them.  We don't support allowing user environment vars
          dnl to set these, as it just makes things too complicated and
          dnl error-prone.
          CURL_INC="/usr/include/curl"
          CURL_LIB="/usr/lib"

          dnl If the user specified --with-curl-include=foo, use that...
          AS_IF([test "x$withcurlinc" != "xno"], [CURL_INC="$withcurlinc"])

          dnl If the user specified --with-curl-lib=foo, use that...
          AS_IF([test "x$withcurllib" != "xno"], [CURL_LIB="$withcurllib"])

          # Initialize eventual Makefile/config.h substitution variables
          CURL_INCLUDE=""
          CURL_LIBRARY=""

          # Actually test for the existence of headers and libs.
          AS_IF([test "x$enablecurl" = "xyes"],
                [
                  dnl Check for existence of a header file in the specified location
                  curlincFound=no;
                  AC_CHECK_HEADERS($CURL_INC/curl.h, curlincFound=yes)

                  AS_IF([test "x$curlincFound" = "xno"],
                        [
                          AC_MSG_RESULT([CURL header files not found!])
                          enablecurl=no
                        ])

                  dnl Curl versions are defined in curlver.h (not curl.h!) and look like this:
                  #define LIBCURL_VERSION_MAJOR 7
                  #define LIBCURL_VERSION_MINOR 37
                  #define LIBCURL_VERSION_PATCH 1
                  AS_IF([test "x$enablecurl" = "xyes"],
                        [
                          curlmajor=`grep "define LIBCURL_VERSION_MAJOR" $CURL_INC/curlver.h | sed -e "s/#define LIBCURL_VERSION_MAJOR[ ]*//g"`
                          curlminor=`grep "define LIBCURL_VERSION_MINOR" $CURL_INC/curlver.h | sed -e "s/#define LIBCURL_VERSION_MINOR[ ]*//g"`
                          curlversion=$curlmajor.$curlminor
                          AC_MSG_RESULT([<<< Configuring library with CURL version $curlversion support >>>])
                        ])

                  AS_IF([test "x$enablecurl" = "xyes"],
                        [
                          dnl Also Check for existence of required libraries.
                          dnl Save original value of LIBS, then append $CURL_LIB
                          old_LIBS="$LIBS"
                          LIBS="$old_LIBS -L$CURL_LIB"

                          dnl AC_CHECK_LIB tries to link a test code that calls a certain function.
                          dnl AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
                          AC_CHECK_LIB([curl], main, [enablecurl=yes], [enablecurl=no])

                          dnl Reset $LIBS
                          LIBS="$old_LIBS"
                        ])

                  dnl If both the header file and the required libs were found, continue.
                  AS_IF([test "x$enablecurl" = "xyes"],
                        [
                          CURL_INCLUDE="-I$CURL_INC"
                          CURL_LIBRARY="-L$CURL_LIB -lcurl"
                          AS_IF([test "x$RPATHFLAG" != "x" && test -d $CURL_LIB],
                                [AS_IF([test "$CURL_LIB" != "/usr/lib" && test "$CURL_LIB" != "/usr/lib64"], [CURL_LIBRARY="${RPATHFLAG}${CURL_LIB} $CURL_LIBRARY"])])
                          AC_DEFINE(HAVE_CURL, 1, [Flag indicating whether the library will be compiled with CURL support])
                          AC_MSG_RESULT([<<< Configuring library with CURL support >>>])
                        ])
                ])
        ])

  dnl Substitute the substitution variables
  AC_SUBST(CURL_INCLUDE)
  AC_SUBST(CURL_LIBRARY)
])
