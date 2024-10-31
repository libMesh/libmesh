# ----------------------------------------------------------------
# Test compilation against XDR with $XDRINCLUDES and $XDRLINKLIBS,
# ----------------------------------------------------------------
AC_DEFUN([TEST_XDR],
[
  old_CPPFLAGS="$CPPFLAGS"
  old_LIBS="$LIBS"
  CPPFLAGS="$CPPFLAGS $XDRINCLUDES"
  LIBS="$LIBS $XDRLINKLIBS"

  AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@include <stdio.h>
                                   @%:@include <rpc/rpc.h>
                                   @%:@include <rpc/xdr.h>],
                                  [
                                    XDR * xdr;
                                    FILE * fp;
                                    xdrstdio_create(xdr, fp, XDR_ENCODE);
                                  ])],
                 [
                   AC_MSG_RESULT(yes)
                   AC_DEFINE(HAVE_XDR, 1, [Flag indicating headers and libraries for XDR IO are available])
                   enablexdr=yes
                 ],
                 [
                   AC_MSG_RESULT(no)
                   enablexdr=no
                 ])

  dnl Reset flags after testing
  CPPFLAGS="$old_CPPFLAGS"
  LIBS="$old_LIBS"
])


# --------------------------------------------------------------------
# Attempt to find a working XDR configuration, based on user options
# if supplied or checking common defaults if not.
#
# Return the status of the attempt in $enablexdr.  Return the
# necessary flags in $XDRINCLUDES and $XDRLINKLIBS
# --------------------------------------------------------------------
AC_DEFUN([CONFIGURE_XDR],
[
  dnl If the user sets ${TIRPC_DIR} in their environment, then use that as the
  dnl default value for --with-xdr-include
  AS_IF([test "x${TIRPC_DIR}" != "x"], [withxdrinc_default=${TIRPC_DIR}], [withxdrinc_default=no])

  # Support for --with-xdr-include configure option
  AC_ARG_WITH(xdr-include,
              AS_HELP_STRING([--with-xdr-include=PATH],[Specify the path for RPC headers required by XDR]),
              withxdrinc=$withval,
              withxdrinc=$withxdrinc_default)

  # Support for --with-xdr-libdir configure option.  If the user does
  # not specify this option, we assume the desired RPC libraries are
  # already in the user's linker path.
  AC_ARG_WITH(xdr-libdir,
              AS_HELP_STRING([--with-xdr-libdir=/foo/bar],[Specify XDR library directory to be searched when linking]),
              withxdrlibdir=$withval,
              withxdrlibdir=no)

  # Support for --with-xdr-libname configure option.  If the user does
  # not specify this option, we do check "tirpc" since that seems to
  # be the most common way of linking "external" RPC libraries.
  AC_ARG_WITH(xdr-libname,
              AS_HELP_STRING([--with-xdr-libname=foo],[Specify name of the XDR library to be appended to -l when linking]),
              withxdrlibname=$withval,
              withxdrlibname=no)

  XDRINCLUDES=''
  XDRLINKLIBS=''
  AS_IF([test "x$withxdrlibname" != "xno"],
        [
          XDRLINKLIBS="-l$withxdrlibname"
        ])

  AS_IF([test "x$withxdrlibdir" != "xno"],
        [
          AS_IF([test "x$withxdrlibname" = "xno"],
                [AC_MSG_ERROR([xdr-libdir was set to $withxdrlibdir, but no xdr-libname was set])
                ])
          XDRLINKLIBS="${RPATHFLAG}$withxdrlibdir $XDRLINKLIBS"
        ])

  dnl If there is a user-provided XDR include/library provided, then we only
  dnl try compiling/linking against those, to avoid accidentally bringing in
  dnl any unwanted system RPC headers.
  AS_IF([test "x$withxdrinc" != "xno"],
        [
          XDRINCLUDES="-I$withxdrinc"

          AS_IF([test "x$withxdrlibdir" != "xno"],
                [ AC_MSG_CHECKING([for XDR support in $withxdrinc and $withxdrlibdir]) ],
                [ AC_MSG_CHECKING([for XDR support in $withxdrinc]) ])

          TEST_XDR

          dnl Check again, specifically for libtirpc, if the user
          dnl didn't specify a library name.
          AS_IF([test "x$withxdrlibname" = "xno"],
                [
                  AS_IF([test "x$withxdrlibdir" != "xno"],
                        [ AC_MSG_CHECKING([for XDR support in $withxdrinc and $withxdrlibdir with -ltirpc]) ],
                        [ AC_MSG_CHECKING([for XDR support in $withxdrinc with -ltirpc]) ])
                  XDRLINKLIBS="$XDRLINKLIBS -ltirpc"

                  TEST_XDR
               ])
        ],
        [
          dnl The user did not provide a custom XDR include/header path, so we
          dnl instead check whether the system/compiler has built-in support for
          dnl compiling/linking a code that includes xdr.h, followed by checking
          dnl other known locations.

          AS_IF([test "x$withxdrlibdir" != "xno"],
                [ AC_MSG_CHECKING([for XDR library in $withxdrlibdir]) ],
                [ AC_MSG_CHECKING([for built-in XDR support]) ])

          TEST_XDR

          dnl Check again for an explicitly named library.
          AS_IF([test "x$enablexdr" = "xno" && test "x$withxdrlibname" = "xno"],
                [
                  AC_MSG_CHECKING([for XDR support via libtirpc])
                  XDRLINKLIBS="-ltirpc"

                  TEST_XDR
               ])

          dnl Check again for headers in /usr/include/tirpc. (Fedora 28 does this.)
          AS_IF([test "x$enablexdr" = "xno" && test "x$withxdrlibdir" = "xno"],
                [
                  AC_MSG_CHECKING([for XDR support in /usr/include/tirpc])
                  XDRINCLUDES="-I/usr/include/tirpc"
                  XDRLINKLIBS="-ltirpc"

                  TEST_XDR
               ])
        ])

  dnl If nothing worked, don't report any flags that didn't work.
  AS_IF([test "x$enablexdr" = "xno"],
        [
          XDRINCLUDES=''
          XDRLINKLIBS=''
        ])
])
