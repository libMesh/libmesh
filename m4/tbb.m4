dnl -------------------------------------------------------------
dnl Threading Building Blocks
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TBB],
[
  AC_ARG_ENABLE(tbb,
                AS_HELP_STRING([--disable-tbb],
                               [build without threading support via Threading Building Blocks]),
                [AS_CASE("${enableval}",
                         [yes], [enabletbb=yes],
                         [no],  [enabletbb=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-tbb)])],
                [enabletbb=$enableoptional])

  AS_IF([test "x$enabletbb" = "xyes"],
        [
          AC_ARG_WITH(tbb,
                      AS_HELP_STRING([--with-tbb=PATH],[Specify the path where Threading Building Blocks is installed]),
                      withtbb=$withval,
                      withtbb=$TBB_DIR)

          AC_ARG_WITH(tbb-lib,
                      AS_HELP_STRING([--with-tbb-lib=PATH],[Specify the path to Threading Building Blocks libraries]),
                      withtbblib=$withval,
                      withtbblib=$TBB_LIB_PATH)

          AS_IF([test "$withtbb" != no],
                [
                  AS_IF([test "x$withtbb" = "x"], [withtbb=/usr])

                  dnl Expand a leading ~ to $HOME so that --with-tbb=~/path works.
                  dnl bash does not expand ~ after = in non-assignment arguments.
                  case $withtbb in
                    "~/"*) withtbb=$HOME${withtbb#"~"} ;;
                  esac

                  dnl Detect which TBB era we have.
                  dnl oneTBB (>= 2021) ships tbb/version.h; legacy Intel TBB uses tbb/tbb_stddef.h.
                  dnl tbb/version.h did not exist in the legacy releases, so its presence is a
                  dnl reliable signal for oneTBB.
                  tbb_is_onetbb=no
                  AS_IF([test -r $withtbb/include/tbb/version.h],
                        [TBB_INCLUDE_PATH=$withtbb/include
                         tbb_is_onetbb=yes],
                        [AS_IF([test -r $withtbb/include/tbb/tbb_stddef.h],
                               [TBB_INCLUDE_PATH=$withtbb/include])])

                  AS_IF([test "x$withtbblib" != "x"], [TBB_LIBS=$withtbblib], [TBB_LIBS=$withtbb/lib])
                ])

          AS_IF([test "x$TBB_INCLUDE_PATH" != "x"],
                [
                  TBB_LIBRARY="-L$TBB_LIBS -ltbb -ltbbmalloc"
                  TBB_INCLUDE=-I$TBB_INCLUDE_PATH

                  dnl Add rpath flags to the link line.
                  AS_IF([test "x$RPATHFLAG" != "x" && test -d "$TBB_LIBS"], [TBB_LIBRARY="${RPATHFLAG}${TBB_LIBS} $TBB_LIBRARY"])

                  dnl Extract TBB_VERSION_MAJOR and TBB_VERSION_MINOR from the era-appropriate header.
                  dnl oneTBB: tbb/version.h   Legacy TBB: tbb/tbb_stddef.h
                  AS_IF([test "x$tbb_is_onetbb" = "xyes"],
                        [tbbverfile=$TBB_INCLUDE_PATH/tbb/version.h],
                        [tbbverfile=$TBB_INCLUDE_PATH/tbb/tbb_stddef.h])
                  tbbmajor=`grep "define TBB_VERSION_MAJOR" $tbbverfile | sed -e "s/#define TBB_VERSION_MAJOR[ ]*//g"`
                  tbbminor=`grep "define TBB_VERSION_MINOR" $tbbverfile | sed -e "s/#define TBB_VERSION_MINOR[ ]*//g"`
                ],
                [enabletbb=no])

          dnl If TBB is still enabled, verify with a compile test.
          dnl tbb::blocked_range is present in all TBB versions and works as a universal probe.
          AS_IF([test "x$enabletbb" != "xno"],
                [
                  AC_MSG_CHECKING(for TBB support)
                  AC_LANG_PUSH([C++])

                  saveCXXFLAGS="$CXXFLAGS"
                  CXXFLAGS="$saveCXXFLAGS $TBB_INCLUDE"

                  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
                      @%:@include <tbb/blocked_range.h>
                    ],
                    [
                      tbb::blocked_range<int> r(0, 1);
                      (void)r.size();
                    ])],
                    [
                      AC_MSG_RESULT(yes)
                      enabletbb=yes
                    ],
                    [
                      AC_MSG_RESULT(no)
                      enabletbb=no
                    ])

                  CXXFLAGS=$saveCXXFLAGS
                  AC_LANG_POP([C++])
                ])

          dnl If TBB is still enabled at this point, set all the necessary defines and print
          dnl a success message.
          AS_IF([test "x$enabletbb" != "xno"],
                [
                  AC_DEFINE_UNQUOTED(DETECTED_TBB_VERSION_MAJOR, [$tbbmajor], [TBB's major version number, as detected by tbb.m4])
                  AC_DEFINE_UNQUOTED(DETECTED_TBB_VERSION_MINOR, [$tbbminor], [TBB's minor version number, as detected by tbb.m4])

                  AS_IF([test "x$tbb_is_onetbb" = "xyes"],
                        [AC_DEFINE(HAVE_ONETBB, 1, [Defined if the installed TBB is oneTBB (>= 2021)])])

                  AC_SUBST(TBB_LIBRARY)
                  AC_SUBST(TBB_INCLUDE)
                  AC_DEFINE(USING_THREADS, 1, [Flag indicating whether the library shall be compiled to use any particular thread API.])
                  AC_DEFINE(HAVE_TBB_API, 1, [Flag indicating whether the library shall be compiled to use the Threading Building Blocks])
                  AC_MSG_RESULT(<<< Configuring library with Intel TBB threading support >>>)

                  dnl look for thread-local storage
                  AX_TLS
                ])
        ])
])
