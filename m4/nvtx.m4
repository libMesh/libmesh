dnl -------------------------------------------------------------
dnl NVIDIA Tools Extension Library (NVTX)
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_NVTX],
[
  AC_ARG_ENABLE(nvtx,
                AS_HELP_STRING([--disable-nvtx],
                               [build without annotation support via NVTX]),
                [AS_CASE("${enableval}",
                         [yes], [enablenvtx=yes],
                         [no],  [enablenvtx=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-nvtx)])],
                [enablenvtx=$enableoptional])

  AS_IF([test "x$enablenvtx" = "xyes"],
        [
          AC_ARG_WITH(nvtx,
                      AS_HELP_STRING([--with-nvtx=PATH],[Specify the path where NVTX is installed]),
                      withnvtx=$withval,
                      withnvtx=$NVTX_DIR)

          AS_IF([test "$withnvtx" != no],
                [
                  AS_IF([test "x$withnvtx" = "x"], [withnvtx=/usr])
                  AC_CHECK_HEADER($withnvtx/include/nvtx3/nvtx3.hpp, NVTX_INCLUDE_PATH=$withnvtx/include)
                ])

          AS_IF([test -r $NVTX_INCLUDE_PATH/nvtx3/nvtx3.hpp],
                [
                  NVTX_INCLUDE=-I$NVTX_INCLUDE_PATH
                ],
                [enablenvtx=no])

          dnl If NVTX is still enabled at this point, make sure we can compile
          dnl a test code which uses nvtx3::mark
          AS_IF([test "x$enablenvtx" != "xno"],
                [
                  AC_MSG_CHECKING(for nvtx3::mark support)
                  AC_LANG_PUSH([C++])

                  dnl Add NVTX headers to CXXFLAGS, which will be used by AC_COMPILE_IFELSE.
                  saveCXXFLAGS="$CXXFLAGS"
                  CXXFLAGS="$saveCXXFLAGS $NVTX_INCLUDE"

                  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
                      @%:@include <nvtx3/nvtx3.hpp>
                    ],
                    [
                      nvtx3::mark("Hello world!");
                    ])],
                    [
                      AC_MSG_RESULT(yes)
                      enablenvtx=yes
                    ],
                    [
                      AC_MSG_RESULT(no)
                      enablenvtx=no
                    ])

                  dnl Restore original flags
                  CXXFLAGS=$saveCXXFLAGS

                  AC_LANG_POP([C++])
                ])

          dnl If NVTX is still enabled at this point, set the necessary define and print
          dnl a success message.
          AS_IF([test "x$enablenvtx" != "xno"],
                [
                  AC_SUBST(NVTX_INCLUDE)
                  AC_DEFINE(HAVE_NVTX_API, 1, [Flag indicating whether the library shall be compiled to use NVTX annotations])
                  AC_MSG_RESULT(<<< Configuring library with NVTX annotation support >>>)
                ])
        ])
])
