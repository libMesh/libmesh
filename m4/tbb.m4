dnl -------------------------------------------------------------
dnl Threading Building Blocks
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TBB],
[
  AC_ARG_ENABLE(tbb,
                AS_HELP_STRING([--disable-tbb],
                               [build without threading support via Threading Building Blocks]),
                [case "${enableval}" in
                  yes)  enabletbb=yes ;;
                  no)  enabletbb=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-tbb) ;;
                esac],
                [enabletbb=$enableoptional])

  if (test $enabletbb = yes); then

    AC_ARG_WITH(tbb,
                AS_HELP_STRING([--with-tbb=PATH],[Specify the path where Threading Building Blocks is installed]),
                withtbb=$withval,
                withtbb=$TBB_DIR)

    AC_ARG_WITH(tbb-lib,
                AS_HELP_STRING([--with-tbb-lib=PATH],[Specify the path to Threading Building Blocks libraries]),
                withtbblib=$withval,
                withtbblib=$TBB_LIB_PATH)

    if test "$withtbb" != no ; then
      if test "x$withtbb" = x ; then
        withtbb=/usr
      fi
      AC_CHECK_HEADER($withtbb/include/tbb/task_scheduler_init.h,
                      TBB_INCLUDE_PATH=$withtbb/include)
      if test "x$withtbblib" != "x" ; then
        TBB_LIBS=$withtbblib
      else
        TBB_LIBS=$withtbb/lib
      fi
    fi

    if (test -r $TBB_INCLUDE_PATH/tbb/task_scheduler_init.h) ; then
      TBB_LIBRARY="-L$TBB_LIBS -ltbb -ltbbmalloc"
      TBB_INCLUDE=-I$TBB_INCLUDE_PATH

      dnl Add rpath flags to the link line.
      if (test "x$RPATHFLAG" != "x" -a -d $TBB_LIBS); then
        TBB_LIBRARY="${RPATHFLAG}${TBB_LIBS} $TBB_LIBRARY"
      fi

      dnl Extract TBB_VERSION_MAJOR and TBB_VERSION_MINOR from
      dnl tbb_stddef.h.  This will allow us to set up a
      dnl TBB_VERSION_LESS_THAN macro.
      tbbmajor=`grep "define TBB_VERSION_MAJOR" $TBB_INCLUDE_PATH/tbb/tbb_stddef.h | sed -e "s/#define TBB_VERSION_MAJOR[ ]*//g"`
      tbbminor=`grep "define TBB_VERSION_MINOR" $TBB_INCLUDE_PATH/tbb/tbb_stddef.h | sed -e "s/#define TBB_VERSION_MINOR[ ]*//g"`
    else
      enabletbb=no
    fi

    # If TBB is still enabled at this point, make sure we can compile
    # a test code which uses tbb::tbb_thread
    if test "$enabletbb" != no ; then

      AC_MSG_CHECKING(for tbb::tbb_thread support)
      AC_LANG_PUSH([C++])

      # Add TBB headers to CXXFLAGS, which will be used by AC_COMPILE_IFELSE.
      saveCXXFLAGS="$CXXFLAGS"
      CXXFLAGS="$saveCXXFLAGS $TBB_INCLUDE"

      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
          @%:@include <tbb/tbb_thread.h>
        ],
        [
          tbb::tbb_thread t;
          t.join();
        ])],
        [
          AC_MSG_RESULT(yes)
          enabletbb=yes
        ],
        [
          AC_MSG_RESULT(no)
          enabletbb=no
        ])

      # Restore original flags
      CXXFLAGS=$saveCXXFLAGS

      AC_LANG_POP([C++])
    fi


    # If TBB is still enabled at this point, set all the necessary defines and print
    # a success message.
    if test "$enabletbb" != no ; then
      AC_DEFINE_UNQUOTED(DETECTED_TBB_VERSION_MAJOR, [$tbbmajor],
        [TBB's major version number, as detected by tbb.m4])

      AC_DEFINE_UNQUOTED(DETECTED_TBB_VERSION_MINOR, [$tbbminor],
        [TBB's minor version number, as detected by tbb.m4])


      AC_SUBST(TBB_LIBRARY)
      AC_SUBST(TBB_INCLUDE)
      AC_DEFINE(USING_THREADS, 1,
                [Flag indicating whether the library shall be compiled to use any particular thread API.])
      AC_DEFINE(HAVE_TBB_API, 1,
                [Flag indicating whether the library shall be compiled to use the Threading Building Blocks])
      AC_MSG_RESULT(<<< Configuring library with Intel TBB threading support >>>)

      dnl look for thread-local storage
      AX_TLS
    fi
  fi
])
