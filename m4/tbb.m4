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

      AC_SUBST(TBB_LIBRARY)
      AC_SUBST(TBB_INCLUDE)
      AC_DEFINE(USING_THREADS, 1,
                [Flag indicating whether the library shall be compiled to use any particular thread API.])
      AC_DEFINE(HAVE_TBB_API, 1,
                [Flag indicating whether the library shall be compiled to use the Threading Building Blocks])
      AC_MSG_RESULT(<<< Configuring library with Intel TBB threading support >>>)

      dnl look for thread-local storage
      AX_TLS
    else
      enabletbb=no
    fi
  fi
])
