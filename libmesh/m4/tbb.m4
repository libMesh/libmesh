dnl -------------------------------------------------------------
dnl Threading Building Blocks
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TBB],
[
  AC_ARG_WITH(tbb,
              AC_HELP_STRING([--with-tbb=PATH],[Specify the path where Threading Building Blocks is installed]),
              withtbb=$withval,
              withtbb=$TBB_DIR)

  AC_ARG_WITH(tbb-lib,
              AC_HELP_STRING([--with-tbb-lib=PATH],[Specify the path to Threading Building Blocks libraries]),
              withtbblib=$withval,
              withtbblib=$TBB_LIB_PATH)

  if test "$withtbb" != no ; then
    if test "x$withtbb" = x ; then
      withtbb=/usr
    fi

    AC_LANG_PUSH([C++])
    OLD_CPPFLAGS=$CPPFLAGS
    AC_CHECK_HEADER(tbb/task_scheduler_init.h,
      [TBB_INCLUDE=''
       withtbb=builtin],
      [CPPFLAGS="-I$withtbb/include $CPPFLAGS"
       AC_CHECK_HEADER(tbb/task_scheduler_init.h,
                       TBB_INCLUDE='-I$withtbb/include',
                       withtbb=no)
       CPPFLAGS=$OLD_CPPFLAGS])
    AC_LANG_POP([C++])
  fi

  if test "$withtbb" != no ; then
    if test "x$withtbblib" != "x" ; then
      TBB_LIBRARY="-L$withtbblib -ltbb -ltbbmalloc -lpthread"
    else	
      if test "$withtbb" != "builtin" ; then
        TBB_LIBRARY="-L$withtbb/lib -ltbb -ltbbmalloc -lpthread"
      else
        TBB_LIBRARY="-ltbb -ltbbmalloc -lpthread"
      fi
    fi

    AC_SUBST(TBB_LIBRARY)
    AC_SUBST(TBB_INCLUDE)
    AC_DEFINE(HAVE_TBB_API, 1,
              [Flag indicating whether the library shall be compiled to use the Threading Building Blocks])
    AC_MSG_RESULT(<<< Configuring library with Intel TBB threading support >>>)
  fi
])
