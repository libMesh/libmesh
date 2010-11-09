dnl -------------------------------------------------------------
dnl $Id: backtrace.m4 3808 2010-05-03 18:30:34Z benkirk $
dnl -------------------------------------------------------------

dnl ----------------------------------------------------------------------------
dnl Accept options for single-precision or triple-precision scalars
dnl ----------------------------------------------------------------------------

AC_DEFUN([ACX_CHOOSE_PRECISION],
[
AC_ARG_ENABLE(singleprecision,
              AC_HELP_STRING([--enable-singleprecision],
                             [Use single-precision scalars]),
              enablesingleprecision=$enableval,
              enablesingleprecision=no)

AC_ARG_ENABLE(tripleprecision,
              AC_HELP_STRING([--enable-tripleprecision],
                             [Use triple-precision scalars]),
              enabletripleprecision=$enableval,
              enabletripleprecision=no)

if test "$enablesingleprecision" != no ; then
  if test "$enabletripleprecision" != no ; then
    AC_MSG_ERROR(<<< Cannot simultaneously default to single and triple precision >>>)
  else
    AC_DEFINE(SINGLE_PRECISION, 1,
              [Flag indicating if single-precision (float) should be used for most floating-point calculations])
    AC_DEFINE(SCALAR_TYPE, float,
              [Data type to be used for most floating-point calculations])
    AC_MSG_RESULT(<<< Default floating point is single precision (float) >>>)
  fi
elif test "$enabletripleprecision" != no ; then
  AC_DEFINE(TRIPLE_PRECISION, 1,
            [Flag indicating if triple-precision (long double) should be used for most floating-point calculations])
  AC_DEFINE(SCALAR_TYPE, [long double],
            [Data type to be used for most floating-point calculations])
  AC_MSG_RESULT(<<< Default floating point is triple precision (long double) >>>)
else
  AC_DEFINE(DOUBLE_PRECISION, 1,
            [Flag indicating if double-precision (double) should be used for most floating-point calculations])
  AC_DEFINE(SCALAR_TYPE, double,
            [Data type to be used for most floating-point calculations])
  AC_MSG_RESULT(<<< Default floating point is double precision (double) >>>)
fi
])
