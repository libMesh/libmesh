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

AC_ARG_ENABLE(quadrupleprecision,
              AC_HELP_STRING([--enable-quadrupleprecision],
                             [Use quadruple-precision scalars]),
              enablequadrupleprecision=$enableval,
              enablequadrupleprecision=no)

if test "$enablesingleprecision" != no ; then
  if test "$enabletripleprecision" != no ; then
    AC_MSG_ERROR(<<< Cannot simultaneously default to single and triple precision >>>)
  elif test "$enablequadrupleprecision" != no ; then
    AC_MSG_ERROR(<<< Cannot simultaneously default to single and quadruple precision >>>)
  else
    AC_DEFINE(DEFAULT_SINGLE_PRECISION, 1,
              [Flag indicating if single-precision (float) should be used for most floating-point calculations])
    AC_DEFINE(DEFAULT_SCALAR_TYPE, float,
              [Data type to be used for most floating-point calculations])
    AC_MSG_RESULT(<<< Default floating point is single precision (float) >>>)
  fi
elif test "$enabletripleprecision" != no ; then
  if test "$enablequadrupleprecision" != no ; then
    AC_MSG_ERROR(<<< Cannot simultaneously default to single and quadruple precision >>>)
  else
    AC_DEFINE(DEFAULT_TRIPLE_PRECISION, 1,
              [Flag indicating if triple-precision (long double) should be used for most floating-point calculations])
    AC_DEFINE(DEFAULT_SCALAR_TYPE, [long double],
              [Data type to be used for most floating-point calculations])
    AC_MSG_RESULT(<<< Default floating point is triple precision (long double) >>>)
  fi
elif test "$enablequadrupleprecision" != no ; then
    AC_DEFINE(DEFAULT_QUADRUPLE_PRECISION, 1,
              [Flag indicating if quadruple-precision (__float128) should be used for most floating-point calculations])
    AC_DEFINE(DEFAULT_SCALAR_TYPE, [__float128],
              [Data type to be used for most floating-point calculations])
    AC_MSG_RESULT(<<< Default floating point is quadruple precision (__float128) >>>)
else
  AC_DEFINE(DEFAULT_DOUBLE_PRECISION, 1,
            [Flag indicating if double-precision (double) should be used for most floating-point calculations])
  AC_DEFINE(DEFAULT_SCALAR_TYPE, double,
            [Data type to be used for most floating-point calculations])
  AC_MSG_RESULT(<<< Default floating point is double precision (double) >>>)
fi
])
