dnl ----------------------------------------------------------------------------
dnl Accept options for single-precision or triple-precision scalars
dnl ----------------------------------------------------------------------------

AC_DEFUN([ACX_CHOOSE_PRECISION],
[
AC_ARG_ENABLE(singleprecision,
              AS_HELP_STRING([--enable-singleprecision],
                             [Use single-precision scalars]),
              enablesingleprecision=$enableval,
              enablesingleprecision=no)

AC_ARG_ENABLE(tripleprecision,
              AS_HELP_STRING([--enable-tripleprecision],
                             [Use triple-precision scalars]),
              enabletripleprecision=$enableval,
              enabletripleprecision=no)

AC_ARG_ENABLE(quadrupleprecision,
              AS_HELP_STRING([--enable-quadrupleprecision],
                             [Use quadruple-precision scalars]),
              enablequadrupleprecision=$enableval,
              enablequadrupleprecision=no)

libmesh_precision_LIBS=""

AS_IF([test "x$enablesingleprecision" != "xno"],
      [
        AS_IF([test "x$enabletripleprecision" != "xno"], [AC_MSG_ERROR(<<< Cannot simultaneously default to single and triple precision >>>)],
              [test "x$enablequadrupleprecision" != "xno"], [AC_MSG_ERROR(<<< Cannot simultaneously default to single and quadruple precision >>>)],
              [
                AC_DEFINE(DEFAULT_SINGLE_PRECISION, 1, [Flag indicating if single-precision (float) should be used for most floating-point calculations])
                AC_DEFINE(DEFAULT_SCALAR_TYPE, float, [Data type to be used for most floating-point calculations])
                AC_MSG_RESULT(<<< Default floating point is single precision (float) >>>)
              ])
      ],
      [test "x$enabletripleprecision" != "xno"],
      [
        AS_IF([test "x$enablequadrupleprecision" != "xno"],
              [AC_MSG_ERROR(<<< Cannot simultaneously default to single and quadruple precision >>>)],
              [
                AC_DEFINE(DEFAULT_TRIPLE_PRECISION, 1, [Flag indicating if triple-precision (long double) should be used for most floating-point calculations])
                AC_DEFINE(DEFAULT_SCALAR_TYPE, [long double], [Data type to be used for most floating-point calculations])
                AC_MSG_RESULT(<<< Default floating point is triple precision (long double) >>>)
              ])
      ],
      [test "x$enablequadrupleprecision" != "xno"],
      [
        AC_DEFINE(DEFAULT_QUADRUPLE_PRECISION, 1, [Flag indicating if quadruple-precision (boost::multiprecision::float128) should be used for most floating-point calculations])
        AC_DEFINE(DEFAULT_SCALAR_TYPE, [boost::multiprecision::float128], [Data type to be used for most floating-point calculations])
        AC_MSG_RESULT(<<< Default floating point is quadruple precision (boost::multiprecision::float128) >>>)

        AC_MSG_CHECKING(whether we can build a trivial quad precision program)

        saveLIBS="$LIBS"
        AS_IF([test "x$REAL_GXX" != "x"],
              [libmesh_precision_LIBS="-lquadmath"
               LIBS="$saveLIBS -lquadmath"
               AC_LINK_IFELSE([AC_LANG_SOURCE([[
                 @%:@include <quadmath.h>
                 int main(int argc, char **argv)
                 {
                   __float128 f = 1;
                   return isinfq(f);
                 }
               ]])],[
                 AC_MSG_RESULT(yes)
               ],[
                 AC_MSG_RESULT(no)
                 AC_MSG_ERROR([*** Quad precision specified, gcc detected, but quadmath not found.])
                 enablequadrupleprecision=no
               ])
              ],
              [test "x$is_intel_icc" != "x"],
              [libmesh_precision_LIBS=""
               AC_LINK_IFELSE([AC_LANG_SOURCE([[
                 @%:@include <mathimf.h>
                 int main(int argc, char **argv)
                 {
                   _Quad f = 0;
                   return int(__cosq(f));
                 }
               ]])],[
                 AC_MSG_RESULT(yes)
               ],[
                 AC_MSG_RESULT(no)
                 AC_MSG_ERROR([*** Quad precision specified, Intel detected, but mathimf not found.])
                 enablequadrupleprecision=no
               ])
              ])
        LIBS="$saveLIBS"
      ],
      [
        AC_DEFINE(DEFAULT_DOUBLE_PRECISION, 1, [Flag indicating if double-precision (double) should be used for most floating-point calculations])
        AC_DEFINE(DEFAULT_SCALAR_TYPE, double, [Data type to be used for most floating-point calculations])
        AC_MSG_RESULT(<<< Default floating point is double precision (double) >>>)
      ])

  AC_SUBST(libmesh_precision_LIBS)
])
