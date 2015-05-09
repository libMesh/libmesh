# ----------------------------------------------------------------
# The nlopt package by Steven G. Johnson provides a common interface
# for a number of different free optimization routines.  Libmesh wraps
# the nlopt interface in an OptimizationSolver derived object.
# http://ab-initio.mit.edu/wiki/index.php/NLopt
# ----------------------------------------------------------------

AC_DEFUN([CONFIGURE_NLOPT],
[
  AC_ARG_ENABLE(nlopt,
                AS_HELP_STRING([--disable-nlopt],
                               [build without NLOPT support]),
                [case "${enableval}" in
                  yes)  enablenlopt=yes ;;
                  no)  enablenlopt=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-nlopt) ;;
                 esac],
                [enablenlopt=$enableoptional])


  if (test $enablenlopt = yes); then
    # User-specific include path
    AC_ARG_WITH(nlopt-include,
                AS_HELP_STRING([--with-nlopt-include=PATH],[Specify the path for NLOPT header files]),
                withnloptinc=$withval,
                withnloptinc=no)

    # User-specific library path
    AC_ARG_WITH(nlopt-lib,
                AS_HELP_STRING([--with-nlopt-lib=PATH],[Specify the path for NLOPT libs]),
                withnloptlib=$withval,
                withnloptlib=no)

    # Use NLOPT_DIR/include if it exists.
    if (test $withnloptinc != no); then
      NLOPT_INC="$withnloptinc"
    elif test "x$NLOPT_DIR" != x -a -f $NLOPT_DIR/include/nlopt.h; then
      NLOPT_INC="$NLOPT_DIR/include"
    else
      NLOPT_INC=""
    fi

    # Use NLOPT_DIR/lib if it exists.
    if (test $withnloptlib != no); then
      NLOPT_LIB="$withnloptlib"
    elif test "x$NLOPT_DIR" != x; then
      NLOPT_LIB="$NLOPT_DIR/lib"
    else
      NLOPT_LIB=""
    fi

    # Initialize Makefile/config.h substitution variables
    NLOPT_INCLUDE=""
    NLOPT_LIBRARY=""

    if (test $enablenlopt = yes); then
      NLOPT_INCLUDE="-I$NLOPT_INC"
      NLOPT_LIBRARY="-L$NLOPT_LIB -lnlopt"

      # Try to compile and link a trivial nlopt test code.
      AC_MSG_CHECKING(for valid nlopt installation)
      AC_LANG_PUSH([C])

      # Save the original CFLAGS and LIBS contents
      saveCFLAGS="$CFLAGS"
      saveLIBS="$LIBS"

      # Append nlopt include paths to the CFLAGS variables
      CFLAGS="$saveCFLAGS $NLOPT_INCLUDE"
      LIBS="$saveLIBS $NLOPT_LIBRARY -lm"

      AC_LINK_IFELSE(
           [
             AC_LANG_PROGRAM([[
@%:@include <math.h>
@%:@include <nlopt.h>
@%:@include <stdio.h>

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
if (grad) {
  grad[0] = 0.0;
  grad[1] = 0.5 / sqrt(x[1]);
 }
 return sqrt(x[1]);
}

typedef struct
{
  double a, b;
} my_constraint_data;

double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
  my_constraint_data *d = (my_constraint_data *) data;
  double a = d->a, b = d->b;
  if (grad) {
    grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
    grad[1] = -1.0;
  }
  return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}
                              ]],
                             [[
double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
nlopt_opt opt;

opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
nlopt_set_lower_bounds(opt, lb);
nlopt_set_min_objective(opt, myfunc, NULL);

my_constraint_data data[2] = { {2,0}, {-1,1} };
nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);
nlopt_set_xtol_rel(opt, 1e-4);
double x[2] = { 1.234, 5.678 };  /* some initial guess */
double minf; /* the minimum objective value, upon return */

if (nlopt_optimize(opt, x, &minf) < 0) {
  printf("nlopt failed!\n");
}
else {
  printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
}
nlopt_destroy(opt);
return 0;
                             ]])
           ],
           [
             enablenlopt=yes
             AC_MSG_RESULT(yes)
           ],
           [
             enablenlopt=no
             AC_MSG_RESULT(no)
           ])

      # Return CFLAGS and LIBS to their original states.
      CFLAGS="$saveCFLAGS"
      LIBS="$saveLIBS"
      AC_LANG_POP([C])

      # If linking a test program succeeded, continue.
      if (test x$enablenlopt = xyes); then

        # add the NLOPT_LIB dir to the linker run path if it is a directory
        if (test "x$RPATHFLAG" != "x" -a -d $NLOPT_LIB); then
          NLOPT_LIBRARY="${RPATHFLAG}${NLOPT_LIB} $NLOPT_LIBRARY"
        fi

        AC_DEFINE(HAVE_NLOPT, 1, [Flag indicating whether the library will be compiled with NLOPT support])
        AC_MSG_RESULT(<<< Configuring library with NLOPT support >>>)
      fi
    fi
  fi

  # Substitute the substitution variables
  AC_SUBST(NLOPT_INCLUDE)
  AC_SUBST(NLOPT_LIBRARY)
])
