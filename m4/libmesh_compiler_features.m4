# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_COMPILER_FEATURES],
[
AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(------- Configuring compiler features -------)
AC_MSG_RESULT(---------------------------------------------)

# Test using the devel mode flags, which should include any
# user-specified flags
libmesh_compiler_features_save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS_DEVEL"
AC_LANG_SAVE
AC_LANG_CPLUSPLUS

# --------------------------------------------------------------
# Real precision - double by default
# --------------------------------------------------------------
ACX_CHOOSE_PRECISION

# --------------------------------------------------------------
# Determine support for the C99 "restrict" keyword in the C++
# compiler.  This should place a pound-define statement in the
# config file along the lines of:
#
# #define _libmesh_restrict __restrict
#
# The restrict keyword allows the compiler to make
# certain optimizations which require non-overlapping pointer
# addresses, i.e. in a call to:
#
# foo(int* restrict a, int* restrict b)
#
# the compiler may assume that the memory addresses pointed to
# by a and b do not overlap, and therefore certain low-level
# optimizations can be made.  From the Autoconf documentation:
#
# If the C compiler recognizes a variant spelling for the restrict
# keyword (__restrict, __restrict__, or _Restrict), then define restrict
# to that; this is more likely to do the right thing with compilers that
# support language variants where plain restrict is not a
# keyword. Otherwise, if the C compiler recognizes the restrict keyword,
# don't do anything. Otherwise, define restrict to be empty. Thus,
# programs may simply use restrict as if every C compiler supported it;
# for those that do not, the makefile or configuration header defines it
# away.
#
# Although support in C++ for the restrict keyword is not required,
# several C++ compilers do accept the keyword. This macro works for
# them, too.
#
# This macro caches 'no' in the ac_cv_c_restrict variable if restrict is
# not supported, and a supported spelling otherwise.
AC_C_RESTRICT


# --------------------------------------------------------------
# getpwuid - enabled by default
# Some systems, for example Crays, actually have getpwuid on the head-node
# but (if I understand correctly) a dynamically-linked glibc is not available
# on the backend, which is running a reduced operating system like Compute
# Node Linux.  Thus functions like getpwuid cannot be called.  This makes
# automatically testing for the existence of getpwuid on the login node
# difficult.  If you know that you are on such a system, configure with
# --disable-getpwuid.
# --------------------------------------------------------------
AC_ARG_ENABLE(getpwuid,
              AC_HELP_STRING([--enable-getpwuid],
                             [allow calls to getpwuid]),
              enablegetpwuid=$enableval,
              enablegetpwuid=yes)

if test "$enablegetpwuid" != no ; then
  AC_DEFINE(HAVE_GETPWUID, 1,
           [Flag indicating if the library should be built with calls to getpwuid()])
  AC_MSG_RESULT(<<< Configuring library with getpwuid >>>)
fi
# --------------------------------------------------------------



# --------------------------------------------------------------
# C++ exceptions - enabled by default
# --------------------------------------------------------------
AC_ARG_ENABLE(exceptions,
              AC_HELP_STRING([--enable-exceptions],
                             [throw exceptions on unexpected errors]),
              enableexceptions=$enableval,
              enableexceptions=yes)

if test "$enableexceptions" != no ; then
  AC_DEFINE(ENABLE_EXCEPTIONS, 1,
           [Flag indicating if the library should be built to throw C++ exceptions on unexpected errors])
  AC_MSG_RESULT(<<< Configuring library with exception throwing support >>>)
fi
# --------------------------------------------------------------



# --------------------------------------------------------------
# Check for important type sizes
# --------------------------------------------------------------
AC_CHECK_SIZEOF(short int)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(long int)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(void *)
# Check the size of a function pointer.  This is the same
# as a void* on most systems which matter (POSIX).
# It turns out that AC_CHECK_SIZEOF can't do this
# automatically... so we use a kind of hack with the 3rd
# argument to do it.
AC_CHECK_SIZEOF([function_pointer], [], [typedef void (*function_pointer)();])
# AC_CHECK_SIZEOF([void(*)(void)]) <-- Does not work!

# --------------------------------------------------------------
# Check for Run Time Type Identification support
# --------------------------------------------------------------
AC_CXX_RTTI



# --------------------------------------------------------------
# Check for headers
# --------------------------------------------------------------
AC_CHECK_HEADERS(getopt.h)
AC_CHECK_HEADERS(csignal)
AC_CHECK_HEADERS(sys/resource.h)
AC_CXX_HAVE_LOCALE
AC_CXX_HAVE_SSTREAM

# Check for GNU libc feenableexcept/fedisableexcept.
AC_CHECK_HEADERS(fenv.h)
AC_CHECK_HEADERS(xmmintrin.h)

# See if we can actually compile code using the fe{enable,disable}except functions
AC_HAVE_FEEXCEPT

AC_ARG_ENABLE(unordered-containers,
              AC_HELP_STRING([--enable-unordered-containers],
                             [Use unordered_map/unordered_set if available]),
              enableunorderedcontainers=$enableval,
              enableunorderedcontainers=yes)

  if test "$enableunorderedcontainers" != no ; then
    # The following routines, defined in unordered.m4, check to see if the compiler can compile programs using
    # various quasi-standard hash containers.
    ACX_BEST_UNORDERED_MULTIMAP
    ACX_BEST_UNORDERED_MAP
    ACX_BEST_UNORDERED_SET
  else
    ACX_STD_MAP
    ACX_STD_MULTIMAP
    ACX_STD_SET
  fi

AC_CHECK_HEADERS(dlfcn.h)
AX_CXX_GCC_ABI_DEMANGLE
AX_CXX_GLIBC_BACKTRACE



# -------------------------------------------------------------
# OpenMP Support  -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(openmp,
             AC_HELP_STRING([--enable-openmp],
                            [Build with OpenMP Support]),
             enableopenmp=$enableval,
             enableopenmp=yes)
if (test "$enableopenmp" != no) ; then
   AX_OPENMP([],[enableopenmp=no])
   #The above call only sets the flag for C++
   if (test "x$OPENMP_CXXFLAGS" != x) ; then
     AC_MSG_RESULT(<<< Configuring library with OpenMP support >>>)
     OPENMP_CFLAGS=$OPENMP_CXXFLAGS
     OPENMP_FFLAGS=$OPENMP_CXXFLAGS
     CXXFLAGS_OPT="$CXXFLAGS_OPT $OPENMP_CXXFLAGS"
     CXXFLAGS_DBG="$CXXFLAGS_DBG $OPENMP_CXXFLAGS"
     CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL $OPENMP_CXXFLAGS"
     CXXFLAGS_PROF="$CXXFLAGS_PROF $OPENMP_CXXFLAGS"
     CXXFLAGS_OPROF="$CXXFLAGS_OPROF $OPENMP_CXXFLAGS"
     CFLAGS_OPT="$CFLAGS_OPT $OPENMP_CFLAGS"
     CFLAGS_DBG="$CFLAGS_DBG $OPENMP_CFLAGS"
     CFLAGS_DEVEL="$CFLAGS_DEVEL $OPENMP_CFLAGS"
     CFLAGS_PROF="$CFLAGS_PROF $OPENMP_CFLAGS"
     CFLAGS_OPROF="$CFLAGS_OPROF $OPENMP_CFLAGS"
     FFLAGS="$FFLAGS $OPENMP_FFLAGS"
   fi
   AC_SUBST(OPENMP_CXXFLAGS)
   AC_SUBST(OPENMP_CFLAGS)
   AC_SUBST(OPENMP_FFLAGS)
fi
# -------------------------------------------------------------

# Restore original CXXFLAGS for now
AC_LANG_RESTORE
CXXFLAGS="$libmesh_compiler_features_save_CXXFLAGS"

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----- Done configuring compiler features ----)
AC_MSG_RESULT(---------------------------------------------)
])
