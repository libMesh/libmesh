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
enablegetpwuid_default=yes

# We can't use getpwuid if pwd.h is not found.
AC_CHECK_HEADERS(pwd.h, [have_pwd_h=yes], [have_pwd_h=no])
if test "$have_pwd_h" = no ; then
  enablegetpwuid_default=no
fi

# If the user configured with --enable-all-static, certain functions
# like getpwuid cannot be used. You may see warnings like the
# following from the linker:
# warning: Using 'getpwuid' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking
if test "$enableallstatic" = yes ; then
  enablegetpwuid_default=no
fi

# Let the user override the default (at their own risk).
AC_ARG_ENABLE(getpwuid,
              AS_HELP_STRING([--disable-getpwuid],
                             [do not make calls to getpwuid]),
              enablegetpwuid=$enableval,
              enablegetpwuid=$enablegetpwuid_default)

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
              AS_HELP_STRING([--disable-exceptions],
                             [exit rather than throw exceptions on unexpected errors]),
              enableexceptions=$enableval,
              enableexceptions=yes)

if test "$enableexceptions" != no ; then
  AC_DEFINE(ENABLE_EXCEPTIONS, 1,
           [Flag indicating if the library should be built to throw C++ exceptions on unexpected errors])
  AC_MSG_RESULT(<<< Configuring library with exception throwing support >>>)
fi
# --------------------------------------------------------------



# --------------------------------------------------------------
# __TIME__ __DATE__ stamps - enabled by default
# disabling preprocessor timestamps helps compiler caches such
# as ccache to work more effectively.
# --------------------------------------------------------------
AC_ARG_ENABLE(timestamps,
              AS_HELP_STRING([--disable-timestamps],
                             [do not add preprocessor timestamps to the library (helps ccache)]),
              enabletimestamps=$enableval,
              enabletimestamps=yes)

if test "$enabletimestamps" != no ; then
  AC_DEFINE(ENABLE_TIMESTAMPS, 1,
           [Flag indicating if the library should be built with compile time and date timestamps])
  AC_MSG_RESULT(<<< Configuring library with compile timestamps >>>)
fi
# --------------------------------------------------------------



# --------------------------------------------------------------
# Check for important type sizes
# --------------------------------------------------------------
AC_CHECK_SIZEOF(short int)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(size_t)
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

# Check if we have new-style signal handlers
AC_CHECK_DECLS([sigaction], [], [], [#include <signal.h>])

# Check for mkdir
AC_MSG_CHECKING([for mkdir with two arguments])
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([
        #include <sys/stat.h>
        #include <sys/types.h> ],[
        mkdir("test-dir", 0); ])],
    [AC_DEFINE(HAVE_MKDIR, [], "The make directory command")
     AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)])
AC_CHECK_HEADERS(direct.h)
AC_CHECK_DECLS([_mkdir], [], [], [
    #include <direct.h>
])

# Check for uname header.
AC_CHECK_HEADERS(sys/utsname.h)

AC_ARG_ENABLE(unordered-containers,
              AS_HELP_STRING([--disable-unordered-containers],
                             [Use map/set instead of unordered_map/unordered_set]),
              enableunorderedcontainers=$enableval,
              enableunorderedcontainers=yes)

  if test "$enableunorderedcontainers" != no ; then
    # The following routines, defined in unordered.m4, check to see if the compiler can compile programs using
    # various quasi-standard hash containers.
    ACX_BEST_UNORDERED_MULTIMAP
    ACX_BEST_UNORDERED_MAP
    ACX_BEST_UNORDERED_MULTISET
    ACX_BEST_UNORDERED_SET
  else
    ACX_STD_MAP
    ACX_STD_MULTIMAP
    ACX_STD_SET
    ACX_STD_MULTISET
  fi

# Determine which of std::hash, std::tr1::hash, or __gnu_cxx::hash is available
ACX_BEST_HASH

AX_CXX_DLOPEN

dnl Set preprocessor macro if the test code succeeded
if (test "$ac_cv_cxx_dlopen" = yes); then
  AC_DEFINE(HAVE_DLOPEN, 1, [define if the compiler supports dlopen/dlsym/dlclose])
fi

AX_CXX_GCC_ABI_DEMANGLE
AX_CXX_GLIBC_BACKTRACE



# See if this compiler has a broken errno_t in a very specific
# situation (this is not common).
CHECK_FOR_BROKEN_ERRNO_T


# Restore original CXXFLAGS for now
AC_LANG_RESTORE
CXXFLAGS="$libmesh_compiler_features_save_CXXFLAGS"

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----- Done configuring compiler features ----)
AC_MSG_RESULT(---------------------------------------------)
])
