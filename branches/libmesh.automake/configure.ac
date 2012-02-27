# --------------------------------------------------------------
# $Id$
# --------------------------------------------------------------



AC_INIT([libmesh], [0.7.3rc1], [libmesh-users@lists.sourceforge.net], [libmesh], [http://libmesh.sourceforge.net])

# by default CFLAGS and CXXFLAGS are set to '-g -O2' on systems that support them.
# this causes a problem if we then declare a different optimization level. So
# default them to empty, regardless of what the system supports.
: ${CFLAGS=""}
: ${CXXFLAGS=""}

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----------- Configuring libMesh -------------)
AC_MSG_RESULT(---------------------------------------------)

AC_REVISION($Revision$)

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

# --------------------------------------------------------------
# Compute the canonical target system type variable "target", and its
# three individual parts: "target_cpu", "target_vendor", and "target_os"
# --------------------------------------------------------------
AC_CANONICAL_TARGET

AC_CONFIG_FILES([Make.common
                 Makefile
		 include/Makefile
		 examples/Makefile
		 libmesh.pc
		 ])
		 
AC_CONFIG_FILES(doc/Doxyfile)
AC_CONFIG_FILES(contrib/bin/libmesh-config,         [chmod +x contrib/bin/libmesh-config])
AC_CONFIG_FILES(contrib/bin/xda2mgf,                [chmod +x contrib/bin/xda2mgf])
AC_CONFIG_FILES(contrib/bin/create_libmesh_release, [chmod +x contrib/bin/create_libmesh_release])
AC_CONFIG_FILES(doc/dot,                            [chmod +x doc/dot])
# Create a temporary file with e.g. HAVE_MPI
AC_CONFIG_HEADER(include/base/libmesh_config.h.tmp:include/base/libmesh_config.h.in)
# Prefixes everything in the temporary file with LIBMESH_, e.g. LIBMESH_HAVE_MPI
AX_PREFIX_CONFIG_H(include/base/libmesh_config.h,LIBMESH,include/base/libmesh_config.h.tmp)

AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar color-tests parallel-tests])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES(yes)])




# --------------------------------------------------------------
# Print configuration info to the header file
# --------------------------------------------------------------
AC_DEFINE_UNQUOTED(CONFIGURE_INFO, ["$0 run on `date`"], [Configuration information.])	

# --------------------------------------------------------------
# Do compilation tests using the compiler, preprocessor, and file extensions
# of C++.  This uses the .C extension for test programs.
# --------------------------------------------------------------
AC_LANG([C++])


# --------------------------------------------------------------
# Set our preferred compilers
# --------------------------------------------------------------
LIBMESH_SET_COMPILERS

# --------------------------------------------------------------
# See aclocal.m4 for the definition of this as well.  It uses the
# output of DETERMINE_CXX_BRAND to set special flags for idiosyncratic
# compilers.  For example, newer versions of GCC have additional C++ debugging
# flags like _GLIBCXX_DEBUG_PEDANTIC which were not present in older versions.
# --------------------------------------------------------------
LIBMESH_SET_CXX_FLAGS
LIBMESH_SET_METHOD

LT_INIT([disable-static])

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
# Shared libraries - enabled by default
# --------------------------------------------------------------
AC_ARG_ENABLE(shared,
              AC_HELP_STRING([--enable-shared],
                             [build shared libraries]),
              enableshared=$enableval,
              enableshared=yes)
AC_SUBST(enableshared)	
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
# The AC_SUBST(foo) command creates an output variable from a shell variable.
# This command also causes AC_OUTPUT to replace all instances of @VAR@ with
# the value of the shell variable VAR in the output files (e.g. one or more Makefiles).
# --------------------------------------------------------------
AC_SUBST(host)
AC_SUBST(GXX_VERSION)
AC_SUBST(CXXFLAGS_OPT)
AC_SUBST(CXXFLAGS_DBG)
AC_SUBST(CXXFLAGS_DVL)
AC_SUBST(CFLAGS_OPT)
AC_SUBST(CFLAGS_DBG)
AC_SUBST(CFLAGS_DVL)
AC_SUBST(LDFLAGS)
AC_SUBST(NODEPRECATEDFLAG)
AC_SUBST(CXXSHAREDFLAG)
AC_SUBST(RPATHFLAG)
AC_SUBST(PROFILING_FLAGS)
AC_SUBST(OPROFILE_FLAGS)
AC_SUBST(FLIBS)
AC_SUBST(METHOD)

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
# Portable symbolic links
# --------------------------------------------------------------
AC_PROG_LN_S

# --------------------------------------------------------------
# Check for headers
# --------------------------------------------------------------
AC_CHECK_HEADERS(getopt.h)
AC_CHECK_HEADERS(csignal)
AC_CHECK_HEADERS(sys/resource.h)
AC_CXX_HAVE_LOCALE
AC_CXX_HAVE_SSTREAM

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


# --------------------------------------------------------------
# This is a special test for the mpCC (IBM) compiler.
# mpCC fails the AC_LANG tests. Set the compiler to 
# mpCC if it is there instead of xlC, but only do this
# _after_ configure is finished using CXX!
# --------------------------------------------------------------
if test "$CXX" = xlC_r; then
  #           Var   | look for | name if found | name if not | where
  AC_CHECK_PROG(MPCXX,  mpCC_r,      mpCC_r,           $CXX,         $PATH)
  CXX=$MPCXX
  if test "$CXX" = mpCC_r; then
    AC_MSG_RESULT(<<< Using $CXX instead of xlC_r to compile C++ >>>)
  fi
fi
if test "$CXX" = xlC; then
  #           Var   | look for | name if found | name if not | where
  AC_CHECK_PROG(MPCXX,  mpCC,      mpCC,           $CXX,         $PATH)
  CXX=$MPCXX
  if test "$CXX" = mpCC; then
    AC_MSG_RESULT(<<< Using $CXX instead of xlC to compile C++ >>>)
  fi
fi
# --------------------------------------------------------------


# --------------------------------------------------------------
# Get the source tree directory
# --------------------------------------------------------------
PWD=`pwd`
AC_SUBST(PWD)



# --------------------------------------------------------------
# Look for a user or system-provided boost.  If there is not
# one available then use the minimal ./contrib/boost provided.
# --------------------------------------------------------------
AX_BOOST_BASE([], AC_MSG_RESULT(<<< Using external boost installation >>>))
  # If that did not work, try using our builtin boost
  if test "x$BOOST_CPPFLAGS" = "x" ; then
    BOOST_ROOT="$PWD/contrib"
    AX_BOOST_BASE([], AC_MSG_RESULT(<<< Using libmesh-provided boost in ./contrib >>>))
  fi
# --------------------------------------------------------------



# --------------------------------------------------------------
# configure parts of the libmesh core
# --------------------------------------------------------------
LIBMESH_CONFIGURE_FEATURES



# --------------------------------------------------------------
# configure optional, external packages and libraries
# --------------------------------------------------------------
LIBMESH_CONFIGURE_OPTIONAL_PACKAGES



# -------------------------------------------------------------
# Perl - look for perl.  Perl is used to generate dependencies
# -------------------------------------------------------------
AC_PATH_PROG(PERL, perl)
AC_SUBST(PERL)
# -------------------------------------------------------------

# Query configuration environment
AX_SUMMARIZE_ENV


AC_OUTPUT

# Final summary
AX_SUMMARIZE_CONFIG



AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(--------- Done Configuring libMesh ----------)
AC_MSG_RESULT(---------------------------------------------)
# rm -f -r conftest* config.cache include/base/libmesh_config.h.tmp
