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


LIBMESH_COMPILER_FEATURES




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
# Portable symbolic links
# --------------------------------------------------------------
AC_PROG_LN_S


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
LIBMESH_CORE_FEATURES



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
