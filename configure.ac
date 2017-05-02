# --------------------------------------------------------------
# configure.ac file for libmesh.  Used by Automake to generate
# aclocal.m4 and Makefile.in. See also: http://tinyurl.com/hatuyfn
# --------------------------------------------------------------
AC_INIT([libmesh], [1.3.0-pre], [libmesh-users@lists.sourceforge.net], [libmesh], [http://libmesh.github.io])
AC_PREREQ([2.69])

# by default CFLAGS and CXXFLAGS are set to '-g -O2' on systems that support them.
# this causes a problem if we then declare a different optimization level. So
# default them to empty, regardless of what the system supports.
: ${CFLAGS=""}
: ${CXXFLAGS=""}

# Infer the source directory as the path to the ./configure script
srcdir=`dirname $0`
top_srcdir=`dirname $0`


AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----------- Configuring libMesh -------------)
AC_MSG_RESULT(---------------------------------------------)

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

# --------------------------------------------------------------
# Compute the canonical target system type variable "target", and its
# three individual parts: "target_cpu", "target_vendor", and "target_os"
# --------------------------------------------------------------
AC_CANONICAL_TARGET

AC_CONFIG_FILES([Makefile
                 include/Makefile
                 include/libmesh/Makefile
                 contrib/Makefile
                 contrib/utils/Makefile
                 contrib/utils/Make.common
                 tests/Makefile
                 contrib/utils/libmesh.pc:contrib/utils/libmesh-opt.pc.in
                 contrib/utils/libmesh-opt.pc
                 contrib/utils/libmesh-dbg.pc
                 contrib/utils/libmesh-devel.pc
                 contrib/utils/libmesh-prof.pc
                 contrib/utils/libmesh-oprof.pc
                 doc/Doxyfile
                 doc/Makefile
                 doc/html/Makefile
                 ])

AC_CONFIG_FILES(contrib/bin/libmesh-config,          [chmod +x contrib/bin/libmesh-config])
#AC_CONFIG_FILES(contrib/bin/strip_dup_libs.pl,       [chmod +x contrib/bin/strip_dup_libs.pl])
#AC_CONFIG_FILES(contrib/bin/strip_dup_incl_paths.pl, [chmod +x contrib/bin/strip_dup_incl_paths.pl])
AC_CONFIG_FILES(contrib/bin/create_libmesh_release,  [chmod +x contrib/bin/create_libmesh_release])
AC_CONFIG_FILES(doc/dot,                             [chmod +x doc/dot])
AC_CONFIG_FILES(tests/run_unit_tests.sh,             [chmod +x tests/run_unit_tests.sh])

# Create a temporary file with e.g. HAVE_MPI
AC_CONFIG_HEADER(include/libmesh_config.h.tmp:include/libmesh_config.h.in)
# Prefixes everything in the temporary file with LIBMESH_, e.g. LIBMESH_HAVE_MPI
AX_PREFIX_CONFIG_H(include/libmesh_config.h,LIBMESH,include/libmesh_config.h.tmp)

# Print configuration info to the header file
AC_DEFINE_UNQUOTED(CONFIGURE_INFO, ["$0 $ac_configure_args"], [Configuration information.])

# require automake 1.11 - color tests
AM_INIT_AUTOMAKE([dist-xz dist-bzip2 tar-ustar color-tests serial-tests 1.12])

# use silent rules - automake 1.11
AM_SILENT_RULES(yes)

# We originally passed [enable] to AM_MAINTAINER_MODE, but this caused
# ordinary users to have to configure with --disable-maintainer-mode,
# which is easy to forget unless you are closely following the
# instructions, or are aware of the issue and have run the provided
# ./bootstrap script.
#
# We've now switched this to simply AM_MAINTAINER_MODE with no
# arguments, which is equivalent to AM_MAINTAINER_MODE([disable]).  In
# this configuration, developers who modify Makefile.am files (or
# other build system files) *must* configure with
# --enable-maintainer-mode in order for dependent files (Makefile.in,
# etc) to be properly rebuilt whenever changes to the build system are
# made.
#
# Note that this approach is not without its drawbacks and detractors:
# https://blogs.gnome.org/desrt/2011/09/08/am_maintainer_mode-is-not-cool
# https://www.gnu.org/software/automake/manual/html_node/maintainer_002dmode.html
# but the benefit is that casual and first-time users are not
# surprised by the default behavior of simply running ./configure.
#
# Note that we also do things a bit differently than most other
# projects regarding the number of auto-generated build files
# (Makefile.in) that are checked in, and it seems like the new
# approach makes more sense for the way our project works --
# rebuilding the build files are not part of a normal build for most
# users.
AM_MAINTAINER_MODE


# --------------------------------------------------------------
# Determine if we are doing a local or VPATH build. In the case
# of the latter we will need to go through some extra effort to
# create links for the example source files.
# --------------------------------------------------------------
AC_CHECK_FILE(src/base/libmesh.C, [libmesh_vpath_build=no], [libmesh_vpath_build=yes])
if (test $libmesh_vpath_build = yes); then
    AC_MSG_RESULT(<<< Configuring build directory for VPATH build >>>)
fi
AM_CONDITIONAL(LIBMESH_VPATH_BUILD, test x$libmesh_vpath_build = xyes)

# --------------------------------------------------------------
# Do compilation tests using the compiler, preprocessor, and file extensions
# of C++.  This uses the .C extension for test programs.
# --------------------------------------------------------------
AC_LANG([C++])

# -------------------------------------------------------------
# Perl - look for perl.  Perl is used to generate dependencies.  We
# also use it to (optionally) parse compiler version strings, so $PERL
# needs to be set before LIBMESH_SET_COMPILERS is called.
# -------------------------------------------------------------
AC_PATH_PROG(PERL, perl)
AC_SUBST(PERL)

# --------------------------------------------------------------
# Set our preferred compilers
# --------------------------------------------------------------
LIBMESH_SET_COMPILERS

# --------------------------------------------------------------
# Allow the user to explicitly turn off C++11 support by saying
# --disable-cxx11.  Otherwise, C++11 support is enabled by default
# although libmesh does not (currently) use any C++11 features.
# --------------------------------------------------------------
AC_ARG_ENABLE(cxx11,
              AC_HELP_STRING([--disable-cxx11],
                             [build without C++11 support]),
                             [case "${enableval}" in
                   yes) enablecxx11=yes ;;
                   no)  enablecxx11=no ;;
                   *)   AC_MSG_ERROR(bad value ${enableval} for --disable-cxx11) ;;
                   esac],
                   [enablecxx11=yes])


# --------------------------------------------------------------
# Setting --enable-cxx11-required causes an error to be emitted during
# configure if the compiler does not support a specific subset of
# C++11 features.  This is useful for app codes which require C++11,
# since it prevents situations where libmesh is built without C++11
# support (which may take a very long time), and then the app fails to
# compile, requiring you to redo everything.
# --------------------------------------------------------------
AC_ARG_ENABLE(cxx11-required,
              AC_HELP_STRING([--enable-cxx11-required],
                             [Error if compiler does not support C++11]),
                             [case "${enableval}" in
                   yes) cxx11required=yes ;;
                   no)  cxx11required=no ;;
                   *)   AC_MSG_ERROR(bad value ${enableval} for --enable-cxx11-required) ;;
                   esac],
                   [cxx11required=no])

# --------------------------------------------------------------
# Standard autoconf macro for determining the proper -std=c++11
# flag for the current compiler.  Adds the result to CXXFLAGS if
# one is found.  See ax_cxx_compile_stdcxx_11.m4 for details.
# --------------------------------------------------------------
if (test $cxx11required = yes); then
  AX_CXX_COMPILE_STDCXX_11([],[mandatory])
else
  AX_CXX_COMPILE_STDCXX_11([],[optional])
fi

# --------------------------------------------------------------
# See compiler.m4 for the definition of this function.  It uses the
# output of DETERMINE_CXX_BRAND to set special flags for idiosyncratic
# compilers.  For example, newer versions of GCC have additional C++ debugging
# flags like _GLIBCXX_DEBUG_PEDANTIC which were not present in older versions.
# --------------------------------------------------------------
LIBMESH_SET_CXX_FLAGS

#------------------------------------------------------
# Checks for code coverage
# this must be done after we look for our own compilers
#------------------------------------------------------
AX_CODE_COVERAGE
libmesh_CFLAGS="$GCOV_FLAGS $libmesh_CFLAGS"
libmesh_CXXFLAGS="$GCOV_FLAGS $libmesh_CXXFLAGS"

#-----------------------------------------------------
# Set compiler flags for devel, opt, etc. methods
#-----------------------------------------------------
LIBMESH_SET_METHODS

# --------------------------------------------------------------
# Test for specific C++11 features.  Note, we still run these
# tests even when C++11 has been explicitly disabled even though
# we know they will fail, since they contain AM_CONDITIONALs.
# --------------------------------------------------------------
LIBMESH_TEST_CXX11_MOVE
LIBMESH_TEST_CXX11_MOVE_CONSTRUCTORS
LIBMESH_TEST_CXX11_RANGEFOR
LIBMESH_TEST_CXX11_DECLTYPE
LIBMESH_TEST_CXX11_RVALUE_REFERENCES
LIBMESH_TEST_CXX11_AUTO
LIBMESH_TEST_CXX11_LAMBDA
LIBMESH_TEST_CXX11_CONSTEXPR
LIBMESH_TEST_CXX11_ALIAS_DECLARATIONS
LIBMESH_TEST_CXX11_SHARED_PTR
LIBMESH_TEST_CXX11_UNIQUE_PTR
LIBMESH_TEST_CXX14_MAKE_UNIQUE
LIBMESH_TEST_CXX11_MAKE_UNIQUE_WORKAROUND
LIBMESH_TEST_CXX11_REGEX
LIBMESH_TEST_CXX11_OVERRIDE
LIBMESH_TEST_CXX11_INITIALIZER_LIST
LIBMESH_TEST_CXX11_VARIADIC_TEMPLATES
LIBMESH_TEST_CXX11_THREAD
LIBMESH_TEST_CXX11_CONDITION_VARIABLE
LIBMESH_TEST_CXX11_TYPE_TRAITS
LIBMESH_TEST_CXX11_INVERSE_HYPERBOLIC_FUNCS
LIBMESH_TEST_CXX11_DELETED_FUNCTIONS
LIBMESH_TEST_CXX11_FINAL
LIBMESH_TEST_CXX11_NULLPTR
LIBMESH_TEST_CXX11_TO_STRING
LIBMESH_TEST_CXX11_NOEXCEPT

#-----------------------------------------------------
# Initialize libtool.  By default, we will build
# only shared libraries on platforms that support them
#-----------------------------------------------------
LT_INIT([disable-static])

#-----------------------------------------------------
# See if the user requested --enable-all-static
#-----------------------------------------------------
AX_ALL_STATIC

# --------------------------------------------------------------
# Release versioning - after we find the C++ compiler
# --------------------------------------------------------------
AX_SPLIT_VERSION

AC_DEFINE_UNQUOTED(MAJOR_VERSION,            [$AX_MAJOR_VERSION],  [libMesh major version number])
AC_DEFINE_UNQUOTED(MINOR_VERSION,            [$AX_MINOR_VERSION],  [libMesh minor version number])
AC_DEFINE_UNQUOTED(MICRO_VERSION,            [$AX_POINT_VERSION],  [libMesh micro version number])
AC_DEFINE_UNQUOTED(LIB_VERSION,              ["$VERSION"],         [libMesh version number])
AC_DEFINE_UNQUOTED(LIB_RELEASE,              ["$BUILD_DEVSTATUS"], [libMesh source code version])
AC_DEFINE_UNQUOTED(CXX,                      ["$CXX"],             [C++ compiler])
AC_DEFINE_UNQUOTED(IO_COMPATIBILITY_VERSION, ["0.7.4"],            [libMesh I/O file format compatiblity string])



# --------------------------------------------------------------
# Portable symbolic links
# --------------------------------------------------------------
AC_PROG_LN_S
# -------------------------------------------------------------

# ---------------------------------------------------------------
# this seems to get done automatically for C++, but not C?
# Required for QHull rules as implemented
# --------------------------------------------------------------
AM_PROG_CC_C_O
# --------------------------------------------------------------

# --------------------------------------------------------------
# parent directories -- portably and without race conditions
# --------------------------------------------------------------
AC_PROG_MKDIR_P
# -------------------------------------------------------------

# --------------------------------------------------------------
# nice sed that doesn't truncate output
# --------------------------------------------------------------
AC_PROG_SED
AC_SUBST(SED)
# --------------------------------------------------------------


# -------------------------------------------------------------
# pkg-config - used to configure installed pakages.  We can
# use it to query our installed targets, if it exists.
# Otherwise we caln fall back to libmesh-config.
#
# Not required to build libmesh, but we can install a config
# file so that users can use 'pkg-config libmesh ...'
# -------------------------------------------------------------
if (test -z "$PKG_CONFIG"); then
  AC_PATH_PROG(PKG_CONFIG, pkg-config, no)
fi
# -------------------------------------------------------------



LIBMESH_COMPILER_FEATURES



# --------------------------------------------------------------
# The AC_SUBST(foo) command creates an output variable from a shell
# variable.  This command also causes AC_OUTPUT to replace all
# instances of @VAR@ with the value of the shell variable VAR in the
# output files (e.g. one or more Makefiles).
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
# configure parts of the libmesh core
# --------------------------------------------------------------
LIBMESH_CORE_FEATURES



# --------------------------------------------------------------
# configure optional, external packages and libraries
# --------------------------------------------------------------
LIBMESH_CONFIGURE_OPTIONAL_PACKAGES


# Query configuration environment
AX_SUMMARIZE_ENV

# Each example gets a makefile
if (test "x$enableexamples" = "xyes"); then
   AC_CONFIG_FILES([
                examples/introduction/introduction_ex1/Makefile
                examples/introduction/introduction_ex2/Makefile
                examples/introduction/introduction_ex3/Makefile
                examples/introduction/introduction_ex4/Makefile
                examples/introduction/introduction_ex5/Makefile
                examples/adaptivity/adaptivity_ex1/Makefile
                examples/adaptivity/adaptivity_ex2/Makefile
                examples/adaptivity/adaptivity_ex3/Makefile
                examples/adaptivity/adaptivity_ex4/Makefile
                examples/adaptivity/adaptivity_ex5/Makefile
                examples/adjoints/adjoints_ex1/Makefile
                examples/adjoints/adjoints_ex2/Makefile
                examples/adjoints/adjoints_ex3/Makefile
                examples/adjoints/adjoints_ex4/Makefile
                examples/adjoints/adjoints_ex5/Makefile
                examples/eigenproblems/eigenproblems_ex1/Makefile
                examples/eigenproblems/eigenproblems_ex2/Makefile
                examples/eigenproblems/eigenproblems_ex3/Makefile
                examples/fem_system/fem_system_ex1/Makefile
                examples/fem_system/fem_system_ex2/Makefile
                examples/fem_system/fem_system_ex3/Makefile
                examples/fem_system/fem_system_ex4/Makefile
                examples/solution_transfer/solution_transfer_ex1/Makefile
                examples/miscellaneous/miscellaneous_ex1/Makefile
                examples/miscellaneous/miscellaneous_ex2/Makefile
                examples/miscellaneous/miscellaneous_ex3/Makefile
                examples/miscellaneous/miscellaneous_ex4/Makefile
                examples/miscellaneous/miscellaneous_ex5/Makefile
                examples/miscellaneous/miscellaneous_ex6/Makefile
                examples/miscellaneous/miscellaneous_ex7/Makefile
                examples/miscellaneous/miscellaneous_ex8/Makefile
                examples/miscellaneous/miscellaneous_ex9/Makefile
                examples/miscellaneous/miscellaneous_ex10/Makefile
                examples/miscellaneous/miscellaneous_ex11/Makefile
                examples/miscellaneous/miscellaneous_ex12/Makefile
                examples/optimization/optimization_ex1/Makefile
                examples/optimization/optimization_ex2/Makefile
                examples/subdomains/subdomains_ex1/Makefile
                examples/subdomains/subdomains_ex2/Makefile
                examples/subdomains/subdomains_ex3/Makefile
                examples/systems_of_equations/systems_of_equations_ex1/Makefile
                examples/systems_of_equations/systems_of_equations_ex2/Makefile
                examples/systems_of_equations/systems_of_equations_ex3/Makefile
                examples/systems_of_equations/systems_of_equations_ex4/Makefile
                examples/systems_of_equations/systems_of_equations_ex5/Makefile
                examples/systems_of_equations/systems_of_equations_ex6/Makefile
                examples/systems_of_equations/systems_of_equations_ex7/Makefile
                examples/systems_of_equations/systems_of_equations_ex8/Makefile
                examples/reduced_basis/reduced_basis_ex1/Makefile
                examples/reduced_basis/reduced_basis_ex2/Makefile
                examples/reduced_basis/reduced_basis_ex3/Makefile
                examples/reduced_basis/reduced_basis_ex4/Makefile
                examples/reduced_basis/reduced_basis_ex5/Makefile
                examples/reduced_basis/reduced_basis_ex6/Makefile
                examples/reduced_basis/reduced_basis_ex7/Makefile
                examples/transient/transient_ex1/Makefile
                examples/transient/transient_ex2/Makefile
                examples/vector_fe/vector_fe_ex1/Makefile
                examples/vector_fe/vector_fe_ex2/Makefile
                examples/vector_fe/vector_fe_ex3/Makefile
                examples/vector_fe/vector_fe_ex4/Makefile
                examples/Makefile
                ])
fi


# hackery.  If we are supporting nested autoconf packages and we want to specify
# additional arguments to be passed to those packages, do that here.
# Specifically, we append libmesh_subpackage_arguments to ac_configure_args
# before AC_OUTPUT recurses into our subpackages
if (test "x$enablenested" = "xyes"); then
  ac_configure_args_SAVE="$ac_configure_args"
  ac_configure_args="$ac_configure_args $libmesh_subpackage_arguments CXX='$CXX' CC='$CC' F77='$F77' FC='$FC' CPPFLAGS='$HDF5_CPPFLAGS $CPPFLAGS' LIBS='$HDF5_LIBS $LIBS'"
fi

# Create output files. Also configures any subpackages
AC_OUTPUT

# end hackery
if (test "x$enablenested" = "xyes"); then
  ac_configure_args="$ac_configure_args_SAVE"
fi

# Final summary
AX_SUMMARIZE_CONFIG



AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(--------- Done Configuring libMesh ----------)
AC_MSG_RESULT(---------------------------------------------)
# rm -f -r conftest* config.cache include/libmesh/libmesh_config.h.tmp
