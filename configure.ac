# --------------------------------------------------------------
# configure.ac file for libmesh.  Used by Automake to generate
# aclocal.m4 and Makefile.in. See also: http://tinyurl.com/hatuyfn
# --------------------------------------------------------------
AC_INIT([libmesh],[1.8.0-pre],[https://github.com/libMesh/libmesh/discussions],[libmesh],[http://libmesh.github.io])
AC_PREREQ([2.68])

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
AC_CONFIG_HEADERS([include/libmesh_config.h.tmp:include/libmesh_config.h.in])
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
AS_IF([test $libmesh_vpath_build = yes],
      [AC_MSG_RESULT(<<< Configuring build directory for VPATH build >>>)])
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

#--------------------------------------------------------------------
# Check for whether to enable MPI and PETSc which will influence our
# choice of compiler
#--------------------------------------------------------------------
ACSM_COMPILER_CONTROL_ARGS

#-----------------------------------------------------------------------
# Scrape PETSc configure information for their CXX, MPI_INCLUDE, MPI_LIB,
# PETSCLINKLIBS, PETSCINCLUDEDIRS, and other variables
#-----------------------------------------------------------------------
ACSM_SCRAPE_PETSC_CONFIGURE

# --------------------------------------------------------------
# Set our preferred compilers
# --------------------------------------------------------------
LIBMESH_SET_COMPILERS

# By default we do not set any -march flag, but this can be turned on
# by configuring with --enable-march
AC_ARG_ENABLE(march,
              AS_HELP_STRING([--enable-march],[Determine the correct -march flag for this system]),
              [AS_CASE("${enableval}",
                       [yes], [enablemarch=yes],
                       [no],  [enablemarch=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --enable-march)])],
              [enablemarch=no])

# Set the value of the cache variable $ax_cv_gcc_archflag based
# on the current architecture. For the first argument:
# [yes] = generate portable code (-mtune)
# [no]  = generate architecture-specific code (-march)
# This script is specific to GCC; the flags it determines might not
# work for other compilers. There is also an optional third
# argument which controls what happens when determining an arch flag
# fails, we don't do anything in that case.
AS_IF([test "x$enablemarch" != "xno"],
      [
        AX_GCC_ARCHFLAG([no],
                        [
                          dnl Action if success - append $ax_cv_gcc_archflag to CXXFLAGS flags for
                          dnl all but debugging METHODS. I don't think it is needed with -O0 when
                          dnl debugging.
                          CXXFLAGS_OPT="$CXXFLAGS_OPT $ax_cv_gcc_archflag"
                          CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL $ax_cv_gcc_archflag"
                          CXXFLAGS_PROF="$CXXFLAGS_PROF $ax_cv_gcc_archflag"
                          CXXFLAGS_OPROF="$CXXFLAGS_OPROF $ax_cv_gcc_archflag"
                        ])
      ])

# --------------------------------------------------------------
# Autoconf macro for determining the proper -std=c++??
# flag, for the current compiler, for the user's requested C++
# standards level.  Adds the required flag to CXXFLAGS if
# one is found.  Exits if no acceptable flag is found.
# --------------------------------------------------------------
ACSM_CXX_COMPILER_STANDARD([2017], [2017])

# --------------------------------------------------------------
# See compiler.m4 for the definition of this function.  It uses the
# output of DETERMINE_CXX_BRAND to set special flags for idiosyncratic
# compilers.  For example, newer versions of GCC have additional C++ debugging
# flags like _GLIBCXX_DEBUG_PEDANTIC which were not present in older versions.
# --------------------------------------------------------------
LIBMESH_SET_BUILD_FLAGS

#------------------------------------------------------
# Checks for code coverage
# this must be done after we look for our own compilers
#------------------------------------------------------
ACSM_CODE_COVERAGE
libmesh_CFLAGS="$GCOV_FLAGS $libmesh_CFLAGS"
libmesh_CXXFLAGS="$GCOV_FLAGS $libmesh_CXXFLAGS"

#-----------------------------------------------------
# Set compiler flags for devel, opt, etc. methods
#-----------------------------------------------------
LIBMESH_SET_METHODS

AC_CONFIG_FILES([contrib/utils/libmesh.pc:$LIBMESH_PC_IN])

ACSM_ENABLE_WERROR

ACSM_ENABLE_PARANOID

# --------------------------------------------------------------
# Test for optional modern C++ features, which libMesh offers shims
# for and/or which libMesh applications may want to selectively use.
# --------------------------------------------------------------
LIBMESH_TEST_CXX14_MAKE_UNIQUE
LIBMESH_TEST_CXX11_REGEX
LIBMESH_TEST_CXX11_THREAD
LIBMESH_TEST_CXX11_CONDITION_VARIABLE
LIBMESH_TEST_CXX11_TYPE_TRAITS
LIBMESH_TEST_CXX11_MATH_FUNCS
LIBMESH_TEST_CXX17_FALLTHROUGH_ATTRIBUTE

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

# Strip any '-pre' suffix from version; our code expects a number
LIBMESH_NUMERIC_POINT_VERSION=`echo "$AX_POINT_VERSION" | $SED 's/-.*//'`

AC_DEFINE_UNQUOTED(MAJOR_VERSION,            [$AX_MAJOR_VERSION],  [libMesh major version number])
AC_DEFINE_UNQUOTED(MINOR_VERSION,            [$AX_MINOR_VERSION],  [libMesh minor version number])
AC_DEFINE_UNQUOTED(MICRO_VERSION,            [$LIBMESH_NUMERIC_POINT_VERSION],  [libMesh micro version number])
AC_DEFINE_UNQUOTED(LIB_VERSION,              ["$VERSION"],         [libMesh version number])
AC_DEFINE_UNQUOTED(LIB_RELEASE,              ["$BUILD_DEVSTATUS"], [libMesh source code version])
AC_DEFINE_UNQUOTED(CXX,                      ["$CXX"],             [C++ compiler])
AC_DEFINE_UNQUOTED(IO_COMPATIBILITY_VERSION, ["1.7.0"],            [libMesh I/O file format compatibility string])



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
# pkg-config - used to configure installed packages.  We can
# use it to query our installed targets, if it exists.
# Otherwise we can fall back to libmesh-config.
#
# Not required to build libmesh, but we can install a config
# file so that users can use 'pkg-config libmesh ...'
# -------------------------------------------------------------
AS_IF([test -z "$PKG_CONFIG"],
      [AC_PATH_PROG(PKG_CONFIG, pkg-config, no)])
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
AS_IF([test "$CXX" = xlC_r],
      [
        AC_CHECK_PROG(MPCXX,mpCC_r,mpCC_r,$CXX,$PATH)
        CXX=$MPCXX
        AS_IF([test "$CXX" = mpCC_r],
              [AC_MSG_RESULT(<<< Using $CXX instead of xlC_r to compile C++ >>>)])
      ])

AS_IF([test "$CXX" = xlC],
      [
        AC_CHECK_PROG(MPCXX,mpCC,mpCC,$CXX,$PATH)
        CXX=$MPCXX
        AS_IF([test "$CXX" = mpCC],
              [AC_MSG_RESULT(<<< Using $CXX instead of xlC to compile C++ >>>)])
      ])
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
# configure TIMPI, which is practically libMesh core
#
# Use AX_SUBDIRS_CONFIGURE here so we get timpi_config.h before
# contrib/metaphysicl might need it later.
# --------------------------------------------------------------
AX_SUBDIRS_CONFIGURE([contrib/timpi])


# --------------------------------------------------------------
# configure optional, external packages and libraries
# --------------------------------------------------------------
LIBMESH_CONFIGURE_OPTIONAL_PACKAGES


# --------------------------------------------------------------
# With both configure-requested (LIBMESH_CORE_FEATURES) and
# PETSc-provided DOF_ID_BYTES information, we can finally set that
# value with no risk we'll have to redefine it later - any
# redefinition via AC_DEFINE breaks configure tests with e.g.
# g++ and user-specified CXXFLAGS=-Werror
# --------------------------------------------------------------
AS_CASE("$dof_bytes",
        [1], [AC_DEFINE(DOF_ID_BYTES, 1, [size of dof_id])],
        [2], [AC_DEFINE(DOF_ID_BYTES, 2, [size of dof_id])],
        [4], [AC_DEFINE(DOF_ID_BYTES, 4, [size of dof_id])],
        [8], [AC_DEFINE(DOF_ID_BYTES, 8, [size of dof_id])],
        [
          AC_MSG_ERROR([>>> unrecognized dof_id size: $dof_bytes - only 1,2,4,8 are supported])
        ])

# Query configuration environment
ACSM_SUMMARIZE_ENV

# Each example gets a makefile
AS_IF([test "x$enableexamples" = "xyes"],
      [
        AC_CONFIG_FILES([examples/introduction/introduction_ex1/Makefile
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
                         examples/adjoints/adjoints_ex6/Makefile
                         examples/adjoints/adjoints_ex7/Makefile
                         examples/eigenproblems/eigenproblems_ex1/Makefile
                         examples/eigenproblems/eigenproblems_ex2/Makefile
                         examples/eigenproblems/eigenproblems_ex3/Makefile
                         examples/eigenproblems/eigenproblems_ex4/Makefile
                         examples/fem_system/fem_system_ex1/Makefile
                         examples/fem_system/fem_system_ex2/Makefile
                         examples/fem_system/fem_system_ex3/Makefile
                         examples/fem_system/fem_system_ex4/Makefile
                         examples/fem_system/fem_system_ex5/Makefile
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
                         examples/miscellaneous/miscellaneous_ex13/Makefile
                         examples/miscellaneous/miscellaneous_ex14/Makefile
                         examples/miscellaneous/miscellaneous_ex15/Makefile
                         examples/miscellaneous/miscellaneous_ex16/Makefile
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
                         examples/systems_of_equations/systems_of_equations_ex9/Makefile
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
                         examples/vector_fe/vector_fe_ex5/Makefile
                         examples/vector_fe/vector_fe_ex6/Makefile
                         examples/vector_fe/vector_fe_ex7/Makefile
                         examples/vector_fe/vector_fe_ex8/Makefile
                         examples/vector_fe/vector_fe_ex9/Makefile
                         examples/Makefile])
     ])

# Create output files. Also configures any subpackages
AC_OUTPUT

# Final summary
AX_SUMMARIZE_CONFIG



AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(--------- Done Configuring libMesh ----------)
AC_MSG_RESULT(---------------------------------------------)
# rm -f -r conftest* config.cache include/libmesh/libmesh_config.h.tmp
