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

# --------------------------------------------------------------
# Set our preferred compilers
# --------------------------------------------------------------
LIBMESH_SET_COMPILERS

# --------------------------------------------------------------
# Allow the user to explicitly turn off C++11 support by saying
# --disable-cxx11?  This is no longer supported.
#
# We'll leave the "option" in place for a release, though, just in
# case anyone is specifying --enable-cxx11 in their build scripts
# --------------------------------------------------------------
AC_ARG_ENABLE(cxx11,
              AC_HELP_STRING([--disable-cxx11],
                             [build without C++11 support (no longer supported!)]),
              [AS_CASE("${enableval}",
                       [yes], [enablecxx11=yes],
                       [no],  [AC_MSG_ERROR(libMesh now requires C++11)],
                       [AC_MSG_ERROR(bad value ${enableval} for --disable-cxx11)])],
              [enablecxx11=irrelevant])

AS_IF([test $enablecxx11 != irrelevant],
      [AC_MSG_WARN([--enable/disable-cxx11 are now deprecated])])
enablecxx11=yes

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
              [AS_CASE("${enableval}",
                       [yes], [cxx11required=yes],
                       [no],  [cxx11required=no],
                       [AC_MSG_ERROR(bad value ${enableval} for --enable-cxx11-required)])],
              [cxx11required=irrelevant])

AS_IF([test $cxx11required != irrelevant],
      [AC_MSG_WARN([C++11 is now always required for libMesh])])

# --------------------------------------------------------------
# Standard autoconf macro for determining the proper -std=c++11
# flag for the current compiler.  Adds the result to CXXFLAGS if
# one is found.  See ax_cxx_compile_stdcxx_11.m4 for details.
# --------------------------------------------------------------
AX_CXX_COMPILE_STDCXX_11([],[mandatory])

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
# Test for mandatory C++11 features, those which libMesh now requires
# for internal use.
# --------------------------------------------------------------
LIBMESH_TEST_CXX11_AUTO
AS_IF([test "x$have_cxx11_auto" != "xyes"],
      [AC_MSG_ERROR([libMesh requires C++11 auto variable support])])

LIBMESH_TEST_CXX11_RANGEFOR
AS_IF([test "x$have_cxx11_rangefor" != "xyes"],
      [AC_MSG_ERROR([libMesh requires C++11 ranged-for loop support])])

# Technically we require uniform initialization, but this is the
# closest test we have
LIBMESH_TEST_CXX11_INITIALIZER_LIST
AS_IF([test "x$have_cxx11_initializer_list" != "xyes"],
      [AC_MSG_ERROR([libMesh requires C++11 initializer list support])])

LIBMESH_TEST_CXX11_UNIQUE_PTR
AS_IF([test "x$have_cxx11_unique_ptr" != "xyes"],
      [AC_MSG_ERROR([libMesh requires C++11 std::unique_ptr support])])

LIBMESH_TEST_CXX11_MAKE_UNIQUE_WORKAROUND
AS_IF([test "x$have_cxx11_make_unique_workaround" != "xyes"],
      [AC_MSG_ERROR([libMesh requires C++11 variadic template support])])

# --------------------------------------------------------------
# Test for optional modern C++ features, which libMesh offers shims
# for and/or which libMesh applications may want to selectively use.
# --------------------------------------------------------------
LIBMESH_TEST_CXX11_MOVE
LIBMESH_TEST_CXX11_MOVE_CONSTRUCTORS
LIBMESH_TEST_CXX11_DECLTYPE
LIBMESH_TEST_CXX11_RVALUE_REFERENCES
LIBMESH_TEST_CXX11_LAMBDA
LIBMESH_TEST_CXX11_CONSTEXPR
LIBMESH_TEST_CXX11_ALIAS_DECLARATIONS
LIBMESH_TEST_CXX11_SHARED_PTR
LIBMESH_TEST_CXX14_MAKE_UNIQUE
LIBMESH_TEST_CXX11_REGEX
LIBMESH_TEST_CXX11_OVERRIDE
LIBMESH_TEST_CXX11_VARIADIC_TEMPLATES
LIBMESH_TEST_CXX11_THREAD
LIBMESH_TEST_CXX11_CONDITION_VARIABLE
LIBMESH_TEST_CXX11_TYPE_TRAITS
LIBMESH_TEST_CXX11_MATH_FUNCS
LIBMESH_TEST_CXX11_DELETED_FUNCTIONS
LIBMESH_TEST_CXX11_FINAL
LIBMESH_TEST_CXX11_NULLPTR
LIBMESH_TEST_CXX11_TO_STRING
LIBMESH_TEST_CXX11_NOEXCEPT
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

AC_DEFINE_UNQUOTED(MAJOR_VERSION,            [$AX_MAJOR_VERSION],  [libMesh major version number])
AC_DEFINE_UNQUOTED(MINOR_VERSION,            [$AX_MINOR_VERSION],  [libMesh minor version number])
AC_DEFINE_UNQUOTED(MICRO_VERSION,            [$AX_POINT_VERSION],  [libMesh micro version number])
AC_DEFINE_UNQUOTED(LIB_VERSION,              ["$VERSION"],         [libMesh version number])
AC_DEFINE_UNQUOTED(LIB_RELEASE,              ["$BUILD_DEVSTATUS"], [libMesh source code version])
AC_DEFINE_UNQUOTED(CXX,                      ["$CXX"],             [C++ compiler])
AC_DEFINE_UNQUOTED(IO_COMPATIBILITY_VERSION, ["0.7.4"],            [libMesh I/O file format compatibility string])



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
# configure optional, external packages and libraries
# --------------------------------------------------------------
LIBMESH_CONFIGURE_OPTIONAL_PACKAGES


# Query configuration environment
AX_SUMMARIZE_ENV

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
                         examples/miscellaneous/miscellaneous_ex13/Makefile
                         examples/miscellaneous/miscellaneous_ex14/Makefile
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
                         examples/Makefile])
     ])



# hackery.  If we are supporting nested autoconf packages and we want to specify
# additional arguments to be passed to those packages, do that here.
# Specifically, we append libmesh_subpackage_arguments to ac_configure_args
# before AC_OUTPUT recurses into our subpackages
AS_IF([test "x$enablenested" = "xyes"],
      [
        ac_configure_args_SAVE="$ac_configure_args"
        ac_configure_args="$ac_configure_args $libmesh_subpackage_arguments CXX='$CXX' CC='$CC' F77='$F77' FC='$FC' CPPFLAGS='$HDF5_CPPFLAGS $CPPFLAGS' LIBS='$HDF5_LIBS $LIBS'"
      ])

# Create output files. Also configures any subpackages
AC_OUTPUT

# end hackery
AS_IF([test "x$enablenested" = "xyes"],
      [ac_configure_args="$ac_configure_args_SAVE"])

# Final summary
AX_SUMMARIZE_CONFIG



AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(--------- Done Configuring libMesh ----------)
AC_MSG_RESULT(---------------------------------------------)
# rm -f -r conftest* config.cache include/libmesh/libmesh_config.h.tmp
