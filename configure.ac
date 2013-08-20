# --------------------------------------------------------------
# $Id$
# --------------------------------------------------------------



AC_INIT([libmesh], [0.9.3-pre], [libmesh-users@lists.sourceforge.net], [libmesh], [http://libmesh.sourceforge.net])
AC_PREREQ([2.61]) # require autoconf 2.61 - I havent tried anything older

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

AC_REVISION($Revision$)

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
                 doc/Doxyfile
                 contrib/utils/libmesh.pc:contrib/utils/libmesh-opt.pc.in
                 contrib/utils/libmesh-opt.pc
                 contrib/utils/libmesh-dbg.pc
                 contrib/utils/libmesh-devel.pc
                 contrib/utils/libmesh-prof.pc
                 contrib/utils/libmesh-oprof.pc
                 ])

AC_CONFIG_FILES(contrib/bin/libmesh-config,          [chmod +x contrib/bin/libmesh-config])
#AC_CONFIG_FILES(contrib/bin/strip_dup_libs.pl,       [chmod +x contrib/bin/strip_dup_libs.pl])
#AC_CONFIG_FILES(contrib/bin/strip_dup_incl_paths.pl, [chmod +x contrib/bin/strip_dup_incl_paths.pl])
AC_CONFIG_FILES(contrib/bin/xda2mgf,                 [chmod +x contrib/bin/xda2mgf])
AC_CONFIG_FILES(contrib/bin/create_libmesh_release,  [chmod +x contrib/bin/create_libmesh_release])
AC_CONFIG_FILES(doc/dot,                             [chmod +x doc/dot])

# Create a temporary file with e.g. HAVE_MPI
AC_CONFIG_HEADER(include/libmesh_config.h.tmp:include/libmesh_config.h.in)
# Prefixes everything in the temporary file with LIBMESH_, e.g. LIBMESH_HAVE_MPI
AX_PREFIX_CONFIG_H(include/libmesh_config.h,LIBMESH,include/libmesh_config.h.tmp)

# Print configuration info to the header file
AC_DEFINE_UNQUOTED(CONFIGURE_INFO, ["$0 $ac_configure_args"], [Configuration information.])

# require automake 1.11 - color tests
AM_INIT_AUTOMAKE([dist-xz dist-bzip2 tar-ustar color-tests 1.11])

# use silent rules - automake 1.11
AM_SILENT_RULES(yes)

# by default we want 'maintainer mode', which means automake will generate rules to recreate e.g. Makefile.in from
# Makefile.am. But this adds support for --disable-maintainer-mode, which may be useful for forcing configure to
# work on systems with no autotools
AM_MAINTAINER_MODE([enable])


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

#-----------------------------------------------------
# Initialize libtool.  By default, we will build
# only shared libraries on platforms that support them
#-----------------------------------------------------
LT_INIT([disable-static])



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
# Perl - look for perl.  Perl is used to generate dependencies
# -------------------------------------------------------------
AC_PATH_PROG(PERL, perl)
AC_SUBST(PERL)
# -------------------------------------------------------------



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
                examples/subdomains/subdomains_ex1/Makefile
                examples/subdomains/subdomains_ex2/Makefile
                examples/systems_of_equations/systems_of_equations_ex1/Makefile
                examples/systems_of_equations/systems_of_equations_ex2/Makefile
                examples/systems_of_equations/systems_of_equations_ex3/Makefile
                examples/systems_of_equations/systems_of_equations_ex4/Makefile
                examples/systems_of_equations/systems_of_equations_ex5/Makefile
                examples/systems_of_equations/systems_of_equations_ex6/Makefile
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
