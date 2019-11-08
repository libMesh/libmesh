 dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(timpi, 0.0.0, roystgnr@ices.utexas.edu)
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(src/utilities/include/timpi/timpi_config.h.tmp)
dnl Prefixes everything in the temporary file with TIMPI_
AX_PREFIX_CONFIG_H(src/utilities/include/timpi/timpi_config.h,TIMPI,src/utilities/include/timpi/timpi_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar 1.11 no-define color-tests parallel-tests subdir-objects])
AM_SILENT_RULES(yes)  # use silent rules where available - automake 1.11

dnl Turn off AM_MAINTAINER_MODE for bootstrapped release by
dnl uncommenting this:
AM_MAINTAINER_MODE # git trashes timestamps, forcing autoreconf needlessly

dnl Release versioning
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION
GENERIC_LIBRARY_NAME=timpi

dnl shared library versioning

PACKAGE=$GENERIC_LIBRARY_NAME
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

dnl AC_SUBST(GENERIC_LIBRARY_VERSION)
dnl AC_SUBST(GENERIC_LIBRARY_NAME)

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])



dnl------------------------------
dnl Checks for compilers
dnl------------------------------

TIMPI_SET_COMPILERS

# --------------------------------------------------------------
# Standard autoconf macro for determining the proper -std=c++11
# flag for the current compiler.  Adds the result to CXXFLAGS if
# one is found.  See ax_cxx_compile_stdcxx_11.m4 for details.
# --------------------------------------------------------------
enablecxx11=yes
AX_CXX_COMPILE_STDCXX_11([],[mandatory])

#-----------------------------------------------------
# Set compiler flags for devel, opt, etc. methods
#-----------------------------------------------------
TIMPI_SET_CXX_FLAGS

SET_BUILD_METHODS

dnl---------------------------------------------------------
dnl Checks for library prerequisites for other libraries...
dnl---------------------------------------------------------

dnl------------------------------------------------------
dnl If the user didn't --disable-mpi, try hard to find it
dnl------------------------------------------------------
timpi_optional_INCLUDES=""
timpi_optional_LIBS=""

AS_IF([test "x$enablempi" = xyes],
      [
        ACX_MPI
        AS_IF([test "x$enablempi" = xno],
              [
dnl Double-check for MPI via PETSc; libMesh users are accustomed to
dnl being able to autodetect MPI support this way.
                CONFIGURE_PETSC_MPI

                AS_IF([test $enablepetsc != no],
                      [
                        timpi_optional_INCLUDES="$PETSCINCLUDEDIRS $timpi_optional_INCLUDES"
                        timpi_optional_LIBS="$PETSCLINKLIBS $timpi_optional_LIBS"
                      ])
              ],
              [
                timpi_optional_INCLUDES="-I$MPI_INCLUDES_PATH $timpi_optional_INCLUDES"
                timpi_optional_LIBS="-L$MPI_LIBS_PATH $LIBS $timpi_optional_LIBS"
              ])
      ])

AC_SUBST(timpi_optional_INCLUDES)
AC_SUBST(timpi_optional_LIBS)

dnl---------------------------------------------------------
dnl Sets library options
dnl---------------------------------------------------------

TIMPI_OPTIONAL_FEATURES

dnl-----------------------------------------------------
dnl Initialize libtool.  By default, we will build
dnl only shared libraries on platforms that support them
dnl-----------------------------------------------------
LT_INIT([disable-static])

dnl---------------------------------
dnl Query configuration environment
dnl---------------------------------
AX_SUMMARIZE_ENV

dnl-----------------------------------------------
dnl Generate header files
dnl-----------------------------------------------
AC_CONFIG_FILES(src/utilities/include/timpi/timpi_version.h)

dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([
  Makefile
  timpi.pc
  src/Makefile
  test/Makefile
  ])

# Must still call AC_OUTPUT() after generating all the files
AC_OUTPUT()

dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_ENV

AX_SUMMARIZE_CONFIG
