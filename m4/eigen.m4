# ----------------------------------------------------------------
# Locate header files for the C++ linear algebra library Eigen.
# Eigen is a header-only template library. By default we check for the
# Eigen files in the --with-eigen-include=xxx argument provided to
# configure, or if those don't exist in the $EIGEN_INC/Eigen directory,
# or in /usr/include.
#
# If an argument to this macro is specified, it should give the minimum useful
# version number.
#
# If a second argument to this macro is specified as 'yes', then this is a
# required dependency and configure will exit with an error if it is not
# satisfied.
#
# Note: Eigen is installed (by default) at the location
# /path/to/eigen/Eigen, i.e. with path ending in capital 'Eigen'.
# You should specify --with-eigen-include=/path/to/eigen
# during configure, or set your $EIGEN_INC environment variable
# to /path/to/eigen.
# ----------------------------------------------------------------

AC_DEFUN([CONFIGURE_EIGEN],
[
  AC_ARG_ENABLE(eigen,
                AS_HELP_STRING([--disable-eigen],
                               [build without Eigen linear algebra support]),
                [AS_CASE("${enableval}",
                         [yes], [enableeigen=yes],
                         [no],  [enableeigen=no],
                         [AC_MSG_ERROR(bad value ${enableval} for --enable-eigen)])],
                [enableeigen=$enableoptional])

  dnl package requirement; if not specified, the default is to assume that
  dnl the package is optional
  is_package_required=ifelse([$2], ,no, $2 )

  install_internal_eigen=no
  AS_IF([test "x$enableeigen" = "xyes"],
        [
          dnl User-specific include path
          AC_ARG_WITH(eigen-include,
                      AS_HELP_STRING([--with-eigen-include=PATH],[Specify the path for EIGEN header files]),
                      witheigeninc=$withval,
                      witheigeninc=no)

          dnl Fall back on default paths to Eigen's include files
          AS_IF([test "x$witheigeninc" != "xno"],                                        [EIGEN_INC="$witheigeninc"],
                [test "x$EIGEN_INC" != "x" && test -f $EIGEN_INC/Eigen/Eigen],           [AS_ECHO(["Environment EIGEN_INC=$EIGEN_INC"])],
                [test "x$EIGEN3_INCLUDE" != "x" && test -f $EIGEN3_INCLUDE/Eigen/Eigen], [EIGEN_INC=$EIGEN3_INCLUDE
                                                                                          AS_ECHO(["Environment EIGEN3_INCLUDE=$EIGEN_INC"])],
                [test "x$EIGEN_INCLUDE" != "x" && test -f $EIGEN_INCLUDE/Eigen/Eigen],   [EIGEN_INC=$EIGEN_INCLUDE
                                                                                          AS_ECHO(["Environment EIGEN_INCLUDE=$EIGEN_INC"])],
                [test -f /usr/include/eigen3/Eigen/Eigen],                               [EIGEN_INC="/usr/include/eigen3"
                                                                                          AS_ECHO(["System EIGEN_INC=$EIGEN_INC"])],
                [EIGEN_INC="/usr/include"
                 AS_ECHO(["Testing EIGEN_INC=$EIGEN_INC"])])

          dnl Initialize Makefile/config.h substitution variables
          EIGEN_INCLUDE=""

          dnl Check for existence of a header file in the specified location.  Note: here
          dnl we are checking for the header file "Eigen" in the Eigen directory.
          AC_LANG_SAVE
          AC_LANG_CPLUSPLUS

          externaleigenincFound=no;
          AC_CHECK_HEADERS($EIGEN_INC/Eigen/Eigen, externaleigenincFound=yes)

          dnl Check to make sure the external header files are sufficiently up
          dnl to date - this fixes our Eigen detection on Scientific Linux 6
          AS_IF([test "x$externaleigenincFound" = "xyes"],
                [
                  enableeigenincFound=yes
                  ac_eigen_save_CPPFLAGS="$CPPFLAGS"
                  CPPFLAGS="-I${EIGEN_INC} ${CPPFLAGS}"
                  AC_CHECK_HEADERS([Eigen/Dense],[],[enableeigenincFound=no])
                  AS_IF([test "x$enableeigensparse" = "xyes"], [AC_CHECK_HEADERS([Eigen/Sparse],[],[enableeigenincFound=no])])
                  CPPFLAGS="${ac_eigen_save_CPPFLAGS}"
                ])

          AS_IF([test "x$enableeigenincFound" = "xyes"],   [EIGEN_INCLUDE="-I$EIGEN_INC"],
                [test -d $top_srcdir/contrib/eigen/eigen], [AC_MSG_RESULT([<<< external Eigen header files not found, using Eigen from ./contrib >>>])
                                                            EIGEN_INC=$top_srcdir/contrib/eigen/eigen
                                                            EIGEN_INCLUDE="-I\$(top_srcdir)/contrib/eigen/eigen"
                                                            install_internal_eigen=yes],
                [enableeigen=no])

          dnl OK, we have a usable eigen path, make sure the headers we want are good.
          AS_IF([test "x$enableeigen" = "xyes"],
                [
                  ac_eigen_save_CPPFLAGS="$CPPFLAGS"
                  CPPFLAGS="-I${EIGEN_INC} ${CPPFLAGS}"

                  dnl Do not use cached results for the header checks
                  AS_UNSET([ac_cv_header_Eigen_Dense])
                  AS_UNSET([ac_cv_header_Eigen_Sparse])

                  AC_CHECK_HEADERS([Eigen/Dense],[],[enableeigen=no])

                  AS_IF([test "x$enableeigensparse" = "xyes"],
                        [AC_CHECK_HEADERS([Eigen/Sparse],[],[enableeigen=no])])

                    dnl -----------------------
                    dnl  Minimum version check
                    dnl ----------------------

                    min_eigen_version=ifelse([$1], ,0.0.0, $1)

                    dnl looking for major.minor.micro (which Eigen calls world.major.minor) style versioning
                    MAJOR_VER=`echo $min_eigen_version | sed 's/^\([[0-9]]*\).*/\1/'`
                    AS_IF([test "x${MAJOR_VER}" = "x"], [MAJOR_VER=0])

                    MINOR_VER=`echo $min_eigen_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
                    AS_IF([test "x${MINOR_VER}" = "x"], [MINOR_VER=0])

                    MICRO_VER=`echo $min_eigen_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
                    AS_IF([test "x${MICRO_VER}" = "x"], [MICRO_VER=0])

                    AC_MSG_CHECKING(for eigen - version >= $min_eigen_version)

                    AC_LANG_PUSH([C++])
                    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
                    @%:@include "Eigen/Core"
                        ]], [[
                        @%:@if EIGEN_WORLD_VERSION > $MAJOR_VER
                        @%:@elif (EIGEN_WORLD_VERSION >= $MAJOR_VER) && (EIGEN_MAJOR_VERSION > $MINOR_VER)
                        @%:@elif (EIGEN_WORLD_VERSION >= $MAJOR_VER) && (EIGEN_MAJOR_VERSION >= $MINOR_VER) && (EIGEN_MINOR_VERSION >= $MICRO_VER)
                        @%:@else
                        @%:@  error version is too old
                        @%:@endif
                    ]])],[
                        AC_MSG_RESULT(yes)
                    ],[
                        AC_MSG_RESULT(no)
                        enableeigen=no
                    ])
                    AC_LANG_POP([C++])

                    CPPFLAGS="${ac_eigen_save_CPPFLAGS}"

                  dnl if we survived, we really have Eigen
                  AS_IF([test "x$enableeigen" = "xyes"],         [HAVE_EIGEN=1
                                                                  AC_DEFINE(HAVE_EIGEN, 1, [Flag indicating whether the library will be compiled with Eigen support])
                                                                  AC_MSG_RESULT(<<< Configuring library with Eigen support >>>)],
                        [test "x$is_package_required" = "xyes"], [AC_MSG_ERROR([Your EIGEN version ($EIGEN_INC) does not meet the minimum versioning
                                                                                requirements ($min_eigen_version).  Please use --with-eigen-include to
                                                                                specify the location of an updated installation.])])
                ])
          AC_LANG_RESTORE
        ])

  dnl Substitute the substitution variables
  AC_SUBST(EIGEN_INCLUDE)
])
