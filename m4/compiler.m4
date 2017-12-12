# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_SET_COMPILERS],
[
  # --------------------------------------------------------------
  # look for a decent C++ compiler or honor --with-cxx=...
  CXX_TRY_LIST="g++ icpc icc pgCC c++"

     # -------------------------------------------------------------------
     # MPI -- enabled by default.  Check for it now so we can be somewhat
     #                             smart about which compilers to look for
     # -------------------------------------------------------------------
     AC_ARG_ENABLE(mpi,
                   AS_HELP_STRING([--disable-mpi],
                                  [build without MPI message passing support]),
                   [case "${enableval}" in
                     yes)  enablempi=yes ;;
                     no)  enablempi=no ;;
                     *)  AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
                   esac],
                   [enablempi=yes])

  if  (test "$enablempi" != no) ; then
    CXX_TRY_LIST="mpicxx mpiCC mpicc $CXX_TRY_LIST"
  else
    AC_MSG_RESULT(>>> Disabling MPI per user request <<<)
  fi

  AC_ARG_WITH([cxx],
              AS_HELP_STRING([--with-cxx=CXX], [C++ compiler to use]),
              [CXX="$withval"],
              [])

  # --------------------------------------------------------------
  # Determines a C++ compiler to use.  First checks if the variable CXX is
  # already set.  If not, then searches under g++, c++, and other names.
  # --------------------------------------------------------------
  AC_PROG_CXX([$CXX_TRY_LIST])
  # --------------------------------------------------------------



  # --------------------------------------------------------------
  # See below for the definition of this function.  It can
  # figure out which version of a particular compiler, e.g. GCC 4.0,
  # you are using.
  # --------------------------------------------------------------
  DETERMINE_CXX_BRAND



  # --------------------------------------------------------------
  # look for a decent C compiler or honor --with-cc=...
  CC_TRY_LIST="gcc icc pgcc cc"
  if  (test "$enablempi" != no) ; then
    CC_TRY_LIST="mpicc $CC_TRY_LIST"
  fi
  AC_ARG_WITH([cc],
              AS_HELP_STRING([--with-cc=CC], [C compiler to use]),
              [CC="$withval"],
              [])

  # --------------------------------------------------------------
  # Determine a C compiler to use.  If CC is not already set, checks for
  # gcc, cc, and other C compilers.  Then sets the CC variable to the result.
  # --------------------------------------------------------------
  AC_PROG_CC([$CC_TRY_LIST])
  # --------------------------------------------------------------


  # --------------------------------------------------------------
  # libMesh itself is not written in any Fortran and does not need
  # a Fortran compiler.  Many optional packages however are and
  # we need the compiler to figure out how to link those libraries
  #
  # note though than on OSX for example the XCode tools provide
  # a 'mpif77' which will be detected below but is actually an
  # empty shell script wrapper.  Then the compiler will fail to
  # make executables and we will wind up with static libraries
  # due to a bizarre chain of events.  So, add support for
  # --disable-fortran
  # --------------------------------------------------------------
  AC_ARG_ENABLE(fortran,
                AS_HELP_STRING([--disable-fortran],
                               [build without Fortran language support]),
                [case "${enableval}" in
                  yes)  enablefortran=yes ;;
                  no)  enablefortran=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-fortran) ;;
                esac],
                [enablefortran=yes])

  if (test "x$enablefortran" = xyes); then

    # look for a decent F90+ compiler or honor --with-fc=...
    FC_TRY_LIST="gfortran ifort pgf90 xlf95"
    if  (test "$enablempi" != no) ; then
      FC_TRY_LIST="mpif90 $FC_TRY_LIST"
    fi
    AC_ARG_WITH([fc],
                AS_HELP_STRING([--with-fc=FC], [Fortran compiler to use]),
                [FC="$withval"],
                [])

    # --------------------------------------------------------------
    # Determine a F90+ compiler to use.
    # --------------------------------------------------------------
    AC_PROG_FC([$FC_TRY_LIST])

    if (test "x$FC" = "x"); then
      AC_MSG_RESULT(>>> No valid Fortran compiler <<<)
      FC=no
      enablefortran=no
    fi
    # --------------------------------------------------------------



    # --------------------------------------------------------------
    # look for a decent F77 compiler or honor --with-77=...
    F77_TRY_LIST="gfortran g77 ifort f77 xlf frt pgf77 fort77 fl32 af77 f90 xlf90 pgf90 epcf90 f95 fort xlf95 ifc efc pgf95 lf95"
    if  (test "$enablempi" != no) ; then
      F77_TRY_LIST="mpif77 $F77_TRY_LIST"
    fi
    AC_ARG_WITH([f77],
                AS_HELP_STRING([--with-f77=F77], [Fortran compiler to use]),
                [F77="$withval"],
                [])

    # --------------------------------------------------------------
    # Determine a F77 compiler to use.
    # --------------------------------------------------------------
    AC_PROG_F77([$F77_TRY_LIST])

    if (test "x$F77" = "x"); then
      AC_MSG_RESULT(>>> No valid Fortran 77 compiler <<<)
      F77=no
      enablefortran=no
    fi

    # --------------------------------------------------------------
  else
      # when --disable-fortran is specified, explicitly set these
      # to "no" to instruct libtool not to bother with them.
      AC_MSG_RESULT(>>> Disabling Fortran language support per user request <<<)
      FC=no
      F77=no
  fi # end enablefortran
])



dnl -------------------------------------------------------------
dnl Determine the C++ compiler in use. Return the name and possibly
dnl version of this compiler in GXX_VERSION.
dnl
dnl Usage: DETERMINE_CXX_BRAND
dnl
dnl -------------------------------------------------------------
AC_DEFUN([DETERMINE_CXX_BRAND],
[
  dnl Set this flag once the compiler "brand" has been detected to skip remaining tests.
  compiler_brand_detected=no

  dnl First check for gcc version, prevents intel's icc from
  dnl pretending to be gcc
  REAL_GXX=`($CXX -v 2>&1) | grep "gcc version"`

  dnl Do not allow Intel to masquerade as a "real" GCC.
  is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'Compiler'`"
  if test "x$is_intel_icc" != "x" ; then
    REAL_GXX=""
  fi

  if (test "$GXX" = yes -a "x$REAL_GXX" != "x" ) ; then
    dnl find out the right version
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "gcc version"`
    case "$GXX_VERSION_STRING" in
      *gcc\ version\ 7.*)
        AC_MSG_RESULT(<<< C++ compiler is gcc-7.x >>>)
        GXX_VERSION=gcc7
        ;;
      *gcc\ version\ 6.*)
        AC_MSG_RESULT(<<< C++ compiler is gcc-6.x >>>)
        GXX_VERSION=gcc6
        ;;
      *gcc\ version\ 5.*)
        AC_MSG_RESULT(<<< C++ compiler is gcc-5.x >>>)
        GXX_VERSION=gcc5
        ;;
      *4.9.*)
        AC_MSG_RESULT(<<< C++ compiler is gcc-4.9 >>>)
        GXX_VERSION=gcc4.9
        ;;
      *4.8.*)
        AC_MSG_RESULT(<<< C++ compiler is gcc-4.8 >>>)
        GXX_VERSION=gcc4.8
        ;;
      *4.7.*)
        AC_MSG_RESULT(<<< C++ compiler is gcc-4.7 >>>)
        GXX_VERSION=gcc4.7
        ;;
      *4.6.*)
        AC_MSG_RESULT(<<< C++ compiler is gcc-4.6 >>>)
        GXX_VERSION=gcc4.6
        ;;
      *)
        AC_MSG_RESULT(<<< C++ compiler is unknown but accepted gcc version >>>)
        GXX_VERSION=gcc-other
        ;;
    esac

    dnl Detection was successful, so set the flag.
    compiler_brand_detected=yes
  fi

  dnl Clang/LLVM C++?
  if (test "x$compiler_brand_detected" = "xno"); then
    clang_version="`($CXX --version 2>&1)`"
    is_clang="`echo $clang_version | grep 'clang'`"

    if test "x$is_clang" != "x" ; then
      dnl Detect if clang is the version built by
      dnl Apple, because then the version number means
      dnl something different...
      is_apple_clang="`echo $clang_version | grep 'Apple'`"
      clang_vendor="clang"
      if test "x$is_apple_clang" != "x" ; then
        clang_vendor="Apple clang"
      fi

      dnl If we have perl, we can also pull out the clang version number using regexes.
      dnl Note that @S|@ is a quadrigraph for the dollar sign.
      clang_major_minor=unknown

      if test "x$PERL" != "x" ; then
         clang_major_minor=`echo $clang_version | $PERL -ne 'print @S|@1 if /version\s(\d+\.\d+)/'`
         if test "x$clang_major_minor" = "x" ; then
           clang_major_minor=unknown
         fi
      fi

      AC_MSG_RESULT([<<< C++ compiler is ${clang_vendor}, version ${clang_major_minor} >>>])
      GXX_VERSION=clang
      compiler_brand_detected=yes
    fi
  fi

  dnl Intel's ICC C++ compiler?
  if (test "x$compiler_brand_detected" = "xno"); then
    is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'Compiler'`"
    if test "x$is_intel_icc" != "x" ; then
      GXX_VERSION_STRING="`($CXX -V 2>&1) | grep 'Version '`"
      case "$GXX_VERSION_STRING" in
        *18.*)
          AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 18 >>>)
          GXX_VERSION=intel_icc_v18.x
          ;;
        *17.*)
          AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 17 >>>)
          GXX_VERSION=intel_icc_v17.x
          ;;
        *16.*)
          AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 16 >>>)
          GXX_VERSION=intel_icc_v16.x
          ;;
        *15.*)
          AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 15 >>>)
          GXX_VERSION=intel_icc_v15.x
          ;;
        *14.*)
          AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 14 >>>)
          GXX_VERSION=intel_icc_v14.x
          ;;
        *13.*)
          AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 13 >>>)
          GXX_VERSION=intel_icc_v13.x
          ;;
      esac
      compiler_brand_detected=yes
    fi
  fi

  dnl Check for IBM xlC. Depending on environment
  dnl variables, moon position, and other reasons unknown to me, this
  dnl compiler can display different names in the first line of output, so
  dnl check various possibilities.  Calling xlC with no arguments displays
  dnl the man page.  Grepping for case-sensitive xlc is not enough if the
  dnl user wants xlC, so we use case-insensitive grep instead.
  if (test "x$compiler_brand_detected" = "xno"); then
    is_ibm_xlc="`($CXX 2>&1) | egrep -i 'xlc'`"
    if test "x$is_ibm_xlc" != "x"  ; then
      AC_MSG_RESULT(<<< C++ compiler is IBM xlC >>>)
      GXX_VERSION=ibm_xlc
      compiler_brand_detected=yes
    fi
  fi

  dnl Cray C++?
  if (test "x$compiler_brand_detected" = "xno"); then
    is_cray_cc="`($CXX -V 2>&1) | grep 'Cray '`"
    if test "x$is_cray_cc" != "x" ; then
      AC_MSG_RESULT(<<< C++ compiler is Cray C++ >>>)
      GXX_VERSION=cray_cc
      compiler_brand_detected=yes
    fi
  fi

  dnl Portland Group C++?
  if (test "x$compiler_brand_detected" = "xno"); then
    is_pgcc="`($CXX -V 2>&1) | grep 'Portland Group'`"
    if test "x$is_pgcc" != "x" ; then
      AC_MSG_RESULT(<<< C++ compiler is Portland Group C++ >>>)
      GXX_VERSION=portland_group
      compiler_brand_detected=yes
    fi
  fi

  dnl Could not recognize the compiler. Warn the user and continue.
  if (test "x$compiler_brand_detected" = "xno"); then
    AC_MSG_RESULT( WARNING:)
    AC_MSG_RESULT( >>> Unrecognized compiler: "$CXX" <<<)
    AC_MSG_RESULT( You will likely need to modify)
    AC_MSG_RESULT( Make.common directly to specify)
    AC_MSG_RESULT( proper compiler flags)
    GXX_VERSION=unknown
  fi
])





# -------------------------------------------------------------
# Set C++ compiler flags to their default values. They will be
# modified according to other options in later steps of
# configuration
#
# CXXFLAGS_OPT    : flags for optimized mode
# CXXFLAGS_DEVEL  : flags for development mode
# CXXFLAGS_DBG    : flags for debug mode
# CPPFLAGS_OPT    : preprocessor flags for optimized mode
# CPPFLAGS_DEVEL  : preprocessor flags for development mode
# CPPFLAGS_DBG    : preprocessor flags for debug mode
# PROFILING_FLAGS : flags to enable code profiling
# ASSEMBLY_FLAGS  : flags to enable assembly language output
#
# Usage: SET_CXX_FLAGS
#
# (Note the CXXFLAGS and the CPPFLAGS used for further tests may
#  be augmented)
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_SET_CXX_FLAGS],
[
  # method-specific preprocessor flags, independent of compiler.
  CPPFLAGS_OPT="-DNDEBUG"
  CPPFLAGS_DBG="-DDEBUG"
  CPPFLAGS_DEVEL=""

  # Flag to add directories to the dynamic library search path; can
  # be changed at a later stage
  RPATHFLAG="-Wl,-rpath,"

  # Flag for profiling mode; can me modified at a later stage
  PROFILING_FLAGS="-pg"

  # Flag for assembly-output mode; can me modified at a later stage
  ASSEMBLY_FLAGS="-S"

  # The -g flag is necessary for OProfile to produce annotations
  # -fno-omit-frame-pointer flag turns off an optimization that
  # interferes with OProfile callgraphs
  OPROFILE_FLAGS="-g -fno-omit-frame-pointer"

  # For compilers that support it (clang >= 3.5.0 and GCC >= 4.8) the
  # user can selectively enable "sanitize" flags for different METHODs
  # with e.g.
  #
  # --enable-sanitize="dbg opt"
  #
  # These flags generally slow down code execution, so you don't
  # necessarily want to turn them on all the time or in all METHODs.

  # Declaring something AC_ARG_VAR does several things, see
  # http://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/Setting-Output-Variables.html
  # for more information... not sure we need it in this case.
  # AC_ARG_VAR([SANITIZE_METHODS], [methods we apply sanitizer flags to, e.g. "opt dbg devel"])

  AC_ARG_ENABLE([sanitize],
              AS_HELP_STRING([--enable-sanitize="opt dbg devel prof oprof"],
                             [turn on sanitizer flags for the given methods]),
              [for method in ${enableval} ; do
                 # make sure each method specified makes sense
                 case "${method}" in
                     optimized|opt)      ;;
                     debug|dbg)          ;;
                     devel)              ;;
                     profiling|pro|prof) ;;
                     oprofile|oprof)     ;;
                     *)
                         AC_MSG_ERROR(bad value ${method} for --enable-sanitize)
                         ;;
                 esac
               done
               # If we made it here, the case statement didn't detect any errors
               SANITIZE_METHODS=${enableval}],
               [])

  if test "x$SANITIZE_METHODS" != x; then
    AC_MSG_RESULT([<<< Testing sanitizer flags for method(s) "$SANITIZE_METHODS" >>>])

    # Both Clang and GCC docs suggest using "-fsanitize=address -fno-omit-frame-pointer".
    # The Clang documentation further suggests using "-O1 -g -fno-optimize-sibling-calls".
    # Since these flags also work in GCC, we'll use them there as well...
    COMMON_SANITIZE_OPTIONS="-fsanitize=address -fno-omit-frame-pointer -O1 -g -fno-optimize-sibling-calls"

    # Test the sanitizer flags.  Currently Clang and GCC are the only
    # compilers that support the address sanitizer, and they use the
    # same set of flags.  If the set of flags used by Clang and GCC ever
    # diverges, we'll need to set up separate flags and test them in the
    # case blocks below...  The LIBMESH_TEST_SANITIZE_FLAGS function sets
    # the variable have_address_sanitizer to either "no" or "yes"
    LIBMESH_TEST_SANITIZE_FLAGS([$COMMON_SANITIZE_OPTIONS])

    # Enable the address sanitizer stuff if the test code compiled
    if test "x$have_address_sanitizer" = xyes; then
      # As of clang 3.9.0 or so, we also need to pass the sanitize flag to the linker
      # if it's being used during compiling. It seems that we do not need to pass all
      # the flags above, just the sanitize flag itself.
      libmesh_LDFLAGS="$libmesh_LDFLAGS -Wc,-fsanitize=address"

      for method in ${SANITIZE_METHODS}; do
          case "${method}" in
              optimized|opt)
                SANITIZE_OPT_FLAGS=$COMMON_SANITIZE_OPTIONS
                ;;

              debug|dbg)
                SANITIZE_DBG_FLAGS=$COMMON_SANITIZE_OPTIONS
                ;;

              devel)
                SANITIZE_DEVEL_FLAGS=$COMMON_SANITIZE_OPTIONS
                ;;

              profiling|pro|prof)
                SANITIZE_PROF_FLAGS=$COMMON_SANITIZE_OPTIONS
                ;;

              oprofile|oprof)
                SANITIZE_OPROF_FLAGS=$COMMON_SANITIZE_OPTIONS
                ;;

              *)
                AC_MSG_ERROR(bad value ${method} for --enable-sanitize)
                ;;
          esac
      done
    fi
  fi

  # in the case blocks below we may add GLIBCXX-specific pedantic debugging preprocessor
  # definitions. however, allow the knowing user to preclude that if they need to.
  AC_ARG_ENABLE(glibcxx-debugging,
                [AS_HELP_STRING([--disable-glibcxx-debugging],
                                [omit -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC even in dbg mode])],
                [case "${enableval}" in
                  yes)  enableglibcxxdebugging=yes ;;
                  no)  enableglibcxxdebugging=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-glibcxx-debugging) ;;
                esac],
                [enableglibcxxdebugging=yes])

  # GLIBCXX debugging causes untold woes on mac machines - so disable it
  if (test `uname` = "Darwin"); then
    AC_MSG_RESULT(<<< Disabling GLIBCXX debugging on Darwin >>>)
    enableglibcxxdebugging=no
  fi
  AM_CONDITIONAL(LIBMESH_ENABLE_GLIBCXX_DEBUGGING, test x$enableglibcxxdebugging = xyes)


  # GLIBCXX-specific debugging flags conflict with cppunit on many of
  # our users' systems.  However, being able to override this allows
  # us to increase our unit test coverage.
  AC_ARG_ENABLE(glibcxx-debugging-cppunit,
	 [AS_HELP_STRING([--enable-glibcxx-debugging-cppunit],
	                 [Use GLIBCXX debugging flags for unit tests])],
	 [case "${enableval}" in
	   yes)  enableglibcxxdebuggingcppunit=yes ;;
	    no)  enableglibcxxdebuggingcppunit=no ;;
 	     *)  AC_MSG_ERROR(bad value ${enableval} for --enable-glibcxx-debugging-cppunit) ;;
	  esac],
	[enableglibcxxdebuggingcppunit=no])

  AM_CONDITIONAL(LIBMESH_ENABLE_GLIBCXX_DEBUGGING_CPPUNIT, test x$enableglibcxxdebuggingcppunit = xyes)


  # First the flags for gcc compilers
  if (test "$GXX" = yes -a "x$REAL_GXX" != "x" ) ; then
    CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization"
    CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -O2 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -Wuninitialized -funroll-loops -fstrict-aliasing -Woverloaded-virtual -Wdisabled-optimization"
    CXXFLAGS_DBG="$CXXFLAGS_DBG -O0 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -Woverloaded-virtual"
    NODEPRECATEDFLAG="-Wno-deprecated"

    CFLAGS_OPT="-O2 -funroll-loops -fstrict-aliasing"
    CFLAGS_DEVEL="$CFLAGS_OPT -g -Wimplicit -funroll-loops -fstrict-aliasing"
    CFLAGS_DBG="-g -Wimplicit"
    ASSEMBLY_FLAGS="$ASSEMBLY_FLAGS -fverbose-asm"

    if (test "x$enableglibcxxdebugging" = "xyes"); then
      CPPFLAGS_DBG="$CPPFLAGS_DBG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"
    fi

    # GCC 4.6.3 warns about variadic macros but supports them just
    # fine, so let's turn off that warning.
    case "$GXX_VERSION" in
      gcc4.6 | gcc5)
        CXXFLAGS_OPT="$CXXFLAGS_OPT -Wno-variadic-macros"
        CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -Wno-variadic-macros"
        CXXFLAGS_DBG="$CXXFLAGS_DBG -Wno-variadic-macros"
        ;;

      *)
        ;;
    esac


    # Set OS-specific flags for linkers & other stuff
    case "$target" in

      # For Solaris we need to pass a different flag to the linker for specifying the
      # dynamic library search path and add -lrpcsvc to use XDR
      *solaris*)
          RPATHFLAG="-Wl,-R,"
          LIBS="-lrpcsvc $LIBS"
          ;;

      *)
          ;;
    esac


  else
    # Non-gcc compilers

    case "$GXX_VERSION" in
      ibm_xlc)
          CXXFLAGS_OPT="-O3 -qmaxmem=-1 -w -qansialias -Q=10 -qrtti=all -qstaticinline"
          CXXFLAGS_DBG="-qmaxmem=-1 -qansialias -qrtti=all -g -qstaticinline"
          CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
          NODEPRECATEDFLAG=""
          CFLAGS_OPT="-O3 -qmaxmem=-1 -w -qansialias -Q=10"
          CFLAGS_DBG="-qansialias -g"
          CFLAGS_DEVEL="$CFLAGS_DBG"
          ;;

      # All Intel ICC/ECC flavors
      intel_*)

        # Intel understands the gcc-like no-deprecated flag
        NODEPRECATEDFLAG="-Wno-deprecated"

        # Intel compilers use -qp for profiling
        PROFILING_FLAGS="-qp"

        # Intel options for annotated assembly
        ASSEMBLY_FLAGS="$ASSEMBLY_FLAGS -fverbose-asm -fsource-asm"

        # The -g flag is all OProfile needs to produce annotations
        OPROFILE_FLAGS="-g"

        # Specific flags for specific versions
        case "$GXX_VERSION" in

          # Intel ICC >= v13.x
          intel_icc_v13.x | intel_icc_v14.x | intel_icc_v15.x | intel_icc_v16.x | intel_icc_v17.x | intel_icc_v18.x)
              # Disable some warning messages:
              # #161: 'unrecognized #pragma
              #       #pragma GCC diagnostic warning "-Wdeprecated-declarations"'
              #       I don't understand this, pragmas for other compilers should
              #       just be silently ignored, isn't that the whole point of pragmas?
              # #175: 'subscript out of range'
              #       FIN-S application code causes many false
              #       positives with this
              # #266: 'function declared implicitly'
              #       Metis function "GKfree" caused this error
              #       in almost every file.
              # #488: 'template parameter "Scalar1" is not used in declaring the
              #       parameter types of function template'
              #       This warning was generated from one of the type_vector.h
              #       constructors that uses some SFINAE tricks.
              # #1476: 'field uses tail padding of a base class'
              # #1505: 'size of class is affected by tail padding'
              #        simply warns of a possible incompatibility with
              #        the g++ ABI for this case
              # #1572: 'floating-point equality and inequality comparisons are unreliable'
              #        Well, duh, when the tested value is computed...  OK when it
              #        was from an assignment.
              PROFILING_FLAGS="-p"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -w1 -g -wd175 -wd1476 -wd1505 -wd1572 -wd488 -wd161"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -O3 -unroll -w0 -ftz"
              CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -w1 -g -wd175 -wd1476 -wd1505 -wd1572 -wd488 -wd161"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -g -wd266 -wd1572 -wd488 -wd161"
              CFLAGS_OPT="$CFLAGS_OPT -O3 -unroll -w0 -ftz"
              CFLAGS_DEVEL="$CFLAGS_DBG"
              ;;

          *)
              AC_MSG_RESULT(Unknown Intel compiler, "$GXX_VERSION")
              ;;
        esac
      ;;

      portland_group)
          CXXFLAGS_DBG="-g --no_using_std"
          CXXFLAGS_OPT="-O2 --no_using_std -fast -Minform=severe"
          CXXFLAGS_DEVEL="$CXXFLAGS_DBG"

          # PG C++ definitely doesn't understand -Wno-deprecated...
          NODEPRECATEDFLAG=""

          CFLAGS_DBG="-g"
          CFLAGS_OPT="-O2"
          CFLAGS_DEVEL="$CFLAGS_DBG"

          # Disable exception handling if we don't use it
          if test "$enableexceptions" = no ; then
            CXXFLAGS_DBG="$CXXFLAGS_DBG --no_exceptions"
            CXXFLAGS_OPT="$CXXFLAGS_OPT --no_exceptions"
          fi
          ;;

      cray_cc)
          CXXFLAGS_DBG="-h conform,one_instantiation_per_object,instantiate=used,noimplicitinclude -G n"
          CXXFLAGS_OPT="-h conform,one_instantiation_per_object,instantiate=used,noimplicitinclude -G n"
          CXXFLAGS_DEVEL="-h conform,one_instantiation_per_object,instantiate=used,noimplicitinclude -G n"
          NODEPRECATEDFLAG=""
          CFLAGS_DBG="-G n"
          CFLAGS_OPT="-G n"
          CFLAGS_DEVEL="-G n"
          ;;

      clang)
          CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -felide-constructors -Qunused-arguments -Wunused-parameter -Wunused"
          CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -O2 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long -Wunused-parameter -Wunused -Wpointer-arith -Wformat -Wparentheses -Wuninitialized -Qunused-arguments -Woverloaded-virtual -fno-limit-debug-info"
          CXXFLAGS_DBG="$CXXFLAGS_DBG -O0 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long -Wunused-parameter -Wunused -Wpointer-arith -Wformat -Wparentheses -Qunused-arguments -Woverloaded-virtual -fno-limit-debug-info"
          NODEPRECATEDFLAG="-Wno-deprecated"

          CFLAGS_OPT="-O2 -Qunused-arguments -Wunused"
          CFLAGS_DEVEL="$CFLAGS_OPT -g -Wimplicit -fno-limit-debug-info -Wunused"
          CFLAGS_DBG="-g -Wimplicit -Qunused-arguments -fno-limit-debug-info -Wunused"
          ;;

      *)
          AC_MSG_RESULT(No specific options for this C++ compiler known)
          CXXFLAGS_DBG="$CXXFLAGS"
          CXXFLAGS_OPT="$CXXFLAGS"
          CXXFLAGS_DEVEL="$CXXFLAGS"
          NODEPRECATEDFLAG=""

          CFLAGS_DBG="$CFLAGS"
          CFLAGS_OPT="$CFLAGS"
          CFLAGS_DEVEL="$CFLAGS"
          ;;
    esac
  fi
])
