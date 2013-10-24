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
                   AC_HELP_STRING([--enable-mpi],
                                  [build with MPI message passing support]),
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
  	    AC_HELP_STRING([--with-cxx=CXX],
                             [C++ compiler to use]),
  	    [CXX="$withval"],
  	    [])

  # --------------------------------------------------------------
  # Determines a C++ compiler to use.  First checks if the variable CXX is
  # already set.  If not, then searches under g++, c++, and other names.
  # --------------------------------------------------------------
  AC_PROG_CXX([$CXX_TRY_LIST])
  # --------------------------------------------------------------



  # --------------------------------------------------------------
  # See aclocal.m4 for the definition of this function.  It can
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
  	    AC_HELP_STRING([--with-cc=CC],
                             [C compiler to use]),
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
  # emtpy shell script wrapper.  Then the compiler will fail to
  # make executables and we will wind up with static libraries
  # due to a bizarre chain of events.  So, add support for
  # --disable-fortran
  # --------------------------------------------------------------
  AC_ARG_ENABLE(fortran,
                AC_HELP_STRING([--enable-fortran],
                               [build with Fortran language support]),
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
    	    AC_HELP_STRING([--with-fc=FC],
                               [Fortran compiler to use]),
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
    	    AC_HELP_STRING([--with-f77=F77],
                               [Fortran compiler to use]),
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



# -------------------------------------------------------------
# Determine the C++ compiler in use. Return the name and possibly
# version of this compiler in GXX_VERSION.
#
# Usage: DETERMINE_CXX_BRAND
#
# -------------------------------------------------------------
AC_DEFUN([DETERMINE_CXX_BRAND],
[
  # First check for gcc version, avoids intel's icc from
  # pretending to be gcc
  REAL_GXX=`($CXX -v 2>&1) | grep "gcc version"`

  # Intel's v12.1 does this:
  # $ icpc -v
  #   icpc version 12.1.0 (gcc version 4.4.4 compatibility)
  # cath that and do not interpret it as 'REAL_GXX' compiler
  is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'Compiler'`"
  if test "x$is_intel_icc" != "x" ; then
    REAL_GXX=""
  fi

  if (test "$GXX" = yes -a "x$REAL_GXX" != "x" ) ; then
    # find out the right version
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "gcc version"`
    case "$GXX_VERSION_STRING" in
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
      *4.5.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-4.5 >>>)
  	GXX_VERSION=gcc4.5
  	;;
      *4.4.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-4.4 >>>)
  	GXX_VERSION=gcc4.4
  	;;
      *4.3.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-4.3 >>>)
  	GXX_VERSION=gcc4.3
  	;;
      *4.2.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-4.2 >>>)
  	GXX_VERSION=gcc4.2
  	;;
      *4.1.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-4.1 >>>)
  	GXX_VERSION=gcc4.1
  	;;
      *4.0.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-4.0 >>>)
  	GXX_VERSION=gcc4.0
  	;;
      *3.4.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.4 >>>)
  	GXX_VERSION=gcc3.4
  	;;
      *3.3.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.3 >>>)
  	GXX_VERSION=gcc3.3
  	;;
      *3.2.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.2 >>>)
  	GXX_VERSION=gcc3.2
  	;;
      *3.1.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.1 >>>)
  	GXX_VERSION=gcc3.1
  	;;
      *3.0.*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.0 >>>)
  	GXX_VERSION=gcc3.0
  	;;
      *2.97*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-2.97 >>>)
  	GXX_VERSION=gcc2.97
  	;;
      *2.96*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-2.96 >>>)
  	GXX_VERSION=gcc2.96
	AC_DEFINE(BROKEN_IOSTREAM, 1,
             [This compiler is known not to support some iostream
              functionality])
         AC_MSG_RESULT(<<< Configuring library for broken iostream >>>)
  	;;
      *2.95*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-2.95 >>>)
  	GXX_VERSION=gcc2.95
	AC_DEFINE(BROKEN_IOSTREAM, 1,
             [This compiler is known not to support some iostream
              functionality])
         AC_MSG_RESULT(<<< Configuring library for broken iostream >>>)
  	;;
      *"egcs-1.1"*)
  	AC_MSG_RESULT(<<< C++ compiler is egcs-1.1 >>>)
  	GXX_VERSION=egcs1.1
  	;;
      *2.4* | *2.5* | *2.6* | *2.7* | *2.8*)
  	# These compilers are too old to support a useful subset
  	# of modern C++, so we don't support them
  	AC_MSG_RESULT(<<< C++ compiler is $GXX_VERSION_STRING >>>)
  	AC_MSG_ERROR(<<< C++ compiler is not supported >>>)
  	;;
      *)
  	AC_MSG_RESULT(<<< C++ compiler is unknown but accepted gcc version >>>)
  	GXX_VERSION=gcc-other
  	;;
    esac
    # Check for Apple compilers
    case "$GXX_VERSION_STRING" in
      *Apple*)
        AC_MSG_RESULT(<<< C++ compiler is built by Apple >>>)
        APPLE_GCC=true
        ;;
      *)
        APPLE_GCC=false
        ;;
    esac
  else
    # Check other (non-gcc) compilers

    # Check for IBM xlC. For some reasons, depending on some environment
    # variables, moon position, and other reasons unknown to me, the
    # compiler displays different names in the first line of output, so
    # check various possibilities.  Calling xlC with no arguments displays
    # the man page.  Grepping for case-sensitive xlc is not enough if the
    # user wants xlC, so we used case-insensitive grep instead.
    is_ibm_xlc="`($CXX 2>&1) | egrep -i 'xlc'`"
    if test "x$is_ibm_xlc" != "x"  ; then
      # IBM's C++ compiler.
      AC_MSG_RESULT(<<< C++ compiler is IBM xlC >>>)
      GXX_VERSION=ibm_xlc
    else

      # Check whether we are dealing with the MIPSpro C++ compiler
      is_mips_pro="`($CXX -version 2>&1) | grep MIPSpro`"
      if test "x$is_mips_pro" != "x" ; then
        AC_MSG_RESULT(<<< C++ compiler is MIPSpro C++ compiler >>>)
        GXX_VERSION=MIPSpro
      else

        # Intel's ICC C++ compiler for Itanium?
        is_intel_ecc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'Itanium(R)' | grep 'Compiler'`"
        if test "x$is_intel_ecc" = "x" ; then
          is_intel_ecc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'IA-64' | grep 'Compiler'`"
        fi
        if test "x$is_intel_ecc" != "x" ; then
          GXX_VERSION_STRING="`($CXX -V -help 2>&1) | grep 'Version '`"
          case "$GXX_VERSION_STRING" in
            *10.1*)
              AC_MSG_RESULT(<<< C++ compiler is Intel Itanium ICC 10.1 >>>)
  	      GXX_VERSION=intel_itanium_icc_v10.1
              ;;
            *10.0*)
              AC_MSG_RESULT(<<< C++ compiler is Intel Itanium ICC 10.0 >>>)
  	      GXX_VERSION=intel_itanium_icc_v10.0
              ;;
            *9.1*)
              AC_MSG_RESULT(<<< C++ compiler is Intel Itanium ICC 9.1 >>>)
  	      GXX_VERSION=intel_itanium_icc_v9.1
              ;;
            *9.0*)
              AC_MSG_RESULT(<<< C++ compiler is Intel Itanium ICC 9.0 >>>)
  	      GXX_VERSION=intel_itanium_icc_v9.0
              ;;
            *8.1*)
              AC_MSG_RESULT(<<< C++ compiler is Intel Itanium ICC 8.1 >>>)
  	      GXX_VERSION=intel_itanium_icc_v8.1
              ;;
            *8.0*)
              AC_MSG_RESULT(<<< C++ compiler is Intel Itanium ICC 8.0 >>>)
  	      GXX_VERSION=intel_itanium_icc_v8.0
              ;;
            *7.1*)
              AC_MSG_RESULT(<<< C++ compiler is Intel Itanium ICC 7.1 >>>)
  	      GXX_VERSION=intel_itanium_icc_v7.1
              ;;
            *7.0*)
              AC_MSG_RESULT(<<< C++ compiler is Intel Itanium ICC 7.0 >>>)
  	      GXX_VERSION=intel_itanium_icc_v7.0
              ;;
          esac
        else

          # Intel's ICC C++ compiler?
          is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'Compiler'`"
          if test "x$is_intel_icc" != "x" ; then
            GXX_VERSION_STRING="`($CXX -V 2>&1) | grep 'Version '`"
            case "$GXX_VERSION_STRING" in
              *14.*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 14 >>>)
                GXX_VERSION=intel_icc_v14.x
                ;;
              *13.*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 13 >>>)
                GXX_VERSION=intel_icc_v13.x
                ;;
              *12.1*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 12.1 >>>)
  	        GXX_VERSION=intel_icc_v12.x
                ;;
              *12.*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 12 >>>)
  	        GXX_VERSION=intel_icc_v12.x
                ;;
              *11.*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 11 >>>)
  	        GXX_VERSION=intel_icc_v11.x
                ;;
              *10.1*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 10.1 >>>)
  	        GXX_VERSION=intel_icc_v10.1
                ;;
              *10.0*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 10.0 >>>)
  	        GXX_VERSION=intel_icc_v10.0
                ;;
              *9.1*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 9.1 >>>)
  	        GXX_VERSION=intel_icc_v9.1
                ;;
              *9.0*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 9.0 >>>)
  	        GXX_VERSION=intel_icc_v9.0
                ;;
              *8.1*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 8.1 >>>)
  	        GXX_VERSION=intel_icc_v8.1
                ;;
              *8.0*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 8.0 >>>)
  	        GXX_VERSION=intel_icc_v8.0
                ;;
              *7.1*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 7.1 >>>)
  	        GXX_VERSION=intel_icc_v7.1
                ;;
              *7.0*)
                AC_MSG_RESULT(<<< C++ compiler is Intel(R) icc 7.0 >>>)
  	        GXX_VERSION=intel_icc_v7.0
                ;;
            esac
          else

            # Or Compaq's cxx compiler?
            is_dec_cxx="`($CXX -V 2>&1) | grep 'Compaq C++'`"
            if test "x$is_dec_cxx" != "x" ; then
              AC_MSG_RESULT(<<< C++ compiler is Compaq cxx >>>)
              GXX_VERSION=compaq_cxx
            else

  	      # Sun Studio?
              is_sun_cc="`($CXX -V 2>&1) | grep 'Sun C++'`"
              if test "x$is_sun_cc" != "x" ; then
                AC_MSG_RESULT(<<< C++ compiler is Sun Studio compiler >>>)
                GXX_VERSION=sun_studio
              else

  	        # Sun Forte?
                is_sun_forte_cc="`($CXX -V 2>&1) | grep 'Forte'`"
                if test "x$is_sun_forte_cc" != "x" ; then
                  AC_MSG_RESULT(<<< C++ compiler is Sun Forte compiler >>>)
                  GXX_VERSION=sun_forte
                else

  	          # Cray C++?
  	          is_cray_cc="`($CXX -V 2>&1) | grep 'Cray '`"
  	          if test "x$is_cray_cc" != "x" ; then
  	            AC_MSG_RESULT(<<< C++ compiler is Cray C++ >>>)
  	            GXX_VERSION=cray_cc
  	          else

  	            # Portland Group C++?
  	            is_pgcc="`($CXX -V 2>&1) | grep 'Portland Group'`"
  	            if test "x$is_pgcc" != "x" ; then
  	              AC_MSG_RESULT(<<< C++ compiler is Portland Group C++ >>>)
  	              GXX_VERSION=portland_group
  	            else

                      # HP-UX 11.11 aCC?
                      is_hpux_acc="`($CXX -V 2>&1) | grep 'aCC: HP ANSI C++'`"
  	              if test "x$is_hpux_acc" != "x" ; then
  	                AC_MSG_RESULT(<<< C++ compiler is HP-UX C++ >>>)
    	                GXX_VERSION=hpux_acc
  	              else

		        # Clang LLVM C++?
			is_clang="`($CXX --version 2>&1) | grep 'clang'`"
  	                if test "x$is_clang" != "x" ; then
                          AC_MSG_RESULT(<<< C++ compiler is LLVM Clang C++ >>>)
	                  GXX_VERSION=clang
			else

                          # No recognized compiler found...
			  # warn the user and continue
			  AC_MSG_RESULT( WARNING:)
                          AC_MSG_RESULT( >>> Unrecognized compiler: "$CXX" <<<)
			  AC_MSG_RESULT( You will likely need to modify)
			  AC_MSG_RESULT( Make.common directly to specify)
			  AC_MSG_RESULT( proper compiler flags)
			  GXX_VERSION=unknown
			fi
                      fi
                    fi
                  fi
                fi
              fi
  	    fi
          fi
        fi
      fi
    fi
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



  # in the case blocks below we may add GLIBCXX-specific pedantic debugging preprocessor
  # definitions. however, allow the knowing user to preclude that if they need to.
  AC_ARG_ENABLE(glibcxx-debugging,
	 [AC_HELP_STRING([--enable-glibcxx-debugging],
	                 [add -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC to CXXFLAGS_DBG (yes by default)])],
	 [case "${enableval}" in
	   yes)  enableglibcxxdebugging=yes ;;
	    no)  enableglibcxxdebugging=no ;;
 	     *)  AC_MSG_ERROR(bad value ${enableval} for --enable-glibcxx-debugging) ;;
	  esac],
	[enableglibcxxdebugging=yes])

  # GLIBCXX debugging causes untold woes on mac machines - so disable it
  if (test `uname` = "Darwin"); then
    AC_MSG_RESULT(<< Disabling GLIBCXX debugging on Darwin >>>)
    enableglibcxxdebugging=no
  fi
  AM_CONDITIONAL(LIBMESH_ENABLE_GLIBCXX_DEBUGGING, test x$enableglibcxxdebugging = xyes)


  # First the flags for gcc compilers
  if (test "$GXX" = yes -a "x$REAL_GXX" != "x" ) ; then
    CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -felide-constructors"
    CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -O2 -felide-constructors -g -ansi -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -Wuninitialized"
    CXXFLAGS_DBG="$CXXFLAGS_DBG -O0 -felide-constructors -g -ansi -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses"
    NODEPRECATEDFLAG="-Wno-deprecated"

    CFLAGS_OPT="-O2"
    CFLAGS_DEVEL="$CFLAGS_OPT -g -Wimplicit"
    CFLAGS_DBG="-g -Wimplicit"
    ASSEMBLY_FLAGS="$ASSEMBLY_FLAGS -fverbose-asm"

    # set some flags that are specific to some versions of the
    # compiler:
    # - egcs1.1 yielded incorrect code with some loop unrolling
    # - after egcs1.1, the optimization flag -fstrict-aliasing was
    #   introduced, which enables better optimizations for
    #   well-written C++ code. (Your code *is* well-written, right?)

    case "$GXX_VERSION" in
      egcs1.1)
          ;;

      # All other gcc versions
      *)
          CXXFLAGS_OPT="$CXXFLAGS_OPT -funroll-loops -fstrict-aliasing"
          CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -funroll-loops -fstrict-aliasing"

          CFLAGS_OPT="$CFLAGS_OPT -funroll-loops -fstrict-aliasing"
          CFLAGS_DEVEL="$CFLAGS_DEVEL -funroll-loops -fstrict-aliasing"
          ;;
    esac


    case "$GXX_VERSION" in
      # - after gcc2.95, some flags were deemed obsolete for C++
      #   (and are only supported for C any more), so only define them for
      #   previous compilers
      egcs1.1 | gcc2.95)
         CXXFLAGS_OPT="$CXXFLAGS_OPT -fnonnull-objects"
         CXXFLAGS_DBG="$CXXFLAGS_DBG -Wmissing-declarations -Wbad-function-cast -Wtraditional -Wnested-externs"
         ;;

      # - define additional debug flags for newer versions of gcc which support them.
      #
      # Note:  do not use -Wold-style-cast...  creates a lot of unavoidable warnings
      #        when dealing with C APIs that take void* pointers.
      gcc4.4 | gcc4.5 | gcc4.6)
	 CXXFLAGS_OPT="$CXXFLAGS_OPT -std=c++0x -Wdisabled-optimization"
         CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -std=c++0x -Woverloaded-virtual -Wdisabled-optimization"
	 CXXFLAGS_DBG="$CXXFLAGS_DBG -std=c++0x -Woverloaded-virtual"

         if (test "x$enableglibcxxdebugging" = "xyes"); then
           CPPFLAGS_DBG="$CPPFLAGS_DBG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"
         fi
	 ;;

      gcc3.* | gcc4.*)
	 CXXFLAGS_OPT="$CXXFLAGS_OPT -Wdisabled-optimization"
         CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -Woverloaded-virtual -Wdisabled-optimization"
	 CXXFLAGS_DBG="$CXXFLAGS_DBG -Woverloaded-virtual"

         if (test "x$enableglibcxxdebugging" = "xyes"); then
	   CPPFLAGS_DBG="$CPPFLAGS_DBG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"
	 fi
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

      MIPSpro)
          CXXFLAGS_OPT="-LANG:std -LANG:libc_in_namespace_std -no_auto_include -ansi -O2 -w"
          CXXFLAGS_DBG="-LANG:std -LANG:libc_in_namespace_std -no_auto_include -ansi -g -woff 1460"
	  CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
          NODEPRECATEDFLAG=""
          CFLAGS_OPT="-O2 -w"
          CFLAGS_DBG="-g"
          CFLAGS_DEVEL="$CFLAGS_DBG"

          # For some reason, CC forgets to add the math lib to the
          # linker line, so we do that ourselves
          LDFLAGS="$LDFLAGS -lm"

	  # Augment CXXFLAGS to include -LANG:std if not there.  This is
          # needed to compile the remaining configure tests
          if test "x`echo $CXXFLAGS | grep 'LANG:std'`" = "x" ; then
	    CXXFLAGS="$CXXFLAGS -LANG:std"
          fi
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

          # Intel ICC >= v11.x
 	  intel_icc_v11.x | intel_icc_v12.x | intel_icc_v13.x | intel_icc_v14.x)
              # Disable some warning messages:
              # #175: 'subscript out of range'
              #       FIN-S application code causes many false
              #       positives with this
              # #266: 'function declared implicitly'
              #       Metis function "GKfree" caused this error
              #       in almost every file.
              # #1476: 'field uses tail padding of a base class'
              # #1505: 'size of class is affected by tail padding'
              #        simply warns of a possible incompatibility with
              #        the g++ ABI for this case
              # #1572: 'floating-point equality and inequality comparisons are unreliable'
              #        Well, duh, when the tested value is computed...  OK when it
              #        was from an assignment.
              PROFILING_FLAGS="-p"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -w1 -g -wd175 -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -O3 -unroll -w0 -ftz"
              CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -w1 -g -wd175 -wd1476 -wd1505 -wd1572"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -g -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O3 -unroll -w0 -ftz"
              CFLAGS_DEVEL="$CFLAGS_DBG"
              ;;
          intel_icc_v10.1)
              # Disable some warning messages:
              # #175: 'subscript out of range'
              #       FIN-S application code causes many false
              #       positives with this
              # #266: 'function declared implicitly'
              #       Metis function "GKfree" caused this error
              #       in almost every file.
              # #1476: 'field uses tail padding of a base class'
              # #1505: 'size of class is affected by tail padding'
              #        simply warns of a possible incompatibility with
              #        the g++ ABI for this case
              # #1572: 'floating-point equality and inequality comparisons are unreliable'
              #        Well, duh, when the tested value is computed...  OK when it
              #        was from an assignment.
              PROFILING_FLAGS="-p"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -w1 -g -wd175 -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -O3 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -w1 -g -wd175 -wd1476 -wd1505 -wd1572"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -g -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O3 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CFLAGS_DEVEL="$CFLAGS_DBG"
              ;;

          # Intel ICC >= 10.0
          intel_icc_v10.0)
              # Disable some warning messages:
              # #266: 'function declared implicitly'
              #       Metis function "GKfree" caused this error
              #       in almost every file.
              # #1572: 'floating-point equality and inequality comparisons are unreliable'
              #        Well, duh, when the tested value is computed...  OK when it
              #        was from an assignment.
              # Note: In Version 8.1 (and possibly newer?) the -inline_debug_info
              #       option causes a segmentation fault in libmesh.C, probably others...

              # CPU-specific flags: -axK is for ia32, -xW is for x86_64
              INTEL_AX_FLAG="-tpp6 -axK"
              if test $target_cpu = "x86_64" ; then
                INTEL_AX_FLAG="-xW"
              fi

              PROFILING_FLAGS="-p"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -Kc++eh -Krtti -O1 -w1 -g -wd504 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -Kc++eh -Krtti -O2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -g -inline_debug_info -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CFLAGS_DEVEL="$CFLAGS_DBG"
              ;;

          # Intel ICC >= 8.1
          intel_icc_v8.1 | intel_icc_v9.0 | intel_icc_v9.1 | intel_icc_v10.0)
              # Disable some warning messages:
              # #266: 'function declared implicitly'
              #       Metis function "GKfree" caused this error
              #       in almost every file.
              # #1572: 'floating-point equality and inequality comparisons are unreliable'
              #        Well, duh, when the tested value is computed...  OK when it
              #        was from an assignment.
              # Note: In Version 8.1 (and possibly newer?) the -inline_debug_info
              #       option causes a segmentation fault in libmesh.C, probably others...

              # CPU-specific flags: -axK is for ia32, -xW is for x86_64
              INTEL_AX_FLAG="-tpp6 -axK"
              if test $target_cpu = "x86_64" ; then
                INTEL_AX_FLAG="-xW"
              fi

              CXXFLAGS_DBG="$CXXFLAGS_DBG -Kc++eh -Krtti -O1 -w1 -g -wd504 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -Kc++eh -Krtti -O2 -Ob2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -g -inline_debug_info -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O2 -Ob2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CFLAGS_DEVEL="$CFLAGS_DBG"
              ;;

          # Intel ICC < v8.1
          intel_icc*)
              # Disable some warning messages:
              # #266: 'function declared implicitly'
              #       Metis function "GKfree" caused this error
              #       in almost every file.
              CXXFLAGS_OPT="-Kc++eh -Krtti -O2 -Ob2 -tpp6 -axiMK -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CXXFLAGS_DEVEL="-Kc++eh -Krtti -O1 -w1 -inline_debug_info -g -wd504"
              CXXFLAGS_DBG="-Kc++eh -Krtti -O0 -w1 -inline_debug_info -g -wd504"
              CFLAGS_DBG="-w1 -inline_debug_info -wd266"
              CFLAGS_OPT="-O2 -Ob2 -tpp6 -axiMK -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CFLAGS_DEVEL="$CFLAGS_DBG"
              ;;

          # Intel Itanium ICC >= v10.1
          intel_itanium_icc_v10.1)
              # Disable some warning messages:
              # #266: 'function declared implicitly'
              #       Metis function "GKfree" caused this error
              #       in almost every file.
              # #1476: 'field uses tail padding of a base class'
              # #1505: 'size of class is affected by tail padding'
              #        simply warns of a possible incompatibility with
              #        the g++ ABI for this case
              # #1572: 'floating-point equality and inequality comparisons are unreliable'
              #        Well, duh, when the tested value is computed...  OK when it
              #        was from an assignment.
              CXXFLAGS_DBG="$CXXFLAGS_DBG -w1 -inline_debug_info -g -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -inline_debug_info -g -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CFLAGS_DEVEL="$CFLAGS_DBG"
              ;;

          intel_itanium_icc_v8.1 | intel_itanium_icc_v9.0 | intel_itanium_icc_v9.1 | intel_itanium_icc_v10.0)
              # Disable some warning messages:
              # #266: 'function declared implicitly'
              #       Metis function "GKfree" caused this error
              #       in almost every file.
              # #1476: 'field uses tail padding of a base class'
              # #1505: 'size of class is affected by tail padding'
              #        simply warns of a possible incompatibility with
              #        the g++ ABI for this case
              # #1572: 'floating-point equality and inequality comparisons are unreliable'
              #        Well, duh, when the tested value is computed...  OK when it
              #        was from an assignment.
              CXXFLAGS_DBG="$CXXFLAGS_DBG -Kc++eh -Krtti -w1 -inline_debug_info -g -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -Kc++eh -Krtti -O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -inline_debug_info -g -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CFLAGS_DEVEL="$CFLAGS_DBG"
              ;;

          # Intel Itanium ICC < v8.1
          intel_itanium_icc*)
              # Disable some warning messages:
              # #266: 'function declared implicitly'
              #       Metis function "GKfree" caused this error
              #       in almost every file.
              CXXFLAGS_DBG="-Kc++eh -Krtti -w1 -inline_debug_info -g"
              CXXFLAGS_OPT="-Kc++eh -Krtti -O2 -unroll -w0 -ftz"
              CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
              CFLAGS_DBG="-w1 -inline_debug_info -g -wd266"
              CFLAGS_OPT="-O2 -unroll -w0 -ftz"
              CFLAGS_DEVEL="$CFLAGS_DBG"
              ;;

          *)
	      AC_MSG_RESULT(Unknown Intel complier, "$GXX_VERSION")
              ;;
        esac
      ;;

      compaq_cxx)
          # Disable some warning messages:
          # #175: `subscript out of range' (detected when instantiating a
          #       template and looking at the indices of an array of
          #       template dependent size, this error is triggered in a
          #       branch that is not taken for the present space dimension)
          # #236 and
          # #237: `controlling expression is constant' (in while(true), or
          #       switch(dim))
          # #487: `Inline function ... cannot be explicitly instantiated'
          #       (also reported when we instantiate the entire class)
          # #1136:`conversion to integral type of smaller size could lose data'
          #       (occurs rather often in addition of int and x.size(),
          #       because the latter is size_t=long unsigned int on Alpha)
          # #1156:`meaningless qualifiers not compatible with "..." and "..."'
          #       (cause unknown, happens when taking the address of a
          #       template member function)
          # #111 and
          # #1182:`statement either is unreachable or causes unreachable code'
          #       (happens in switch(dim) clauses for other dimensions than
          #       the present one)
          #
          # Also disable the following error:
          # #265: `class "..." is inaccessible' (happens when we try to
          #       initialize a static member variable in terms of another
          #       static member variable of the same class if the latter is
          #       not public and therefore not accessible at global scope in
          #       general. I nevertheless think that this is valid.)
          #
          # Besides this, choose the most standard conforming mode of the
          # compiler, i.e. -model ansi and -std strict_ansi. Unfortunately,
          # we have to also add the flag -implicit_local (generating implicit
          # instantiations of template with the `weak' link flag) since
          # otherwise not all templates are instantiated (also some from the
          # standards library are missing).

          CXXFLAGS_DBG="-nousing_std -nocurrent_include -model ansi -std strict_ansi -w1 -msg_display_number -timplicit_local"
          CXXFLAGS_OPT="-nousing_std -nocurrent_include -model ansi -std strict_ansi -w2 -msg_display_number -timplicit_local -O2 -fast"
	  CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
          CFLAGS_DBG="-w1 -msg_display_number -timplicit_local"
          CFLAGS_OPT="-w2 -msg_display_number -timplicit_local -O2 -fast"
  	  CFLAGS_DEVEL="$CFLAGS_DBG"

          NODEPRECATEDFLAG=""

          for i in 175 236 237 487 1136 1156 111 1182 265 ; do
            CXXFLAGS_DBG="$CXXFLAGS_DBG -msg_disable $i"
            CXXFLAGS_OPT="$CXXFLAGS_OPT -msg_disable $i"
            CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -msg_disable $i"
          done

          # For some reason, cxx also forgets to add the math lib to the
          # linker line, so we do that ourselves
          LDFLAGS="$LDFLAGS -lm"
          ;;

      sun_studio | sun_forte)
          CXXFLAGS_DBG="-library=stlport4 -g"
          CXXFLAGS_OPT="-library=stlport4 -fast -xO4"
	  CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
          NODEPRECATEDFLAG=""
          CFLAGS_DBG="-g"
          CFLAGS_OPT="-xO4"
          CFLAGS_DEVEL="$CFLAGS_DBG"

          # librpcsvc for XDR
          LIBS="-lrpcsvc $LIBS"
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

      hpux_acc)
          # This is for aCC A.03.31
          # +DA2.0W requires that the code is only working on
          #  PA-RISC 2.0 systems, i.e. for 64bit
          # -ext allows various C++ extensions
          # +z Cause the compiler to generate position independent
          #  code (PIC) for use in building shared libraries.
          #  However, currently only static lib seems to work.
          # for aCC:
          #  -Aa turns on  newly supported ANSI C++ Standard features
          #  -AA turns on full new ANSI C++ (this includes -Aa)
          # for cc:
          #  -Aa compiles under true ANSI mode
          #  -Ae turns on ANSI C with some HP extensions
          CXXFLAGS_DBG="+DA2.0W -AA +z -ext -g"
          CXXFLAGS_OPT="+DA2.0W -AA +z -ext -O +Onolimit"
	  CXXFLAGS_DEVEL="$CXXFLAGS_DBG"
          NODEPRECATEDFLAG=""
	  CFLAGS_DBG="+DA2.0W -Aa +z -Ae -g"
          CFLAGS_OPT="+DA2.0W -Aa +z -Ae -O +Onolimit"
          CFLAGS_DEVEL="$CFLAGS_DBG"
          LDFLAGS="$LDFLAGS -I/usr/lib/pa20_64"
          LIBS="$LIBS -lrpcsvc"
          FLIBS="$FLIBS -lF90 -lcl -I/opt/fortran90/lib/pa20_64"
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
          CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -felide-constructors -Qunused-arguments"
          CXXFLAGS_DEVEL="$CXXFLAGS_DEVEL -O2 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -Wuninitialized -Qunused-arguments"
          CXXFLAGS_DBG="$CXXFLAGS_DBG -O0 -felide-constructors -g -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -Qunused-arguments"
          NODEPRECATEDFLAG="-Wno-deprecated"

          CFLAGS_OPT="-O2 -Qunused-arguments"
          CFLAGS_DEVEL="$CFLAGS_OPT -g -Wimplicit"
          CFLAGS_DBG="-g -Wimplicit -Qunused-arguments"
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
