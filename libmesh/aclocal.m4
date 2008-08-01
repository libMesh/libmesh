dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------
dnl



dnl -------------------------------------------------------------
dnl Determine the C++ compiler in use. Return the name and possibly
dnl version of this compiler in GXX_VERSION.
dnl
dnl Usage: DETERMINE_CXX_BRAND
dnl
dnl -------------------------------------------------------------
AC_DEFUN(DETERMINE_CXX_BRAND, dnl
[
  dnl First check for gcc version, avoids intel's icc from
  dnl pretending to be gcc
  REAL_GXX=`($CXX -v 2>&1) | grep "gcc version"`

  if (test "$GXX" = yes -a "x$REAL_GXX" != "x" ) ; then
    dnl find out the right version
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "gcc version"`
    case "$GXX_VERSION_STRING" in
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
  	dnl These compilers are too old to support a useful subset
  	dnl of modern C++, so we don't support them
  	AC_MSG_RESULT(<<< C++ compiler is $GXX_VERSION_STRING >>>)
  	AC_MSG_ERROR(<<< C++ compiler is not supported >>>)
  	;;
      *)
  	AC_MSG_RESULT(<<< C++ compiler is unknown but accepted gcc version >>>)
  	GXX_VERSION=gcc-other
  	;;
    esac
  else
    dnl Check other (non-gcc) compilers
  
    dnl Check for IBM xlC. For some reasons, depending on some environment
    dnl variables, moon position, and other reasons unknown to me, the
    dnl compiler displays different names in the first line of output, so
    dnl check various possibilities
    is_ibm_xlc="`($CXX 2>&1) | egrep 'xlc'`"
    if test "x$is_ibm_xlc" != "x"  ; then
      dnl IBM's C++ compiler.
      AC_MSG_RESULT(<<< C++ compiler is IBM xlC >>>)
      GXX_VERSION=ibm_xlc
    else
  
      dnl Check whether we are dealing with the MIPSpro C++ compiler
      is_mips_pro="`($CXX -version 2>&1) | grep MIPSpro`"
      if test "x$is_mips_pro" != "x" ; then
        AC_MSG_RESULT(<<< C++ compiler is MIPSpro C++ compiler >>>)
        GXX_VERSION=MIPSpro
      else

        dnl Intel's ICC C++ compiler for Itanium?
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

          dnl Intel's ICC C++ compiler?
          is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'Compiler'`"
          if test "x$is_intel_icc" != "x" ; then
            GXX_VERSION_STRING="`($CXX -V 2>&1) | grep 'Version '`"
            case "$GXX_VERSION_STRING" in
              *10.1*)
                AC_MSG_RESULT(<<< C++ compiler is Intel ICC 10.1 >>>)
  	        GXX_VERSION=intel_icc_v10.1
                ;;
              *10.0*)
                AC_MSG_RESULT(<<< C++ compiler is Intel ICC 10.0 >>>)
  	        GXX_VERSION=intel_icc_v10.0
                ;;
              *9.1*)
                AC_MSG_RESULT(<<< C++ compiler is Intel ICC 9.1 >>>)
  	        GXX_VERSION=intel_icc_v9.1
                ;;
              *9.0*)
                AC_MSG_RESULT(<<< C++ compiler is Intel ICC 9.0 >>>)
  	        GXX_VERSION=intel_icc_v9.0
                ;;
              *8.1*)
                AC_MSG_RESULT(<<< C++ compiler is Intel ICC 8.1 >>>)
  	        GXX_VERSION=intel_icc_v8.1
                ;;
              *8.0*)
                AC_MSG_RESULT(<<< C++ compiler is Intel ICC 8.0 >>>)
  	        GXX_VERSION=intel_icc_v8.0
                ;;
              *7.1*)
                AC_MSG_RESULT(<<< C++ compiler is Intel ICC 7.1 >>>)
  	        GXX_VERSION=intel_icc_v7.1
                ;;
              *7.0*)
                AC_MSG_RESULT(<<< C++ compiler is Intel ICC 7.0 >>>)
  	        GXX_VERSION=intel_icc_v7.0
                ;;
            esac
          else
  	
            dnl Or Compaq's cxx compiler?
            is_dec_cxx="`($CXX -V 2>&1) | grep 'Compaq C++'`"
            if test "x$is_dec_cxx" != "x" ; then
              AC_MSG_RESULT(<<< C++ compiler is Compaq cxx >>>)
              GXX_VERSION=compaq_cxx
            else
  	
  	      dnl Sun Studio?
              is_sun_cc="`($CXX -V 2>&1) | grep 'Sun C++'`"
              if test "x$is_sun_cc" != "x" ; then
                AC_MSG_RESULT(<<< C++ compiler is Sun Studio compiler >>>)
                GXX_VERSION=sun_studio
              else
  	
  	        dnl Sun Forte?
                is_sun_forte_cc="`($CXX -V 2>&1) | grep 'Forte'`"
                if test "x$is_sun_forte_cc" != "x" ; then
                  AC_MSG_RESULT(<<< C++ compiler is Sun Forte compiler >>>)
                  GXX_VERSION=sun_forte
                else
  	
  	          dnl Cray C++?
  	          is_cray_cc="`($CXX -V 2>&1) | grep 'Cray '`"
  	          if test "x$is_cray_cc" != "x" ; then
  	            AC_MSG_RESULT(<<< C++ compiler is Cray C++ >>>)
  	            GXX_VERSION=cray_cc
  	          else
  	
  	            dnl Portland Group C++?
  	            is_pgcc="`($CXX -V 2>&1) | grep 'Portland Group'`"
  	            if test "x$is_pgcc" != "x" ; then
  	              AC_MSG_RESULT(<<< C++ compiler is Portland Group C++ >>>)
  	              GXX_VERSION=portland_group
  	            else
	
                      dnl HP-UX 11.11 aCC?
                      is_hpux_acc="`($CXX -V 2>&1) | grep 'aCC: HP ANSI C++'`"
  	              if test "x$is_hpux_acc" != "x" ; then
  	                AC_MSG_RESULT(<<< C++ compiler is HP-UX C++ >>>)
    	                GXX_VERSION=hpux_acc
  	              else
	
	
                        dnl No recognized compiler found...
			dnl warn the user and continue
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
])





dnl -------------------------------------------------------------
dnl Set C++ compiler flags to their default values. They will be 
dnl modified according to other options in later steps of
dnl configuration
dnl
dnl CXXFLAGS_OPT    : flags for optimized mode
dnl CXXFLAGS_DVL    : flags for development mode
dnl CXXFLAGS_DBG    : flags for debug mode
dnl PROFILING_FLAGS : flags to enable code profiling
dnl
dnl Usage: SET_CXX_FLAGS
dnl
dnl (Note the CXXFLAGS and the CPPFLAGS used for further tests may
dnl  be augmented)
dnl -------------------------------------------------------------
AC_DEFUN(SET_CXX_FLAGS, dnl
[
  dnl Flag for creating shared objects; can be modified at a later stage
  if test "x$target_os" = "xdarwin9.3.0" ; then
    CXXFLAGS_OPT="-fno-common"
    CXXFLAGS_DVL="-fno-common"
    CXXFLAGS_DBG="-fno-common"
    CXXSHAREDFLAG="-dynamiclib -Wl,-undefined,dynamic_lookup"
    CSHAREDFLAG="-dynamiclib -Wl,-undefined,dynamic_lookup"
  else
    CXXSHAREDFLAG="-shared"
  fi

  dnl Flag to add directories to the dynamic library search path; can
  dnl be changed at a later stage
  RPATHFLAG="-Wl,-rpath,"

  dnl Flag for profiling mode; can me modified at a later stage
  PROFILING_FLAGS="-pg"

  dnl The -g flag is all OProfile needs to produce annotations
  OPROFILE_FLAGS="-g"

  dnl First the flags for gcc compilers
  if (test "$GXX" = yes -a "x$REAL_GXX" != "x" ) ; then
    CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -felide-constructors"
    CXXFLAGS_DVL="$CXXFLAGS_DVL -O2 -felide-constructors -g -ansi -pedantic -W -Wall -Wno-long-long -Wunused -Wpointer-arith -Wimplicit -Wformat -Wparentheses -Wuninitialized"
    CXXFLAGS_DBG="$CXXFLAGS_DBG -O0 -felide-constructors -g -ansi -pedantic -W -Wall -Wno-long-long -Wunused -Wpointer-arith -Wimplicit -Wformat -Wparentheses"

    CFLAGS_OPT="-O2"
    CFLAGS_DVL="$CFLAGS_OPT -g"
    CFLAGS_DBG="-g"

    dnl Position-independent code for shared libraries
    if test "$enableshared" = yes ; then
      CXXFLAGS_OPT="$CXXFLAGS_OPT -fPIC"
      CXXFLAGS_DVL="$CXXFLAGS_DVL -fPIC"
      CXXFLAGS_DBG="$CXXFLAGS_DBG -fPIC"

      CFLAGS_OPT="$CFLAGS_OPT -fPIC"
      CFLAGS_DVL="$CFLAGS_DVL -fPIC"
      CFLAGS_DBG="$CFLAGS_DBG -fPIC"
    fi

    dnl set some flags that are specific to some versions of the
    dnl compiler:
    dnl - egcs1.1 yielded incorrect code with some loop unrolling
    dnl - after egcs1.1, the optimization flag -fstrict-aliasing was
    dnl   introduced, which enables better optimizations for
    dnl   well-written C++ code. (Your code *is* well-written, right?) 
  
    case "$GXX_VERSION" in
      egcs1.1)
          ;;
  
      dnl All other gcc versions
      *)
          CXXFLAGS_OPT="$CXXFLAGS_OPT -funroll-loops -fstrict-aliasing"
          CXXFLAGS_DVL="$CXXFLAGS_DVL -funroll-loops -fstrict-aliasing"

          CFLAGS_OPT="$CFLAGS_OPT -funroll-loops -fstrict-aliasing"
          CFLAGS_DVL="$CFLAGS_DVL -funroll-loops -fstrict-aliasing"
          ;;
    esac
  
  
    case "$GXX_VERSION" in
      dnl - after gcc2.95, some flags were deemed obsolete for C++
      dnl   (and are only supported for C any more), so only define them for
      dnl   previous compilers
      egcs1.1 | gcc2.95)
         CXXFLAGS_OPT="$CXXFLAGS_OPT -fnonnull-objects"
         CXXFLAGS_DBG="$CXXFLAGS_DBG -Wmissing-declarations -Wbad-function-cast -Wtraditional -Wnested-externs"
         ;;

      dnl - define additional debug flags for newer versions of gcc which support them.
      dnl 
      dnl Note:  do not use -Wold-style-cast...  creates a lot of unavoidable warnings
      dnl        when dealing with C APIs that take void* pointers.
      gcc4.3)
	 CXXFLAGS_OPT="$CXXFLAGS_OPT -std=c++0x -Wdisabled-optimization"
         CXXFLAGS_DVL="$CXXFLAGS_DVL -std=c++0x -Woverloaded-virtual -Wdisabled-optimization"
         CXXFLAGS_DBG="$CXXFLAGS_DBG -std=c++0x -Woverloaded-virtual -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"
  	 ;;

      gcc3.* | gcc4.*)
	 CXXFLAGS_OPT="$CXXFLAGS_OPT -Wdisabled-optimization"
         CXXFLAGS_DVL="$CXXFLAGS_DVL -Woverloaded-virtual -Wdisabled-optimization"
         CXXFLAGS_DBG="$CXXFLAGS_DBG -Woverloaded-virtual -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"
  	 ;;
      *)
         ;;
    esac


    dnl Set OS-specific flags for linkers & other stuff
    case "$target" in

      dnl For Solaris we need to pass a different flag to the linker for specifying the
      dnl dynamic library search path and add -lrpcsvc to use XDR
      *solaris*)
          RPATHFLAG="-Wl,-R,"
          LIBS="-lrpcsvc $LIBS"
          ;;
  
      *)
          ;;
    esac
  
  
  else
    dnl Non-gcc compilers
  
    case "$GXX_VERSION" in
      ibm_xlc)
          CXXFLAGS_OPT="-O3 -qmaxmem=-1 -w -qansialias -Q=10 -qrtti=all -qstaticinline"
          CXXFLAGS_DBG="-qmaxmem=-1 -qansialias -qrtti=all -g -qstaticinline"
	  CXXFLAGS_DVL="$CXXFLAGS_DBG"
          CFLAGS_OPT="-O3 -qmaxmem=-1 -w -qansialias -Q=10"
          CFLAGS_DBG="-qansialias -g"
          CFLAGS_DVL="$CFLAGS_DBG"
	  CXXSHAREDFLAG="-G -qmkshrobj -bnoerrmsg"
	  CSHAREDFLAG="-G -qmkshrobj"
	  RPATHFLAG="-Qoption,link,-rpath,"
          ;;
  
      MIPSpro)
          CXXFLAGS_OPT="-LANG:std -LANG:libc_in_namespace_std -no_auto_include -ansi -O2 -w"
          CXXFLAGS_DBG="-LANG:std -LANG:libc_in_namespace_std -no_auto_include -ansi -g -woff 1460"
	  CXXFLAGS_DVL="$CXXFLAGS_DBG"
          CFLAGS_OPT="-O2 -w"
          CFLAGS_DBG="-g"
          CFLAGS_DVL="$CFLAGS_DBG"

          dnl For some reason, CC forgets to add the math lib to the
          dnl linker line, so we do that ourselves
          LDFLAGS="$LDFLAGS -lm"

          dnl Position-independent code for shared libraries
          if test "$enableshared" = yes ; then
            CXXFLAGS_OPT="$CXXFLAGS_OPT -KPIC"
            CXXFLAGS_DBG="$CXXFLAGS_DBG -KPIC"
            CXXFLAGS_DVL="$CXXFLAGS_DVL -KPIC"

            CFLAGS_OPT="$CFLAGS_OPT -KPIC"
            CFLAGS_DBG="$CFLAGS_DBG -KPIC"
            CFLAGS_DVL="$CFLAGS_DVL -KPIC"

            LDFLAGS="$LDFLAGS -KPIC"
          fi

	  dnl Augment CXXFLAGS to include -LANG:std if not there.  This is
          dnl needed to compile the remaining configure tests
          if test "x`echo $CXXFLAGS | grep 'LANG:std'`" = "x" ; then
	    CXXFLAGS="$CXXFLAGS -LANG:std"
          fi
          ;;

      dnl All Intel ICC/ECC flavors
      intel_*)

        dnl Intel compilers use -qp for profiling
        PROFILING_FLAGS="-qp"

	dnl The -g flag is all OProfile needs to produce annotations
        OPROFILE_FLAGS="-g"

        dnl Specific flags for specific versions
        case "$GXX_VERSION" in

          dnl Intel ICC >= v10.1
          intel_icc_v10.1)
              dnl Disable some warning messages:
              dnl #266: 'function declared implicitly'
              dnl       Metis function "GKfree" caused this error
              dnl       in almost every file.
              dnl #1476: 'field uses tail padding of a base class'
              dnl #1505: 'size of class is affected by tail padding'
              dnl        simply warns of a possible incompatibility with
              dnl        the g++ ABI for this case
              dnl #1572: 'floating-point equality and inequality comparisons are unreliable'
              dnl        Well, duh, when the tested value is computed...  OK when it
              dnl        was from an assignment.
              CXXFLAGS_DBG="-w1 -inline_debug_info -g -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="-O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="-w1 -inline_debug_info -wd266 -wd1572"
              CFLAGS_OPT="-O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CFLAGS_DVL="$CFLAGS_DBG"
              ;;
                    dnl Intel ICC >= 10.0	          
          intel_icc_v10.0)		
              dnl Disable some warning messages:
              dnl #266: 'function declared implicitly'
              dnl       Metis function "GKfree" caused this error
              dnl       in almost every file.
              dnl #1572: 'floating-point equality and inequality comparisons are unreliable'
              dnl        Well, duh, when the tested value is computed...  OK when it
              dnl        was from an assignment.
              dnl Note: In Version 8.1 (and possibly newer?) the -inline_debug_info
              dnl       option causes a segmentation fault in libmesh.C, probably others...

              dnl CPU-specific flags: -axK is for ia32, -xW is for x86_64
              INTEL_AX_FLAG="-tpp6 -axK"
              if test $target_cpu = "x86_64" ; then
                INTEL_AX_FLAG="-xW"
              fi

              CXXFLAGS_DBG="-Kc++eh -Krtti -O1 -w1 -g -wd504 -wd1572"
              CXXFLAGS_OPT="-Kc++eh -Krtti -O2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="-w1 -inline_debug_info -wd266 -wd1572"
              CFLAGS_OPT="-O2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CFLAGS_DVL="$CFLAGS_DBG"
              ;;
          
          dnl Intel ICC >= 8.1	
          intel_icc_v8.1 | intel_icc_v9.0 | intel_icc_v9.1 | intel_icc_v10.0)	
              dnl Disable some warning messages:
              dnl #266: 'function declared implicitly'
              dnl       Metis function "GKfree" caused this error
              dnl       in almost every file.
              dnl #1572: 'floating-point equality and inequality comparisons are unreliable'
              dnl        Well, duh, when the tested value is computed...  OK when it
              dnl        was from an assignment.
              dnl Note: In Version 8.1 (and possibly newer?) the -inline_debug_info
              dnl       option causes a segmentation fault in libmesh.C, probably others...

              dnl CPU-specific flags: -axK is for ia32, -xW is for x86_64
              INTEL_AX_FLAG="-tpp6 -axK"
              if test $target_cpu = "x86_64" ; then
                INTEL_AX_FLAG="-xW"
              fi

              CXXFLAGS_DBG="-Kc++eh -Krtti -O1 -w1 -g -wd504 -wd1572"
              CXXFLAGS_OPT="-Kc++eh -Krtti -O2 -Ob2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="-w1 -inline_debug_info -wd266 -wd1572"
              CFLAGS_OPT="-O2 -Ob2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CFLAGS_DVL="$CFLAGS_DBG"
              ;;
          
          dnl Intel ICC < v8.1
          intel_icc*)
              dnl Disable some warning messages:
              dnl #266: 'function declared implicitly'
              dnl       Metis function "GKfree" caused this error
              dnl       in almost every file.
              CXXFLAGS_OPT="-Kc++eh -Krtti -O2 -Ob2 -tpp6 -axiMK -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CXXFLAGS_DVL="-Kc++eh -Krtti -O1 -w1 -inline_debug_info -g -wd504"
              CXXFLAGS_DBG="-Kc++eh -Krtti -O0 -w1 -inline_debug_info -g -wd504"
              CFLAGS_DBG="-w1 -inline_debug_info -wd266"
              CFLAGS_OPT="-O2 -Ob2 -tpp6 -axiMK -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CFLAGS_DVL="$CFLAGS_DBG"
              ;;
          
          dnl Intel Itanium ICC >= v10.1
          intel_itanium_icc_v10.1)
              dnl Disable some warning messages:
              dnl #266: 'function declared implicitly'
              dnl       Metis function "GKfree" caused this error
              dnl       in almost every file.
              dnl #1476: 'field uses tail padding of a base class'
          	  dnl #1505: 'size of class is affected by tail padding'
              dnl        simply warns of a possible incompatibility with
              dnl        the g++ ABI for this case
              dnl #1572: 'floating-point equality and inequality comparisons are unreliable'
              dnl        Well, duh, when the tested value is computed...  OK when it
              dnl        was from an assignment.
              CXXFLAGS_DBG="-w1 -inline_debug_info -g -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="-O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="-w1 -inline_debug_info -wd266 -wd1572"
              CFLAGS_OPT="-O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CFLAGS_DVL="$CFLAGS_DBG"
              ;;
          
          intel_itanium_icc_v8.1 | intel_itanium_icc_v9.0 | intel_itanium_icc_v9.1 | intel_itanium_icc_v10.0)
              dnl Disable some warning messages:
              dnl #266: 'function declared implicitly'
              dnl       Metis function "GKfree" caused this error
              dnl       in almost every file.
              dnl #1476: 'field uses tail padding of a base class'
          	  dnl #1505: 'size of class is affected by tail padding'
              dnl        simply warns of a possible incompatibility with
              dnl        the g++ ABI for this case
              dnl #1572: 'floating-point equality and inequality comparisons are unreliable'
              dnl        Well, duh, when the tested value is computed...  OK when it
              dnl        was from an assignment.
              CXXFLAGS_DBG="-Kc++eh -Krtti -w1 -inline_debug_info -g -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="-Kc++eh -Krtti -O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="-w1 -inline_debug_info -wd266 -wd1572"
              CFLAGS_OPT="-O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CFLAGS_DVL="$CFLAGS_DBG"
              ;;
          
          dnl Intel Itanium ICC < v8.1
          intel_itanium_icc*)
              dnl Disable some warning messages:
              dnl #266: 'function declared implicitly'
              dnl       Metis function "GKfree" caused this error
              dnl       in almost every file.
              CXXFLAGS_DBG="-Kc++eh -Krtti -w1 -inline_debug_info -g"
              CXXFLAGS_OPT="-Kc++eh -Krtti -O2 -unroll -w0 -ftz"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="-w1 -inline_debug_info -wd266"
              CFLAGS_OPT="-O2 -unroll -w0 -ftz"
              CFLAGS_DVL="$CFLAGS_DBG"
              ;;

          *)
	      AC_MSG_RESULT(Unknown Intel complier, "$GXX_VERSION")
              ;;
        esac

        dnl Position-independent code for shared libraries
        if test "$enableshared" = yes ; then
	        
          dnl Specific flags for specific versions
          case "$GXX_VERSION" in

            dnl Intel ICC >= 10.0	
            intel_*_v10.*)	
              CXXFLAGS_OPT="$CXXFLAGS_OPT -fPIC"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -fPIC"
              CXXFLAGS_DVL="$CXXFLAGS_DVL -fPIC"
        
              CFLAGS_OPT="$CFLAGS_OPT -fPIC"
              CFLAGS_DBG="$CFLAGS_DBG -fPIC"
              CFLAGS_DVL="$CFLAGS_DVL -fPIC"

              LDFLAGS="$LDFLAGS -fPIC"
	      ;;
	    *)
              CXXFLAGS_OPT="$CXXFLAGS_OPT -KPIC"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -KPIC"
              CXXFLAGS_DVL="$CXXFLAGS_DVL -KPIC"
        
              CFLAGS_OPT="$CFLAGS_OPT -KPIC"
              CFLAGS_DBG="$CFLAGS_DBG -KPIC"
              CFLAGS_DVL="$CFLAGS_DVL -KPIC"

              LDFLAGS="$LDFLAGS -KPIC"
	      ;;
          esac
        fi
      ;;
  
      compaq_cxx)
          dnl Disable some warning messages:
          dnl #175: `subscript out of range' (detected when instantiating a
          dnl       template and looking at the indices of an array of
          dnl       template dependent size, this error is triggered in a
          dnl       branch that is not taken for the present space dimension)
          dnl #236 and
          dnl #237: `controlling expression is constant' (in while(true), or
          dnl       switch(dim))
          dnl #487: `Inline function ... cannot be explicitly instantiated'
          dnl       (also reported when we instantiate the entire class)
          dnl #1136:`conversion to integral type of smaller size could lose data'
          dnl       (occurs rather often in addition of int and x.size(),
          dnl       because the latter is size_t=long unsigned int on Alpha)
          dnl #1156:`meaningless qualifiers not compatible with "..." and "..."'
          dnl       (cause unknown, happens when taking the address of a
          dnl       template member function)
          dnl #111 and
          dnl #1182:`statement either is unreachable or causes unreachable code'
          dnl       (happens in switch(dim) clauses for other dimensions than
          dnl       the present one)
          dnl
          dnl Also disable the following error:
          dnl #265: `class "..." is inaccessible' (happens when we try to
          dnl       initialize a static member variable in terms of another
          dnl       static member variable of the same class if the latter is
          dnl       not public and therefore not accessible at global scope in
          dnl       general. I nevertheless think that this is valid.)
          dnl
          dnl Besides this, choose the most standard conforming mode of the
          dnl compiler, i.e. -model ansi and -std strict_ansi. Unfortunately,
          dnl we have to also add the flag -implicit_local (generating implicit
          dnl instantiations of template with the `weak' link flag) since
          dnl otherwise not all templates are instantiated (also some from the
          dnl standards library are missing).
  
          CXXFLAGS_DBG="-nousing_std -nocurrent_include -model ansi -std strict_ansi -w1 -msg_display_number -timplicit_local"
          CXXFLAGS_OPT="-nousing_std -nocurrent_include -model ansi -std strict_ansi -w2 -msg_display_number -timplicit_local -O2 -fast"
	  CXXFLAGS_DVL="$CXXFLAGS_DBG"
          CFLAGS_DBG="-w1 -msg_display_number -timplicit_local"
          CFLAGS_OPT="-w2 -msg_display_number -timplicit_local -O2 -fast"
  	  CFLAGS_DVL="$CFLAGS_DBG"

          for i in 175 236 237 487 1136 1156 111 1182 265 ; do
            CXXFLAGS_DBG="$CXXFLAGS_DBG -msg_disable $i"
            CXXFLAGS_OPT="$CXXFLAGS_OPT -msg_disable $i"
            CXXFLAGS_DVL="$CXXFLAGS_DVL -msg_disable $i"
          done
  
          dnl If we use -model ansi to compile the files, we also have to
          dnl specify it for linking
          dnl LDFLAGS="$LDFLAGS -model ansi"
  
          dnl For some reason, cxx also forgets to add the math lib to the
          dnl linker line, so we do that ourselves
          LDFLAGS="$LDFLAGS -lm"


          dnl Position-independent code for shared libraries
          if test "$enableshared" = yes ; then
            CXXFLAGS_OPT="$CXXFLAGS_OPT -shared"
            CXXFLAGS_DBG="$CXXFLAGS_DBG -shared"
            CXXFLAGS_DVL="$CXXFLAGS_DVL -shared"

            CFLAGS_OPT="$CFLAGS_OPT -shared"
            CFLAGS_DBG="$CFLAGS_DBG -shared"
            CFLAGS_DVL="$CFLAGS_DVL -shared"

            LDFLAGS="$LDFLAGS -shared"
          fi
          ;;
  
      sun_studio | sun_forte)
          CXXFLAGS_DBG="-library=stlport4 -g"
          CXXFLAGS_OPT="-library=stlport4 -fast -xO4"
	  CXXFLAGS_DVL="$CXXFLAGS_DBG"
          CFLAGS_DBG="-g"
          CFLAGS_OPT="-xO4"
          CFLAGS_DVL="$CFLAGS_DBG"

          CXXSHAREDFLAG="-G"
          CSHAREDFLAG="-G"

          dnl Linker flags & librpcsvc for XDR
          RPATHFLAG="-R"
          LIBS="-lrpcsvc $LIBS"

          dnl Position-independent code for shared libraries
          if test "$enableshared" = yes ; then
            CXXFLAGS_OPT="$CXXFLAGS_OPT -KPIC"
            CXXFLAGS_DBG="$CXXFLAGS_DBG -KPIC"
            CXXFLAGS_DVL="$CXXFLAGS_DVL -KPIC"

            CFLAGS_OPT="$CFLAGS_OPT -KPIC"
            CFLAGS_DBG="$CFLAGS_DBG -KPIC"
            CFLAGS_DVL="$CFLAGS_DVL -KPIC"

            LDFLAGS="$LDFLAGS -KPIC"
          fi
          ;;
  
      portland_group)
	  CXXFLAGS_DBG="-g --no_using_std"
          CXXFLAGS_OPT="-O2 --no_using_std -fast -Minform=severe"
	  CXXFLAGS_DVL="$CXXFLAGS_DBG"
	  CFLAGS_DBG="-g"
          CFLAGS_OPT="-O2"
          CFLAGS_DVL="$CFLAGS_DBG"

          dnl Disable exception handling if we don't use it
          if test "$enableexceptions" = no ; then
	    CXXFLAGS_DBG="$CXXFLAGS_DBG --no_exceptions"
            CXXFLAGS_OPT="$CXXFLAGS_OPT --no_exceptions"
          fi

          dnl Position-independent code for shared libraries
          if test "$enableshared" = yes ; then
            CXXFLAGS_OPT="$CXXFLAGS_OPT -fpic"
            CXXFLAGS_DBG="$CXXFLAGS_DBG -fpic"
            CXXFLAGS_DVL="$CXXFLAGS_DVL -fpic"
          
            CFLAGS_OPT="$CFLAGS_OPT -fpic"
            CFLAGS_DBG="$CFLAGS_DBG -fpic"
            CFLAGS_DVL="$CFLAGS_DVL -fpic"
          
            LDFLAGS="$LDFLAGS -fpic"
          fi

	  if test $target_cpu = "x86_64" ; then
	    CXXFLAGS_DBG="$CXXFLAGS_DBG -tp amd64"
	    CXXFLAGS_OPT="$CXXFLAGS_OPT -tp amd64"
	    CXXFLAGS_DVL="$CXXFLAGS_DVL -tp amd64"
	    CFLAGS_DBG="$CFLAGS_DBG -tp amd64"
	    CFLAGS_OPT="$CFLAGS_OPT -tp amd64"
	    CFLAGS_DVL="$CFLAGS_DVL -tp amd64"
          fi
          ;;

      hpux_acc)
          dnl This is for aCC A.03.31
          dnl +DA2.0W requires that the code is only working on
          dnl  PA-RISC 2.0 systems, i.e. for 64bit
          dnl -ext allows various C++ extensions
          dnl +z Cause the compiler to generate position independent
          dnl  code (PIC) for use in building shared libraries.
          dnl  However, currently only static lib seems to work.
          dnl for aCC:
          dnl  -Aa turns on  newly supported ANSI C++ Standard features
          dnl  -AA turns on full new ANSI C++ (this includes -Aa) 
          dnl for cc:
          dnl  -Aa compiles under true ANSI mode
          dnl  -Ae turns on ANSI C with some HP extensions
          CXXFLAGS_DBG="+DA2.0W -AA +z -ext -g"
          CXXFLAGS_OPT="+DA2.0W -AA +z -ext -O +Onolimit"
	  CXXFLAGS_DVL="$CXXFLAGS_DBG"
	  CFLAGS_DBG="+DA2.0W -Aa +z -Ae -g"
          CFLAGS_OPT="+DA2.0W -Aa +z -Ae -O +Onolimit"
          CFLAGS_DVL="$CFLAGS_DBG"
          LDFLAGS="$LDFLAGS -I/usr/lib/pa20_64"
          LIBS="$LIBS -lrpcsvc"
          FLIBS="$FLIBS -lF90 -lcl -I/opt/fortran90/lib/pa20_64"
          ;;

      cray_cc)
	  CXXFLAGS_DBG="-h conform,one_instantiation_per_object,instantiate=used,noimplicitinclude -G n"
 	  CXXFLAGS_OPT="-h conform,one_instantiation_per_object,instantiate=used,noimplicitinclude -G n"
 	  CXXFLAGS_DVL="-h conform,one_instantiation_per_object,instantiate=used,noimplicitinclude -G n"
	  CFLAGS_DBG="-G n"
	  CFLAGS_OPT="-G n"
	  CFLAGS_DVL="-G n"

	  CXXSHAREDFLAG=""
	  CSHAREDFLAG=""
	  RPATHFLAG=""
	  ;;

      *)
          AC_MSG_RESULT(No specific options for this C++ compiler known)
	  CXXFLAGS_DBG="$CXXFLAGS"
	  CXXFLAGS_OPT="$CXXFLAGS"
	  CXXFLAGS_DVL="$CXXFLAGS"

	  CFLAGS_DBG="$CFLAGS"
	  CFLAGS_OPT="$CFLAGS"
	  CFLAGS_DVL="$CFLAGS"
	  ;;
    esac
  fi
])



dnl -------------------------------------------------------------
dnl Petsc
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_PETSC, 
[
  AC_CHECK_FILE($PETSC_DIR/include/petsc.h,
                PETSC_H_PATH=$PETSC_DIR/include/petsc.h)

  dnl Test to see if PETSC_ARCH set by user.  If not set, then
  dnl disable petsc.
  if test "x$PETSC_ARCH" = x ; then
    enablepetsc=no
    AC_MSG_RESULT(<<< PETSc disabled.  Please set your "\$PETSC_ARCH" environment variable correctly. >>>)
    dnl PETSc config failed.  Try MPI.
    ACX_MPI

  else
    if (test -r $PETSC_DIR/include/petsc.h) ; then
      AC_ARG_WITH([f77],
      		  AC_HELP_STRING([--with-f77=F77],
                                 [Fortran compiler to use]),
      	          [F77="$withval"],
      	          [])	
      AC_PROG_F77            dnl Petsc requires linking with FORTRAN libraries 
      AC_F77_LIBRARY_LDFLAGS
      AC_SUBST(PETSC_ARCH)
      AC_SUBST(PETSC_DIR)
      AC_DEFINE(HAVE_PETSC, 1,
  	      [Flag indicating whether or not Petsc is available])
      AC_DEFINE(HAVE_MPI, 1,
  	      [Flag indicating whether or not MPI is available])
      MPI_IMPL="petsc_snooped"
  
      dnl Some tricks to discover the version of petsc.
      dnl You have to have grep and sed for this to work.
      petscmajor=`grep "define PETSC_VERSION_MAJOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MAJOR[ ]*//g"`
      petscminor=`grep "define PETSC_VERSION_MINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MINOR[ ]*//g"`
      petscsubminor=`grep "define PETSC_VERSION_SUBMINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_SUBMINOR[ ]*//g"`
      petscversion=$petscmajor.$petscminor.$petscsubminor
      petscmajorminor=$petscmajor.$petscminor.x
      AC_MSG_RESULT(<<< Configuring library with PETSc version $petscversion support >>>)

dnl      PETSCLINKLIBS=`cd $PETSC_DIR ; make getlinklibs`
dnl      PETSCINCLUDEDIRS=`cd $PETSC_DIR ; make getincludedirs`
dnl
dnl      AC_SUBST(PETSCLINKLIBS)
dnl      AC_SUBST(PETSCINCLUDEDIRS)

      AC_SUBST(petscversion)
      AC_SUBST(petscmajorminor)
      AC_SUBST(MPI_IMPL)
  
      else
  
      dnl PETSc config failed.  Try MPI.
      enablepetsc=no
      ACX_MPI
  
    fi
  fi
  AC_SUBST(enablepetsc)
])
dnl -------------------------------------------------------------


dnl -------------------------------------------------------------
dnl SLEPc
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_SLEPC,
[
  AC_CHECK_FILE($SLEPC_DIR/include/slepc.h,
                SLEPC_H_PATH=$SLEPC_DIR/include/slepc.h)
                                                                                       
  if (test -r $SLEPC_DIR/include/slepc.h) ; then
    AC_SUBST(SLEPC_DIR)

    dnl Similar to petsc, we need the slepc version number for.
    dnl Note slepc will most likely only work with the corresponding version of petsc
    slepcmajor=`grep "define SLEPC_VERSION_MAJOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_MAJOR[ ]*//g"`
    slepcminor=`grep "define SLEPC_VERSION_MINOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_MINOR[ ]*//g"`
    slepcsubminor=`grep "define SLEPC_VERSION_SUBMINOR" $SLEPC_DIR/include/slepcversion.h | sed -e "s/#define SLEPC_VERSION_SUBMINOR[ ]*//g"`
    slepcversion=$slepcmajor.$slepcminor.$slepcsubminor
    AC_SUBST(slepcversion)

    if (test $slepcversion != $petscversion) ; then
      AC_MSG_RESULT(WARNING:)
      AC_MSG_RESULT(>>> Different version numbers for SLEPc and PETSc <<<)
      AC_MSG_RESULT(Will skip SLEPc support)
      enableslepc=no
    else	
      AC_DEFINE(HAVE_SLEPC, 1,
                [Flag indicating whether or not SLEPc is available])
      AC_MSG_RESULT(<<< Configuring library with SLEPc version $slepcversion support >>>)
    fi
  else
    enableslepc=no
  fi
                                                                                       
  AC_SUBST(enableslepc)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Trilinos
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_TRILINOS, 
[
  if test "x$TRILINOS_DIR" = "x"; then
    TRILINOS_DIR=no
  fi  	

  AC_ARG_WITH(trilinos,
              AC_HELP_STRING([--with-trilinos=PATH],[Specify the path to Trilinos installation]),
              withtrilinosdir=$withval,
              withtrilinosdir=$TRILINOS_DIR)

  if test "$withtrilinosdir" != no ; then
    AC_CHECK_FILE($withtrilinosdir/include/Makefile.export.aztecoo,
                      AZTECOO_MAKEFILE_EXPORT=$withtrilinosdir/include/Makefile.export.aztecoo,
		      enabletrilinos=no)

    if test "$enabletrilinos" != no ; then
       AC_DEFINE(HAVE_TRILINOS, 1,
                 [Flag indicating whether the library shall be compiled to use the Trilinos solver collection])
       AC_MSG_RESULT(<<< Configuring library with Trilinos support >>>)
    fi
  else
    enabletrilinos=no
  fi

  AC_SUBST(AZTECOO_MAKEFILE_EXPORT)
  AC_SUBST(enabletrilinos)

])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Threading Building Blocks
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_TBB,
[
  AC_ARG_WITH(tbb,
              AC_HELP_STRING([--with-tbb=PATH],[Specify the path where Threading Building Blocks is installed]),
              withtbb=$withval,
              withtbb=no)

  AC_ARG_WITH(tbb-lib,
              AC_HELP_STRING([--with-tbb-lib=PATH],[Specify the path to Threading Building Blocks libraries]),
              withtbblib=$withval,
              withtbblib=no)

  if test "$withtbb" != no ; then
    AC_CHECK_FILE($withtbb/include/tbb/task_scheduler_init.h,
                      TBB_INCLUDE_PATH=$withtbb/include)
    if test "$withtbblib" != no ; then
      TBB_LIBS=$withtbblib
    else	
      TBB_LIBS=""
    fi
  fi

  if (test -r $TBB_INCLUDE_PATH/tbb/task_scheduler_init.h) ; then
    TBB_LIBRARY=$TBB_LIBS
    TBB_INCLUDE=-I$TBB_INCLUDE_PATH
    AC_SUBST(TBB_LIBRARY)
    AC_SUBST(TBB_INCLUDE)
    AC_DEFINE(HAVE_TBB_API, 1,
              [Flag indicating whether the library shall be compiled to use the Threading Building Blocks])
    AC_MSG_RESULT(<<< Configuring library with Intel TBB threading support >>>)
  fi
])
dnl -------------------------------------------------------------
                                                                                       


dnl -------------------------------------------------------------
dnl Space Filling Curves
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_SFC, 
[
  dnl Initialize variables
  SFC_INCLUDE=""
  SFC_LIB=""

  dnl Sanity check: make sure the user really has the contrib directory
  if (test $enablesfc = yes); then
    AC_CHECK_FILE(./contrib/sfcurves/sfcurves.h, [enablesfc=yes], [enablesfc=no])
  fi


  if (test $enablesfc = yes); then
     SFC_INCLUDE="-I$PWD/contrib/sfcurves"
     SFC_LIB="\$(EXTERNAL_LIBDIR)/libsfcurves\$(libext)"
     AC_DEFINE(HAVE_SFCURVES, 1, [Flag indicating whether or not Space filling curves are available])
     AC_MSG_RESULT(<<< Configuring library with SFC support >>>)
  fi

  AC_SUBST(SFC_INCLUDE)
  AC_SUBST(SFC_LIB)	
  AC_SUBST(enablesfc)

])
dnl -------------------------------------------------------------




dnl -------------------------------------------------------------
dnl Read/Write Compressed Streams with gzstream
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_GZ, 
[

dnl Initialize variables
GZSTREAM_INCLUDE=""
GZSTREAM_LIB=""

dnl Sanity check: make sure the user really has the contrib directory
if (test $enablegz = yes); then
  AC_CHECK_FILE(./contrib/gzstream/gzstream.h, [enablegz=yes], [enablegz=no])
fi


if (test $enablegz = yes); then
  dnl First check for the required system headers and libraries
  AC_CHECK_HEADERS(zlib.h, have_zlib_h=yes)
  AC_CHECK_LIB(z, gzopen, have_libz=yes)

  dnl If both tests succeded, continue the configuration process.
  if (test "$have_zlib_h" = yes -a "$have_libz" = yes) ; then
    GZSTREAM_INCLUDE="-I$PWD/contrib/gzstream"
    GZSTREAM_LIB="\$(EXTERNAL_LIBDIR)/libgzstream\$(libext) -lz"
    AC_DEFINE(HAVE_GZSTREAM, 1, [Flag indicating whether or not gzstreams are available])
    AC_MSG_RESULT(<<< Configuring library with gzstreams support >>>)

  dnl Otherwise do not enable gzstreams
  else
    enablegz=no;
  fi
fi

AC_SUBST(GZSTREAM_INCLUDE)
AC_SUBST(GZSTREAM_LIB)	
AC_SUBST(enablegz)

])
dnl -------------------------------------------------------------




dnl -------------------------------------------------------------
dnl LASPACK Iterative Solvers
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_LASPACK, 
[

dnl Initialize variables
LASPACK_INCLUDE=""
LASPACK_LIB=""

dnl Sanity check: make sure the user really has the contrib directory
if (test $enablelaspack = yes); then
  AC_CHECK_FILE(./contrib/laspack/lastypes.h, [enablelaspack=yes], [enablelaspack=no])
fi


if (test $enablelaspack = yes); then

  LASPACK_INCLUDE="-I$PWD/contrib/laspack"
  LASPACK_LIB="\$(EXTERNAL_LIBDIR)/liblaspack\$(libext)"
  AC_DEFINE(HAVE_LASPACK, 1, [Flag indicating whether or not LASPACK iterative solvers are available])
  laspack_version=`grep "define LASPACK_VERSION " $PWD/contrib/laspack/version.h | sed -e "s/[[^0-9.]]*//g"`
  AC_MSG_RESULT(<<< Configuring library with LASPACK version $laspack_version support >>>)

fi

AC_SUBST(LASPACK_INCLUDE)
AC_SUBST(LASPACK_LIB)	
AC_SUBST(enablelaspack)

])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Metis
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_METIS, 
[
  AC_CHECK_FILE(./contrib/metis/Lib/metis.h,
	        [ 
	          METIS_INCLUDE_PATH=$PWD/contrib/metis/Lib
                  METIS_INCLUDE=-I$METIS_INCLUDE_PATH
                  METIS_LIB="\$(EXTERNAL_LIBDIR)/libmetis\$(libext)"
		  AC_SUBST(METIS_INCLUDE)
                  AC_SUBST(METIS_LIB)
                  AC_DEFINE(HAVE_METIS, 1,
	                     [Flag indicating whether or not Metis is available])
                  AC_MSG_RESULT(<<< Configuring library with Metis support >>>)
	          enablemetis=yes
                ],
                [enablemetis=no])
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Parmetis
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_PARMETIS, 
[
  AC_REQUIRE([ACX_MPI])
	
  dnl We require a valid MPI installation for Parmetis
  if (test "x$MPI_IMPL" != x) ; then

    dnl need Metis for Parmetis
    AC_REQUIRE([CONFIGURE_METIS])

    if (test $enablemetis = yes) ; then
      AC_CHECK_FILE(./contrib/parmetis/Lib/parmetis.h,
      	        [
      		  
      	          PARMETIS_INCLUDE_PATH=$PWD/contrib/parmetis/Lib
                      PARMETIS_INCLUDE=-I$PARMETIS_INCLUDE_PATH
                      PARMETIS_LIB="\$(EXTERNAL_LIBDIR)/libparmetis\$(libext)"
      		  AC_SUBST(PARMETIS_INCLUDE)
                      AC_SUBST(PARMETIS_LIB)
                      AC_DEFINE(HAVE_PARMETIS, 1,
      	                     [Flag indicating whether or not ParMetis is available])
                      AC_MSG_RESULT(<<< Configuring library with ParMetis support >>>)
      	          enableparmetis=yes
                    ],
                    [enableparmetis=no])
    else
      enableparmetis=no
    fi  
  else
    enableparmetis=no
  fi
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Tecplot
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_TECPLOT,
[
  AC_ARG_WITH(tecplot,
              AC_HELP_STRING([--with-tecplot=PATH],[Specify the path where Tecplot is installed]),
              withtecplot=$withval,
              withtecplot=no)

  if test "$withtecplot" = no ; then
    AC_CHECK_FILE(./contrib/tecplot/lib/$host/tecio.a,
	  	  TECPLOT_LIBRARY_PATH=$PWD/contrib/tecplot/lib/$host)
    AC_CHECK_FILE(./contrib/tecplot/include/TECIO.h,
 	  	  TECPLOT_INCLUDE_PATH=$PWD/contrib/tecplot/include)
  else
    AC_CHECK_FILE($withtecplot/lib/tecio.a,
	  	  TECPLOT_LIBRARY_PATH=$withtecplot/lib)
    AC_CHECK_FILE($withtecplot/include/TECIO.h,
 	  	  TECPLOT_INCLUDE_PATH=$withtecplot/include)
  fi

  if (test -r $TECPLOT_LIBRARY_PATH/tecio.a -a -r $TECPLOT_INCLUDE_PATH/TECIO.h) ; then
    TECPLOT_LIBRARY=$TECPLOT_LIBRARY_PATH/tecio.a
    TECPLOT_INCLUDE=-I$TECPLOT_INCLUDE_PATH
    AC_SUBST(TECPLOT_LIBRARY)
    AC_SUBST(TECPLOT_INCLUDE)
    AC_DEFINE(HAVE_TECPLOT_API, 1,
              [Flag indicating whether the library shall be compiled to use the Tecplot interface])
    AC_MSG_RESULT(<<< Configuring library with Tecplot API support >>>)
  fi
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl TetGen tetrahedrization library
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_TETGEN, 
[
dnl if TetGen is enabled we need the header path and the lib

  if (test $enabletetgen = yes) ; then
     TETGEN_INCLUDE="-I$PWD/contrib/tetgen"
     TETGEN_LIBRARY="\$(EXTERNAL_LIBDIR)/libtetgen\$(libext)"
     AC_DEFINE(HAVE_TETGEN, 1, [Flag indicating whether the library will be compiled with TetGen support])
     AC_MSG_RESULT(<<< Configuring library with TetGen support >>>)
  else
     TETGEN_INCLUDE=""
     TETGEN_LIBRARY=""
     enabletetgen=no
   fi

  dnl TetGen
  AC_SUBST(TETGEN_INCLUDE)
  AC_SUBST(TETGEN_LIBRARY)	
  AC_SUBST(enabletetgen)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Triangle Delaunay triangulation library by J.R. Shewchuk
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_TRIANGLE, 
[
dnl Triangle is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enabletriangle = yes); then
     TRIANGLE_INCLUDE="-I$PWD/contrib/triangle"
     TRIANGLE_LIBRARY="\$(EXTERNAL_LIBDIR)/libtriangle\$(libext)"
     AC_DEFINE(HAVE_TRIANGLE, 1, [Flag indicating whether the library will be compiled with Triangle support])
     AC_MSG_RESULT(<<< Configuring library with Triangle support >>>)
  else
     TRIANGLE_INCLUDE=""
     TRIANGLE_LIBRARY=""
     enabletriangle=no
  fi

  AC_SUBST(TRIANGLE_INCLUDE)
  AC_SUBST(TRIANGLE_LIBRARY)	
  AC_SUBST(enabletriangle)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl GMV file I/O API for reading GMV files, by Frank Ortega
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_GMV, 
[
dnl The GMV API is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enablegmv = yes); then
     GMV_INCLUDE="-I$PWD/contrib/gmv"
     GMV_LIBRARY="\$(EXTERNAL_LIBDIR)/libgmv\$(libext)"
     AC_DEFINE(HAVE_GMV, 1, [Flag indicating whether the library will be compiled with GMV support])
     AC_MSG_RESULT(<<< Configuring library with GMV support >>>)
  else
     GMV_INCLUDE=""
     GMV_LIBRARY=""
     enablegmv=no
  fi

  AC_SUBST(GMV_INCLUDE)
  AC_SUBST(GMV_LIBRARY)	
  AC_SUBST(enablegmv)
])
dnl -------------------------------------------------------------



dnl ----------------------------------------------------------------
dnl VTK Mesh I/O API (by Wout Ruijter) requires VTK headers and lib
dnl ----------------------------------------------------------------
AC_DEFUN(CONFIGURE_VTK, 
[
  dnl Default path to VTK's include and lib files
  VTK_INC="/usr/include/vtk"
  VTK_LIB="/usr/lib"
  
  dnl User-specific include path
  AC_ARG_WITH(vtk-include,
              AC_HELP_STRING([--with-vtk-include=PATH],[Specify the path for VTK header files]),
              withvtkinc=$withval,
              withvtkinc=no)
	      
  dnl User-specific library path
  AC_ARG_WITH(vtk-lib,
              AC_HELP_STRING([--with-vtk-lib=PATH],[Specify the path for VTK libs]),
              withvtklib=$withval,
              withvtklib=no)

  if (test $withvtkinc != no); then
    VTK_INC="$withvtkinc"
  fi

  if (test $withvtklib != no); then
    VTK_LIB="$withvtklib"
  fi

  dnl Initialize Makefile/config.h substitution variables
  VTK_INCLUDE=""
  VTK_LIBRARY=""

  dnl Properly let the substitution variables
  if (test $enablevtk = yes); then
  
     dnl Check for existence of a header file in the specified location
     dnl AC_CHECK_FILE([$VTK_INC/vtkCommonInstantiator.h], [vtkincFound="OK"], [vtkincFound="FAIL"])
     vtkincFound=no;
     AC_CHECK_HEADERS($VTK_INC/vtkCommonInstantiator.h, vtkincFound=yes)

     if (test $vtkincFound = no); then
       AC_MSG_RESULT(VTK header files not found!)
       enablevtk=no;
     fi

     if (test $enablevtk = yes); then
       dnl Also Check for existence of required libraries
       AC_CHECK_FILE($VTK_LIB/libvtkIO.so, [enablevtk=yes], [enablevtk=no])
       AC_CHECK_FILE($VTK_LIB/libvtkCommon.so, [enablevtk=yes], [enablevtk=no])
     fi
     
     dnl If both the header file and the required libs were found, continue.
     if (test $enablevtk = yes); then
       dnl Since VTK headers use deprecated C++ header files and we don't want to see this
       dnl warning everytime, we can add the -Wno-derecated flag (GCC only?) to disable it.
       dnl Unfortunately, using -Wno-deprecated in the general CFLAGS generates *another*
       dnl warning when you are compiling C files.
       VTK_INCLUDE="-I$VTK_INC"
       VTK_LIBRARY="\$(libmesh_RPATHFLAG)$VTK_LIB -L$VTK_LIB -lvtkIO -lvtkCommon"
       AC_DEFINE(HAVE_VTK, 1, [Flag indicating whether the library will be compiled with VTK support])
       AC_MSG_RESULT(<<< Configuring library with VTK support >>>)
     fi
  fi

  dnl Substitute the substitution variables
  AC_SUBST(VTK_INCLUDE)
  AC_SUBST(VTK_LIBRARY)	
  AC_SUBST(enablevtk)
])
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl netCDF
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_NETCDF, 
[
dnl Netcdf is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enablenetcdf = yes); then
     NETCDF_INCLUDE="-I$PWD/contrib/netcdf/Lib"
     NETCDF_LIBRARY="\$(EXTERNAL_LIBDIR)/libnetcdf\$(libext)"
     AC_DEFINE(HAVE_NETCDF, 1, [Flag indicating whether the library will be compiled with Netcdf support])
     AC_MSG_RESULT(<<< Configuring library with Netcdf support >>>)
     have_netcdf=yes
  else
     NETCDF_INCLUDE=""
     NETCDF_LIBRARY=""
     enablenetcdf=no
     have_netcdf=no
  fi

  AC_SUBST(NETCDF_INCLUDE)
  AC_SUBST(NETCDF_LIBRARY)	
  AC_SUBST(enablenetcdf)
])
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl ExodusII
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_EXODUS, 
[
dnl Exodus is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enablenetcdf = yes -a $enableexodus = yes); then
     EXODUS_INCLUDE="-I$PWD/contrib/exodusii/Lib/include"
     EXODUS_LIBRARY="\$(EXTERNAL_LIBDIR)/libexodusii\$(libext)"
     AC_DEFINE(HAVE_EXODUS_API, 1, [Flag indicating whether the library will be compiled with Exodus support])
     AC_MSG_RESULT(<<< Configuring library with Exodus API support >>>)
  else
     EXODUS_INCLUDE=""
     EXODUS_LIBRARY=""
     enableexodus=no
  fi

  AC_SUBST(EXODUS_INCLUDE)
  AC_SUBST(EXODUS_LIBRARY)	
  AC_SUBST(enableexodus)
])
dnl -------------------------------------------------------------


dnl -------------------------------------------------------------
dnl Nemesis
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_NEMESIS, 
[
dnl Nemesis is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enablenetcdf = yes -a $enableexodus = yes -a $enablenemesis = yes); then
     NEMESIS_INCLUDE="-I$PWD/contrib/nemesis/Lib"
     NEMESIS_LIBRARY="\$(EXTERNAL_LIBDIR)/libnemesis\$(libext)"
     AC_DEFINE(HAVE_NEMESIS_API, 1, [Flag indicating whether the library will be compiled with Nemesis support])
     AC_MSG_RESULT(<<< Configuring library with Nemesis API support >>>)
  else
     NEMESIS_INCLUDE=""
     NEMESIS_LIBRARY=""
     enablenemesis=no
  fi

  AC_SUBST(NEMESIS_INCLUDE)
  AC_SUBST(NEMESIS_LIBRARY)	
  AC_SUBST(enablenemesis)
])
dnl -------------------------------------------------------------


dnl -------------------------------------------------------------
dnl libHilbert
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_LIBHILBERT,
[
dnl libHilbert is distributed with libmesh, so we don't have to guess
dnl where it might be installed...

  if (test $enablelibhilbert = yes); then
     LIBHILBERT_INCLUDE="-I$PWD/contrib/libHilbert/include"
     LIBHILBERT_LIBRARY="\$(EXTERNAL_LIBDIR)/libHilbert\$(libext)"
     AC_DEFINE(HAVE_LIBHILBERT, 1, [Flag indicating whether the library will be compiled with libHilbert support])
     AC_MSG_RESULT(<<< Configuring library with libHilbert support >>>)
  else
     LIBHILBERT_INCLUDE=""
     LIBHILBERT_LIBRARY=""
     enablelibhilbert=no
  fi

  AC_SUBST(LIBHILBERT_INCLUDE)
  AC_SUBST(LIBHILBERT_LIBRARY)	
  AC_SUBST(enablelibhilbert)
])
dnl -------------------------------------------------------------


# dnl -------------------------------------------------------------
# dnl netCDF
# dnl -------------------------------------------------------------
# AC_DEFUN(CONFIGURE_NETCDF,
# [
#   AC_CHECK_FILE(./contrib/netcdf/lib/$host/libnetcdf.a,
# 		NETCDF_LIB=$PWD/contrib/netcdf/lib/$host/libnetcdf.a)
#   AC_CHECK_FILE(./contrib/netcdf/include/netcdf.h,
# 		NETCDF_INCLUDE_PATH=$PWD/contrib/netcdf/include)

#   if (test -r $NETCDF_INCLUDE_PATH/netcdf.h -a "x$NETCDF_LIB" != x) ; then
#     NETCDF_INCLUDE=-I$NETCDF_INCLUDE_PATH
#     AC_SUBST(NETCDF_LIB)
#     AC_SUBST(NETCDF_INCLUDE)
#     AC_DEFINE(HAVE_NETCDF, 1,
#               [Flag indicating whether the library shall be compiled to support netcdf files])
#     AC_MSG_RESULT(<<< Configuring library with netCDF support >>>)
#     have_netcdf=yes
#   fi
# ])
# dnl -------------------------------------------------------------



# dnl -------------------------------------------------------------
# dnl ExodusII 
# dnl -------------------------------------------------------------
# AC_DEFUN(CONFIGURE_EXODUS,
# [
#   AC_CHECK_FILE(./contrib/exodus/lib/$host/libexoIIv2c.a,
# 		EXODUS_LIB=$PWD/contrib/exodus/lib/$host/libexoIIv2c.a)
#   AC_CHECK_FILE(./contrib/exodus/include/exodusII.h,
# 		EXODUS_INCLUDE_PATH=$PWD/contrib/exodus/include)

#   if (test -r $EXODUS_INCLUDE_PATH/exodusII.h -a "x$EXODUS_LIB" != x) ; then
#     EXODUS_INCLUDE=-I$EXODUS_INCLUDE_PATH
#     AC_SUBST(EXODUS_LIB)
#     AC_SUBST(EXODUS_INCLUDE)
#     AC_DEFINE(HAVE_EXODUS_API, 1,
# 	      [Flag indicating whether the library shall be compiled to use the Exodus interface])
#     AC_MSG_RESULT(<<< Configuring library with Exodus API support >>>)
#   fi
# ])
# dnl -------------------------------------------------------------






dnl -------------------------------------------------------------
dnl Contributed tests -- see
dnl http://www.gnu.org/software/ac-archive/htmldoc/index.html
dnl -------------------------------------------------------------




dnl -------------------------------------------------------------
dnl AC_CXX_HAVE_NAMESPACES
dnl -------------------------------------------------------------
AC_DEFUN([AC_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ac_cv_cxx_namespaces,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([namespace Outer { namespace Inner { int i = 0; }}],
                [using namespace Outer::Inner; return i;],
 ac_cv_cxx_namespaces=yes, ac_cv_cxx_namespaces=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])


dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl AC_CXX_HAVE_LOCALE
dnl -------------------------------------------------------------
AC_DEFUN([AC_CXX_HAVE_LOCALE],
[AC_CACHE_CHECK(whether the compiler has locale,
ac_cv_cxx_have_locale,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <locale>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],[locale loc; return 0;],
 ac_cv_cxx_have_locale=yes, ac_cv_cxx_have_locale=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_locale" = yes; then
  AC_DEFINE(HAVE_LOCALE,,[define if the compiler has locale])
fi
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl AC_CXX_HAVE_SSTREAM
dnl -------------------------------------------------------------
AC_DEFUN([AC_CXX_HAVE_SSTREAM],
[AC_CACHE_CHECK(whether the compiler has stringstream,
ac_cv_cxx_have_sstream,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <sstream>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],[stringstream message; message << "Hello"; return 0;],
 ac_cv_cxx_have_sstream=yes, ac_cv_cxx_have_sstream=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_sstream" = yes; then
  AC_DEFINE(HAVE_SSTREAM,,[define if the compiler has the sstream header])
  AC_DEFINE(HAVE_STRINGSTREAM,,[define if the compiler has stringstream functionality])
else
  dnl Some compilers implement strstream instead of sstream.  Check for it.
  AC_CXX_HAVE_STRSTREAM
fi
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl AC_CXX_HAVE_STRSTREAM
dnl -------------------------------------------------------------
AC_DEFUN([AC_CXX_HAVE_STRSTREAM],
[AC_CACHE_CHECK(whether the compiler has strstream,
ac_cv_cxx_have_strstream,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <strstream>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],[strstream message; message << "Hello"; return 0;],
 ac_cv_cxx_have_strstream=yes, ac_cv_cxx_have_strstream=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_strstream" = yes; then
  AC_DEFINE(HAVE_STRSTREAM,,[define if the compiler has the strstream header])
  AC_DEFINE(HAVE_STRINGSTREAM,,[define if the compiler has stringstream functionality])  
fi
])
dnl -------------------------------------------------------------



dnl ----------------------------------------------------------------------------
dnl @synopsis ACX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/).
dnl On success, it sets the BLAS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @version acsite.m4,v 1.3 2002/08/02 09:28:12 steve Exp
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_blas_ok=no
acx_blas_save_LIBS="$LIBS"

AC_ARG_WITH(blas,
            AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>]))
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
        -* | */*) LIBS = "$with_blas $LIBS";;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
	LIBS="$save_LIBS"
fi

# BLAS in Intel MKL libraries?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(mkl, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lmkl -lguide -lpthread"])
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl ACX_BLAS





dnl ----------------------------------------------------------------------------
dnl @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/).
dnl On success, it sets the LAPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_LAPACK.
dnl
dnl @version acsite.m4,v 1.3 2002/08/02 09:28:12 steve Exp
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>

AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no

AC_ARG_WITH(lapack,
            AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>]))
case $with_lapack in
        yes | "") ;;
        no) acx_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
        acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
        AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
        AC_MSG_RESULT($acx_lapack_ok)
        LIBS="$save_LIBS"
        if test acx_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($cheev, [acx_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $acx_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $cheev,
                    [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK





dnl ---------------------------------------------------------------------------
dnl check for the required MPI library
dnl ---------------------------------------------------------------------------
AC_DEFUN([ACX_MPI], [

# if MPIHOME is empty, set it to /usr
if (test "x$MPIHOME" = x) ; then
  MPIHOME="/usr"
fi

AC_ARG_WITH([mpi],
	    AC_HELP_STRING([--with-mpi=PATH],
                           [Prefix where MPI is installed (MPIHOME)]),
	    [MPI="$withval"],
	    [
              echo "note: MPI library path not given... trying prefix=$MPIHOME"
	      MPI=$MPIHOME
            ])

if test -z "$MPI"; then
	MPI="/usr"
fi



MPI_LIBS_PATH="$MPI/lib"
MPI_INCLUDES_PATH="$MPI/include"

# Check that the compiler uses the library we specified...

# look for LAM or other MPI implementation
if (test -e $MPI_LIBS_PATH/libmpi.a || test -e $MPI_LIBS_PATH/libmpi.so) ; then
	echo "note: using $MPI_LIBS_PATH/libmpi(.a/.so)"


	# Ensure the comiler finds the library...
	tmpLIBS=$LIBS
	AC_LANG_SAVE
	AC_LANG_CPLUSPLUS

	LIBS="-L$MPI_LIBS_PATH $LIBS"

	# look for lam_version_show in liblam.(a/so)
	# (this is needed in addition to libmpi.(a/so) for
        # LAM MPI
	AC_CHECK_LIB([lam],
                     [lam_show_version],
                     [
                       LIBS="-llam $LIBS"
                       MPI_LIBS="-llam $MPI_LIBS"
                     ],	
                     [])

	# Quadricss MPI requires the elan library to be included too
	if (nm $MPI_LIBS_PATH/libmpi.* | grep elan > /dev/null); then
	  echo "note: MPI found to use Quadrics switch, looking for elan library"
		 AC_CHECK_LIB([elan],
	                      [elan_init],
	                      [
                                LIBS="-lelan $LIBS"
                                MPI_LIBS="-lelan $MPI_LIBS"
                              ],
	                      [AC_MSG_ERROR( [Could not find elan library... exiting] )] )
	fi

	AC_CHECK_LIB([mpi],
                     [MPI_Init],                     
                     [
		       MPI_LIBS="-lmpi $MPI_LIBS"
	               MPI_LIBS_PATHS="-L$MPI_LIBS_PATH"
	               MPI_IMPL="mpi"
                       AC_MSG_RESULT([Found valid MPI installlaion...])
                     ],
                     [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no] )

	AC_LANG_RESTORE
	LIBS=$tmpLIBS
fi

if (test -e $MPI_LIBS_PATH/libmpich.a || test -e $MPI_LIBS_PATH/libmpich.so) ; then
	echo "note: using $MPI_LIBS_PATH/libmpich(.a/.so)"

	# Ensure the comiler finds the library...
	tmpLIBS=$LIBS
	AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	LIBS="-L$MPI_LIBS_PATH $LIBS"

	# Myricomm MPICH requires the gm library to be included too
	if (nm $MPI_LIBS_PATH/libmpich.* | grep gm_open > /dev/null); then
	  echo "note: MPICH found to use Myricomm's Myrinet, looking for gm library"

          if (test "x$GMHOME" = x) ; then
            GMHOME="/usr"
          fi 
          AC_ARG_WITH([gm],
	              AC_HELP_STRING([--with-gm=PATH],
                                     [Prefix where GM is installed (GMHOME)]),
		      [GM="$withval"],
		      [
                        echo "note: GM library path not given... trying prefix=$MPIHOME"
	                GM=$GMHOME
                      ])

          LIBS="-L$GM/lib $LIBS"
          MPI_LIBS="-L$GM/lib $MPI_LIBS"

          AC_CHECK_LIB([gm],
	               [gm_open],
	               [
                         LIBS="$LIBS -lgm"
                         MPI_LIBS="$MPI_LIBS -lgm"
                       ],
	               [AC_MSG_ERROR( [Could not find gm library... exiting] )] )
	fi

	# look for MPI_Init in libmpich.(a/so)
	AC_CHECK_LIB([mpich],
		     [MPI_Init],
		     [
		       MPI_LIBS="-lmpich $MPI_LIBS"
		       MPI_LIBS_PATHS="-L$MPI_LIBS_PATH"
	               MPI_IMPL="mpich"
                       AC_MSG_RESULT([Found valid MPICH installation...])
                     ],
		     [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no] )

	AC_LANG_RESTORE
	LIBS=$tmpLIBS
fi

if (test "x$MPI_IMPL" != x) ; then

	# Ensure the comiler finds the header file...
	if test -e $MPI_INCLUDES_PATH/mpi.h; then
		echo "note: using $MPI_INCLUDES_PATH/mpi.h"
		tmpCPPFLAGS=$CPPFLAGS
		AC_LANG_SAVE
		AC_LANG_CPLUSPLUS
		CPPFLAGS="-I$MPI_INCLUDES_PATH $CPPFLAGS"
		AC_CHECK_HEADER([mpi.h],
			        [AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])],
			        [AC_MSG_RESULT([Could not compile in the MPI headers...]); enablempi=no] )
		MPI_INCLUDES_PATHS="-I$MPI_INCLUDES_PATH"
		AC_LANG_RESTORE
		CPPFLAGS=$tmpCPPFLAGS
	else
		AC_MSG_RESULT([Could not find MPI header <mpi.h>...])
                enablempi=no
	fi
else
   
	# no MPI install found, see if the compiler supports it
      	AC_TRY_COMPILE([#include <mpi.h>],
	  	       [int np; MPI_Comm_size (MPI_COMM_WORLD, &np);],
                       [
	                 MPI_IMPL="built-in"
                         AC_MSG_RESULT( [$CXX Compiler Supports MPI] )
                         AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                       ],
                       [AC_MSG_RESULT([$CXX Compiler Does NOT Support MPI...]); enablempi=no] )                     	
fi 

# Save variables...
AC_SUBST(MPI)
AC_SUBST(MPI_IMPL)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_LIBS_PATH)
AC_SUBST(MPI_LIBS_PATHS)
AC_SUBST(MPI_INCLUDES_PATH)
AC_SUBST(MPI_INCLUDES_PATHS)
])dnl ACX_MPI



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::tr1::unordered_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_TR1_UNORDERED_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::tr1::unordered_map,
ac_cv_cxx_tr1_unordered_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <tr1/unordered_map>],
[
  std::tr1::unordered_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_tr1_unordered_map=yes, ac_cv_cxx_tr1_unordered_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_tr1_unordered_map" = yes; then
  AC_DEFINE(HAVE_TR1_UNORDERED_MAP,1,
            [define if the compiler supports std::tr1::unordered_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::tr1::unordered_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_TR1_UNORDERED_SET],
[AC_CACHE_CHECK(whether the compiler supports std::tr1::unordered_set,
ac_cv_cxx_tr1_unordered_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <tr1/unordered_set>],
[
  std::tr1::unordered_set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_tr1_unordered_set=yes, ac_cv_cxx_tr1_unordered_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_tr1_unordered_set" = yes; then
  AC_DEFINE(HAVE_TR1_UNORDERED_SET,1,
            [define if the compiler supports std::tr1::unordered_set])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_UNORDERED_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_map,
ac_cv_cxx_unordered_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <unordered_map>],
[
  std::unordered_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_unordered_map=yes, ac_cv_cxx_unordered_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_unordered_map" = yes; then
  AC_DEFINE(HAVE_UNORDERED_MAP,1,
            [define if the compiler supports std::unordered_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::unordered_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_UNORDERED_SET],
[AC_CACHE_CHECK(whether the compiler supports std::unordered_set,
ac_cv_cxx_unordered_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <unordered_set>],
[
  std::unordered_set<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_unordered_set=yes, ac_cv_cxx_unordered_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_unordered_set" = yes; then
  AC_DEFINE(HAVE_UNORDERED_SET,1,
            [define if the compiler supports std::unordered_set])
fi
])




dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl __gnu_cxx::hash_map
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_EXT_HASH_MAP],
[AC_CACHE_CHECK(whether the compiler supports __gnu_cxx::hash_map,
ac_cv_cxx_ext_hash_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <ext/hash_map>],
[
  __gnu_cxx::hash_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_ext_hash_map=yes, ac_cv_cxx_ext_hash_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_ext_hash_map" = yes; then
  AC_DEFINE(HAVE_EXT_HASH_MAP,1,
            [define if the compiler supports __gnu_cxx::hash_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl __gnu_cxx::hash_set
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_EXT_HASH_SET],
[AC_CACHE_CHECK(whether the compiler supports __gnu_cxx::hash_set,
ac_cv_cxx_ext_hash_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <ext/hash_set>],
[
  __gnu_cxx::hash_set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_ext_hash_set=yes, ac_cv_cxx_ext_hash_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_ext_hash_set" = yes; then
  AC_DEFINE(HAVE_EXT_HASH_SET,1,
            [define if the compiler supports __gnu_cxx::hash_set])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::hash_map (This is unlikely...)
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_HASH_MAP],
[AC_CACHE_CHECK(whether the compiler supports std::hash_map,
ac_cv_cxx_hash_map,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <hash_map>],
[
  std::hash_map<int, int> m;
  m.insert(std::make_pair(1, 2));
],
 ac_cv_cxx_hash_map=yes, ac_cv_cxx_hash_map=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_hash_map" = yes; then
  AC_DEFINE(HAVE_HASH_MAP,1,
            [define if the compiler supports std::hash_map])
fi
])



dnl ----------------------------------------------------------------------------
dnl Check to see if the compiler can compile a test program using
dnl std::hash_set (This is unlikely...)
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_HASH_SET],
[AC_CACHE_CHECK(whether the compiler supports std::hash_set,
ac_cv_cxx_hash_set,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <hash_set>],
[
  std::hash_set<int> m;
  m.insert(1);  m.insert(2);
],
 ac_cv_cxx_hash_set=yes, ac_cv_cxx_hash_set=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_hash_set" = yes; then
  AC_DEFINE(HAVE_HASH_SET,1,
            [define if the compiler supports std::hash_set])
fi
])


dnl ----------------------------------------------------------------------------
dnl check for gcc name demangling function
dnl Copyright © 2004 Neil Ferguson <nferguso@eso.org>
dnl
dnl Copying and distribution of this file, with or without
dnl modification, are permitted in any medium without royalty
dnl provided the copyright notice and this notice are preserved.
dnl ----------------------------------------------------------------------------
AC_DEFUN([AX_CXX_GCC_ABI_DEMANGLE],
[AC_CACHE_CHECK(whether the compiler supports GCC C++ ABI name demangling,
ac_cv_cxx_gcc_abi_demangle,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <typeinfo>
#include <cxxabi.h>
#include <string>

template<typename TYPE>
class A {};
],[A<int> instance;
int status = 0;
char* c_name = 0;

c_name = abi::__cxa_demangle(typeid(instance).name(), 0, 0, &status);

std::string name(c_name);
free(c_name);

return name == "A<int>";
],
 ac_cv_cxx_gcc_abi_demangle=yes, ac_cv_cxx_gcc_abi_demangle=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_gcc_abi_demangle" = yes; then
  AC_DEFINE(HAVE_GCC_ABI_DEMANGLE,1,
            [define if the compiler supports GCC C++ ABI name demangling])
fi
])


dnl ----------------------------------------------------------------------------
dnl check for gcc backtrace functions
dnl ----------------------------------------------------------------------------
AC_DEFUN([AX_CXX_GLIBC_BACKTRACE],
[AC_CACHE_CHECK(whether the c++ compiler supports glibc backtrace,
ac_cv_cxx_glibc_backtrace,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <execinfo.h>],
[void *addresses[10];
int size = backtrace(addresses, 10);
char** strings = backtrace_symbols(addresses, size);
],
 ac_cv_cxx_glibc_backtrace=yes, ac_cv_cxx_glibc_backtrace=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_glibc_backtrace" = yes; then
  AC_DEFINE(HAVE_GLIBC_BACKTRACE,1,
            [define if the compiler supports glibc backtrace])
fi
])


dnl ----------------------------------------------------------------------------
dnl check for the required PETSc library
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_PETSc], [
AC_REQUIRE([ACX_MPI])
AC_REQUIRE([ACX_LAPACK])
BLAS_LIBS="$BLAS_LIBS $FLIBS"
LAPACK_LIBS="$LAPACK_LIBS $BLAS_LIBS"
AC_PATH_XTRA
X_LIBS="$X_PRE_LIBS $X_LIBS -lX11 $X_EXTRA_LIBS"

# Set variables...
AC_ARG_WITH([PETSc],
	    AC_ARG_HELP([--with-PETSc=PATH],
                        [Prefix where PETSc is installed (PETSC_DIR)]),
	    [PETSc="$withval"],
	    [
              if test $PETSC_DIR; then
		PETSc="$PETSC_DIR"
		echo "note: assuming PETSc library is in $PETSc (/lib,/include) as specified by environment variable PETSC_DIR"
	      else
		PETSc="/usr/local"
		echo "note: assuming PETSc library is in /usr/local (/lib,/include)"
	      fi
            ])

AC_ARG_WITH([BOPT],
	    AC_ARG_HELP([--with-BOPT=VAL],[BOPT setting for PETSc (BOPT)]),
 	    [BOPT="$withval"],
	    [
              echo "note: assuming BOPT to O"
	      BOPT="O"
            ])

AC_ARG_WITH([PETSc_ARCH],
	    AC_ARG_HELP([--with-PETSc_ARCH=VAL],[PETSc hardware architecture (PETSC_ARCH)]),
	    [PETSc_ARCH="$withval"],
	    [
              if test $PETSC_ARCH; then
		PETSc_ARCH="$PETSC_ARCH"
		echo "note: assuming PETSc hardware architecture to be $PETSc_ARCH as specified by environment variable PETSC_ARCH"
	      else
		PETSc_ARCH=`uname -p`
		echo "note: assuming PETSc hardware architecture to be $PETSc_ARCH"
	      fi
            ])

PETSc_LIBS_PATH="$PETSc/lib/lib$BOPT/$PETSc_ARCH"
PETSc_INCLUDES_PATH="$PETSc/include"

# Check that the compiler uses the library we specified...
if test -e $PETSc_LIBS_PATH/libpetsc.a || test -e $PETSc_LIBS_PATH/libpetsc.so; then
	echo "note: using $PETSc_LIBS_PATH/libpetsc (.a/.so)"
else
	AC_MSG_ERROR( [Could not physically find PETSc library... exiting] )
fi 
if test -e $PETSc_INCLUDES_PATH/petsc.h; then
	echo "note: using $PETSc_INCLUDES_PATH/petsc.h"
else
	AC_MSG_ERROR( [Could not physically find PETSc header file... exiting] )
fi 

# Ensure the comiler finds the library...
tmpLIBS=$LIBS
tmpCPPFLAGS=$CPPFLAGS
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_CHECK_LIB(
	[dl],
	[dlopen],
	[DL_LIBS="-ldl"],
	[DL_LIBS=""; echo "libdl not found, assuming not needed for this architecture"] )
LIBS="-L$PETSc_LIBS_PATH $MPI_LIBS_PATHS $MPI_LIBS $LAPACK_LIBS $X_LIBS $LIBS -lm $DL_LIBS"
CPPFLAGS="$MPI_INCLUDES_PATHS -I$PETSc_INCLUDES_PATH -I$PETSc/bmake/$PETSc_ARCH $CPPFLAGS"
echo "cppflags=$CPPFLAGS"

AC_CHECK_LIB(
	[petsc],
	[PetscError],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc library... exiting] )] )
AC_CHECK_LIB(
	[petscvec],
	[ISCreateGeneral],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscvec library... exiting] )] )
AC_CHECK_LIB(
	[petscmat],
	[MAT_Copy],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscmat library... exiting] )] )
AC_CHECK_LIB(
	[petscdm],
	[DMInitializePackage],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscdm library... exiting] )] )
AC_CHECK_LIB(
	[petscsles],
	[SLESCreate],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscsles library... exiting] )] )
AC_CHECK_LIB(
	[petscsnes],
	[SNESCreate],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscsnes library... exiting] )] )
AC_CHECK_LIB(
	[petscts],
	[TSCreate],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscts library... exiting] )] )
AC_CHECK_LIB(
	[petscmesh],
	[MESH_CreateFullCSR],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscmesh library... exiting] )] )
AC_CHECK_LIB(
	[petscgrid],
	[GridCreate],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscgrid library... exiting] )] )
AC_CHECK_LIB(
	[petscgsolver],
	[GSolverInitializePackage],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petscgsolver library... exiting] )] )
AC_CHECK_LIB(
	[petscfortran],
	[meshcreate_],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc library... exiting] )] )
	AC_CHECK_LIB(
	[petsccontrib],
	[SDACreate1d],
	[],
	[AC_MSG_ERROR( [Could not link in the PETSc petsccontrib library... exiting] )] )
AC_CHECK_HEADER(
	[petsc.h],
	[AC_DEFINE( 
		[HAVE_PETSC],,
		[Define to 1 if you have the <petsc.h> header file.])],
	[AC_MSG_ERROR( [Could not compile in the PETSc headers... exiting] )] )
PETSc_LIBS="-lpetsc -lpetscvec -lpetscmat -lpetscdm -lpetscsles -lpetscsnes \
	-lpetscts -lpetscmesh -lpetscgrid -lpetscgsolver -lpetscfortran -lpetsccontrib \
	$PETSc_ARCH_LIBS"
PETSc_LIBS_PATHS="-L$PETSc_LIBS_PATH"
PETSc_INCLUDES_PATHS="-I$PETSc_INCLUDES_PATH -I$PETSc/bmake/$PETSc_ARCH"

# Save variables...
AC_LANG_RESTORE
LIBS=$tmpLIBS
CPPFLAGS=$tmpCPPFLAGS
AC_SUBST( PETSc )
AC_SUBST( PETSc_IMPL )
AC_SUBST( PETSc_LIBS )
AC_SUBST( PETSc_LIBS_PATH )
AC_SUBST( PETSc_LIBS_PATHS )
AC_SUBST( PETSc_INCLUDES_PATH )
AC_SUBST( PETSc_INCLUDES_PATHS )
])dnl ACX_PETSc ------------------------------------------------------------
