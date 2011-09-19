dnl -------------------------------------------------------------
dnl $Id$
dnl -------------------------------------------------------------

dnl -------------------------------------------------------------
dnl Determine the C++ compiler in use. Return the name and possibly
dnl version of this compiler in GXX_VERSION.
dnl
dnl Usage: DETERMINE_CXX_BRAND
dnl
dnl -------------------------------------------------------------
AC_DEFUN([DETERMINE_CXX_BRAND],
[
  dnl First check for gcc version, avoids intel's icc from
  dnl pretending to be gcc
  REAL_GXX=`($CXX -v 2>&1) | grep "gcc version"`

  dnl Intel's v12.1 does this:
  dnl $ icpc -v
  dnl   icpc version 12.1.0 (gcc version 4.4.4 compatibility)
  dnl cath that and do not interpret it as 'REAL_GXX' compiler
  is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R)' | grep 'Compiler'`"
  if test "x$is_intel_icc" != "x" ; then
    REAL_GXX=""
  fi	
  
  if (test "$GXX" = yes -a "x$REAL_GXX" != "x" ) ; then
    dnl find out the right version
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "gcc version"`
    case "$GXX_VERSION_STRING" in
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
    dnl Check for Apple compilers
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
    dnl Check other (non-gcc) compilers
  
    dnl Check for IBM xlC. For some reasons, depending on some environment
    dnl variables, moon position, and other reasons unknown to me, the
    dnl compiler displays different names in the first line of output, so
    dnl check various possibilities.  Calling xlC with no arguments displays
    dnl the man page.  Grepping for case-sensitive xlc is not enough if the
    dnl user wants xlC, so we used case-insensitive grep instead.
    is_ibm_xlc="`($CXX 2>&1) | egrep -i 'xlc'`"
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
AC_DEFUN([SET_CXX_FLAGS], dnl
[
  dnl Flag for creating shared objects; can be modified at a later stage
  case "$target_os" in
    *darwin*)
      if test "$enableshared" = yes ; then
        CXXFLAGS_OPT="-fno-common"
        CXXFLAGS_DVL="-fno-common"
        CXXFLAGS_DBG="-fno-common"
        LDFLAGS="$LDFLAGS -Wl,-undefined,dynamic_lookup,-flat_namespace"
        if test $APPLE_GCC = true ; then
          case "$GXX_VERSION_STRING" in
            *4.0.* | *3.4.* | *3.3.* | *3.2.* | *3.1.* | *3.0.* | *2.97* | *2.96* | *2.95* | *"egcs-1.1"*)
              CXXSHAREDFLAG="-dynamiclib -Wl,-undefined,dynamic_lookup,-flat_namespace"
              ;;
            *)
              CXXSHAREDFLAG="-ldylib1.o -dynamiclib -Wl,-undefined,dynamic_lookup,-flat_namespace,-no_compact_linkedit"
              ;;
            *)
          esac
        else
          CXXSHAREDFLAG="-dynamiclib -Wl,-undefined,dynamic_lookup,-flat_namespace"
        fi
        CSHAREDFLAG="-dynamiclib -Wl,-undefined,dynamic_lookup,-flat_namespace"
      fi
      ;;
    *)  
      CXXSHAREDFLAG="-shared"
      ;;
  esac

  dnl Flag to add directories to the dynamic library search path; can
  dnl be changed at a later stage
  RPATHFLAG="-Wl,-rpath,"

  dnl Flag for profiling mode; can me modified at a later stage
  PROFILING_FLAGS="-pg"

  dnl The -g flag is necessary for OProfile to produce annotations
  dnl -fno-omit-frame-pointer flag turns off an optimization that
  dnl interferes with OProfile callgraphs
  OPROFILE_FLAGS="-g -fno-omit-frame-pointer"

  dnl First the flags for gcc compilers
  if (test "$GXX" = yes -a "x$REAL_GXX" != "x" ) ; then
    CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -felide-constructors"
    CXXFLAGS_DVL="$CXXFLAGS_DVL -O2 -felide-constructors -g -ansi -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses -Wuninitialized"
    CXXFLAGS_DBG="$CXXFLAGS_DBG -O0 -felide-constructors -g -ansi -pedantic -W -Wall -Wextra -Wno-long-long -Wunused -Wpointer-arith -Wformat -Wparentheses"

    CFLAGS_OPT="-O2"
    CFLAGS_DVL="$CFLAGS_OPT -g -Wimplicit"
    CFLAGS_DBG="-g -Wimplicit"

    dnl Position-independent code for shared libraries
    if test "$enableshared" = yes ; then
      CXXFLAGS_OPT="$CXXFLAGS_OPT -fPIC"
      CXXFLAGS_DVL="$CXXFLAGS_DVL -fPIC"
      CXXFLAGS_DBG="$CXXFLAGS_DBG -fPIC"

      CFLAGS_OPT="$CFLAGS_OPT -fPIC"
      CFLAGS_DVL="$CFLAGS_DVL -fPIC"
      CFLAGS_DBG="$CFLAGS_DBG -fPIC"

      FFLAGS="$FFLAGS -fPIC"
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

         if test $APPLE_GCC = true ; then
           CXXFLAGS_DBG="$CXXFLAGS_DBG -std=c++0x -Woverloaded-virtual"
         else
           CXXFLAGS_DBG="$CXXFLAGS_DBG -std=c++0x -Woverloaded-virtual -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"
         fi
	 ;;

      gcc3.* | gcc4.*)
	 CXXFLAGS_OPT="$CXXFLAGS_OPT -Wdisabled-optimization"
         CXXFLAGS_DVL="$CXXFLAGS_DVL -Woverloaded-virtual -Wdisabled-optimization"
         
         if test $APPLE_GCC = true ; then
	   CXXFLAGS_DBG="$CXXFLAGS_DBG -Woverloaded-virtual"
	 else    
	   CXXFLAGS_DBG="$CXXFLAGS_DBG -Woverloaded-virtual -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"	
	 fi
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

            FFLAGS="$FFLAGS -KPIC"

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
          intel_icc_v10.1 | intel_icc_v11.x | intel_icc_v12.x)
              dnl Disable some warning messages:
              dnl #175: 'subscript out of range'
              dnl       FIN-S application code causes many false
              dnl       positives with this
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
              PROFILING_FLAGS="-p"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -w1 -g -wd175 -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -O3 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DVL -w1 -g -wd175 -wd1476 -wd1505 -wd1572"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -g -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O3 -unroll -w0 -ftz -par_report0 -openmp_report0"
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

              PROFILING_FLAGS="-p"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -Kc++eh -Krtti -O1 -w1 -g -wd504 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -Kc++eh -Krtti -O2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -g -inline_debug_info -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
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

              CXXFLAGS_DBG="$CXXFLAGS_DBG -Kc++eh -Krtti -O1 -w1 -g -wd504 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -Kc++eh -Krtti -O2 -Ob2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -g -inline_debug_info -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O2 -Ob2 $INTEL_AX_FLAG -unroll -w0 -vec_report0 -par_report0 -openmp_report0"
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
              CXXFLAGS_DBG="$CXXFLAGS_DBG -w1 -inline_debug_info -g -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -inline_debug_info -g -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
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
              CXXFLAGS_DBG="$CXXFLAGS_DBG -Kc++eh -Krtti -w1 -inline_debug_info -g -wd1476 -wd1505 -wd1572"
              CXXFLAGS_OPT="$CXXFLAGS_OPT -Kc++eh -Krtti -O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
              CXXFLAGS_DVL="$CXXFLAGS_DBG"
              CFLAGS_DBG="$CFLAGS_DBG -w1 -inline_debug_info -g -wd266 -wd1572"
              CFLAGS_OPT="$CFLAGS_OPT -O2 -unroll -w0 -ftz -par_report0 -openmp_report0"
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
              CFLAGS_DBG="-w1 -inline_debug_info -g -wd266"
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
            intel_*_v1?.*)
              CXXFLAGS_OPT="$CXXFLAGS_OPT -fPIC"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -fPIC"
              CXXFLAGS_DVL="$CXXFLAGS_DVL -fPIC"
        
              CFLAGS_OPT="$CFLAGS_OPT -fPIC"
              CFLAGS_DBG="$CFLAGS_DBG -fPIC"
              CFLAGS_DVL="$CFLAGS_DVL -fPIC"

              FFLAGS="$FFLAGS -fPIC"

              LDFLAGS="$LDFLAGS -fPIC"
	      ;;
	    *)
              CXXFLAGS_OPT="$CXXFLAGS_OPT -KPIC"
              CXXFLAGS_DBG="$CXXFLAGS_DBG -KPIC"
              CXXFLAGS_DVL="$CXXFLAGS_DVL -KPIC"
        
              CFLAGS_OPT="$CFLAGS_OPT -KPIC"
              CFLAGS_DBG="$CFLAGS_DBG -KPIC"
              CFLAGS_DVL="$CFLAGS_DVL -KPIC"

	      FFLAGS="$FFLAGS -KPIC"

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

	    FFLAGS="$FFLAGS -KPIC"

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
