
dnl -------------------------------------------------------------
dnl $Id: aclocal.m4,v 1.65 2004-08-05 14:44:13 jwpeterson Exp $
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
  if test "$GXX" = yes ; then
    dnl find out the right version
    GXX_VERSION_STRING=`($CXX -v 2>&1) | grep "gcc version"`
    case "$GXX_VERSION_STRING" in
      *3.5*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.5 >>>)
  	GXX_VERSION=gcc3.5
  	;;
      *3.4*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.4 >>>)
  	GXX_VERSION=gcc3.4
  	;;
      *3.3*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.3 >>>)
  	GXX_VERSION=gcc3.3
  	;;
      *3.2*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.2 >>>)
  	GXX_VERSION=gcc3.2
  	;;
      *3.1*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.1 >>>)
  	GXX_VERSION=gcc3.1
  	;;
      *3.0*)
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
    is_ibm_xlc="`($CXX 2>&1) | egrep 'xlC'`"
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
  
        dnl Intel's ICC C++ compiler?
        is_intel_icc="`($CXX -V 2>&1) | grep 'Intel(R) C++ Compiler'`"
        if test "x$is_intel_icc" != "x" ; then
          AC_MSG_RESULT(<<< C++ compiler is Intel ICC >>>)
          GXX_VERSION=intel_icc	
        else	
  	
          dnl Intel's ECC C++ compiler for Itanium?
          is_intel_ecc="`($CXX -V 2>&1) | grep 'Intel(R) C++ Itanium(R) Compiler'`"
          if test "x$is_intel_ecc" != "x" ; then
            AC_MSG_RESULT(<<< C++ compiler is Intel Itanium ECC >>>)
            GXX_VERSION=intel_ecc
          else
  
            dnl Or Compaq's cxx compiler?
            is_dec_cxx="`($CXX -V 2>&1) | grep 'Compaq C++'`"
            if test "x$is_dec_cxx" != "x" ; then
              AC_MSG_RESULT(<<< C++ compiler is Compaq cxx >>>)
              GXX_VERSION=compaq_cxx
            else
  
  	      dnl Sun ONE Studio?
              is_sun_cc="`($CXX -V 2>&1) | grep 'Sun C++'`"
              if test "x$is_sun_cc" != "x" ; then
                AC_MSG_RESULT(<<< C++ compiler is Sun ONE Studio compiler >>>)
                GXX_VERSION=sun_studio
              else
  
  	        dnl Sun Forte?
                is_sun_forte_cc="`($CXX -V 2>&1) | grep 'Forte'`"
                if test "x$is_sun_forte_cc" != "x" ; then
                  AC_MSG_RESULT(<<< C++ compiler is Sun Forte compiler >>>)
                  GXX_VERSION=sun_forte
                else
  
  	          dnl KAI C++?
  	          is_kai_cc="`($CXX -V 2>&1) | grep 'KAI C++'`"
  	          if test "x$is_kai_cc" != "x" ; then
  	            AC_MSG_RESULT(<<< C++ compiler is KAI C++ >>>)
  	            GXX_VERSION=kai_cc
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


                        dnl  Aw, nothing suitable found...
                        AC_MSG_ERROR(Unrecognized compiler, sorry)
                        exit 1
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
dnl CXXFLAGSO  : flags for optimized mode
dnl CXXFLAGSG  : flags for debug mode
dnl CXXFLAGSS  : flags for syntax-checking mode
dnl CXXFLAGSP  : flags for profiler mode (only tested with GCC 3.2.1 but
dnl              should work with most versions.)  Note: the -a option
dnl              has been removed as of this version of GCC.
dnl
dnl Usage: SET_CXX_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(SET_CXX_FLAGS, dnl
[
  dnl Flag for creating shared objects; can be modified at a later stage
  CXXSHAREDFLAG="-shared"

  dnl Flag to add directories to the dynamic library search path; can
  dnl be changed at a later stage
  RPATHFLAG="-Wl,-rpath,"


  dnl First the flags for gcc compilers
  if test "$GXX" = yes ; then
    CXXFLAGSO="-O2 -felide-constructors -DNDEBUG"
    CXXFLAGSG=" -O2 -felide-constructors -g -ansi -pedantic -W -Wall -Wunused -Wpointer-arith -Wimplicit -Wformat -Wparentheses -Wuninitialized -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -DDEBUG"
    CXXFLAGSP="$CXXFLAGSO -g -pg"
    CXXFLAGSS="-fsyntax-only"

    CFLAGSO="-O2 -DNDEBUG"
    CFLAGSG="-g -DDEBUG"
    CFLAGSP="$CFLAGSO -g -pg"
    CFLAGSS="-fsyntax-only"

    dnl Position-independent code for shared libraries
    if test "$enableshared" = yes ; then
      CXXFLAGSO="$CXXFLAGSO -fPIC"
      CXXFLAGSG="$CXXFLAGSG -fPIC"
      CXXFLAGSP="$CXXFLAGSP -fPIC"

      CFLAGSO="$CFLAGSO -fPIC"
      CFLAGSG="$CFLAGSG -fPIC"
      CFLAGSP="$CFLAGSP -fPIC"

      dnl LDFLAGS="$LDFLAGS -fPIC"
    fi

    dnl set some flags that are specific to some versions of the
    dnl compiler:
    dnl - egcs1.1 yielded incorrect code with some loop unrolling
    dnl - after egcs1.1, the optimization flag -fstrict-aliasing was
    dnl   introduced, which enables better optimizations for
    dnl   well-written C++ code. we believe that deal.II falls into that
    dnl   category and thus enable the flag 
    dnl - egcs1.1 yielded incorrect code with vtable-thunks. thus disable
    dnl   them for egcs1.1. however, if on Linux, disabling them
    dnl   prevents programs from being linked, so take the risk of broken
    dnl   thunks on this platform
  
    case "$GXX_VERSION" in
      egcs1.1)
          case "$target" in
            *linux*)
                ;;
  
            *)
                CXXFLAGSG = "$CXXFLAGSG -fno-vtable-thunks"
                CXXFLAGSO = "$CXXFLAGSO -fno-vtable-thunks"
                CXXFLAGSP = "$CXXFLAGSP -fno-vtable-thunks"
                ;;
          esac
          ;;
  
      dnl All other gcc versions
      *)
          CXXFLAGSO="$CXXFLAGSO -funroll-loops -fstrict-aliasing"
          CFLAGSO="$CFLAGSO -funroll-loops -fstrict-aliasing"
          CXXFLAGSP="$CXXFLAGSP -funroll-loops -fstrict-aliasing"
          CFLAGSP="$CFLAGSP -funroll-loops -fstrict-aliasing"
          ;;
    esac
  
    dnl - after gcc2.95, some flags were deemed obsolete for C++
    dnl   (and are only supported for C any more), so only define them for
    dnl   previous compilers
  
    case "$GXX_VERSION" in
      egcs1.1 | gcc2.95)
          CXXFLAGSO="$CXXFLAGSO -fnonnull-objects"
          CXXFLAGSG="$CXXFLAGSG -Wmissing-declarations -Wbad-function-cast -Wtraditional -Wnested-externs"
          CXXFLAGSP="$CXXFLAGSP -fnonnull-objects"
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
          CXXFLAGSG="-DDEBUG -qmaxmem=-1 -qansialias -qrtti=all -g -qstaticinline"
          CXXFLAGSO="-DNDEBUG -O3 -qmaxmem=-1 -w -qansialias -Q=10 -qrtti=all -qstaticinline"
          CXXFLAGSP="$CXXFLAGSO -g -pg"
          CFLAGSG="-DDEBUG -qansialias -g"
          CFLAGSO="-DNDEBUG -O3 -qmaxmem=-1 -w -qansialias -Q=10"
          CFLAGSP="$CFLAGSO -g -pg"
	  CXXSHAREDFLAG="-G -qmkshrobj -bnoerrmsg"
	  CSHAREDFLAG="-G -qmkshrobj"
	  RPATHFLAG="-Qoption,link,-rpath,"
          ;;
  
      MIPSpro)
          CXXFLAGSG="-DDEBUG -LANG:std -no_auto_include -ansi -g -woff 1460"
          CXXFLAGSO="-DNDEBUG -LANG:std -no_auto_include -ansi -O2 -w"
          CFLAGSG="-DDEBUG"
          CFLAGSO="-DNDEBUG -O2 -w"

          dnl For some reason, CC forgets to add the math lib to the
          dnl linker line, so we do that ourselves
          LDFLAGS="$LDFLAGS -lm"

          dnl Position-independent code for shared libraries
          if test "$enableshared" = yes ; then
            CXXFLAGSO="$CXXFLAGSO -KPIC"
            CXXFLAGSG="$CXXFLAGSG -KPIC"
            CXXFLAGSP="$CXXFLAGSP -KPIC"

            CFLAGSO="$CFLAGSO -KPIC"
            CFLAGSG="$CFLAGSG -KPIC"
            CFLAGSP="$CFLAGSP -KPIC"

            LDFLAGS="$LDFLAGS -KPIC"
          fi
          ;;
  
      intel_icc)
          dnl Disable some warning messages:
          dnl #266: 'function declared implicitly'
          dnl       Metis function "GKfree" caused this error
          dnl       in almost every file.
          CXXFLAGSG="-Kc++eh -Krtti -O1 -w1 -DDEBUG -inline_debug_info -g -wd504"
          CXXFLAGSO="-Kc++eh -Krtti -O2 -Ob2 -DNDEBUG -tpp6 -axiMK -unroll -w0"
          CXXFLAGSP="$CXXFLAGSO -g -pg"
          CFLAGSG="-w1 -DDEBUG -inline_debug_info -wd266"
          CFLAGSO="-O2 -Ob2 -DNDEBUG -tpp6 -axiMK -unroll -w0 -vec"
          CFLAGSP="$CFLAGSO -g -pg"

          dnl Position-independent code for shared libraries
          if test "$enableshared" = yes ; then
            CXXFLAGSO="$CXXFLAGSO -KPIC"
            CXXFLAGSG="$CXXFLAGSG -KPIC"
            CXXFLAGSP="$CXXFLAGSP -KPIC"

            CFLAGSO="$CFLAGSO -KPIC"
            CFLAGSG="$CFLAGSG -KPIC"
            CFLAGSP="$CFLAGSP -KPIC"

            LDFLAGS="$LDFLAGS -KPIC"
          fi
          ;;
  
      intel_ecc)
          dnl Disable some warning messages:
          dnl #266: 'function declared implicitly'
          dnl       Metis function "GKfree" caused this error
          dnl       in almost every file.
          CXXFLAGSG="-Kc++eh -Krtti -w1 -DDEBUG -inline_debug_info -g"
          CXXFLAGSO="-Kc++eh -Krtti -O2 -DNDEBUG -unroll -w0"
          CXXFLAGSP="$CXXFLAGSO -g -pg"
          CFLAGSG="-w1 -DDEBUG -inline_debug_info -wd266"
          CFLAGSO="-O2 -DNDEBUG -unroll -w0"
          CFLAGSP="$CFLAGSO -g -pg"

          dnl Position-independent code for shared libraries
          if test "$enableshared" = yes ; then
            CXXFLAGSO="$CXXFLAGSO -KPIC"
            CXXFLAGSG="$CXXFLAGSG -KPIC"
            CXXFLAGSP="$CXXFLAGSP -KPIC"

            CFLAGSO="$CFLAGSO -KPIC"
            CFLAGSG="$CFLAGSG -KPIC"
            CFLAGSP="$CFLAGSP -KPIC"

            LDFLAGS="$LDFLAGS -KPIC"
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
  
          CXXFLAGSG="-nousing_std -nocurrent_include -model ansi -std strict_ansi -w1 -msg_display_number -timplicit_local -DDEBUG"
          CXXFLAGSO="-nousing_std -nocurrent_include -model ansi -std strict_ansi -w2 -msg_display_number -timplicit_local -DNDEBUG -O2 -fast"
          CFLAGSG="-w1 -msg_display_number -timplicit_local -DDEBUG"
          CFLAGSO="-w2 -msg_display_number -timplicit_local -DNDEBUG -O2 -fast"
  
          for i in 175 236 237 487 1136 1156 111 1182 265 ; do
            CXXFLAGSG="$CXXFLAGSG -msg_disable $i"
            CXXFLAGSO="$CXXFLAGSO -msg_disable $i"
          done
  
          dnl If we use -model ansi to compile the files, we also have to
          dnl specify it for linking
          dnl LDFLAGS="$LDFLAGS -model ansi"
  
          dnl For some reason, cxx also forgets to add the math lib to the
          dnl linker line, so we do that ourselves
          LDFLAGS="$LDFLAGS -lm"


          dnl Position-independent code for shared libraries
          if test "$enableshared" = yes ; then
            CXXFLAGSO="$CXXFLAGSO -shared"
            CXXFLAGSG="$CXXFLAGSG -shared"
            CXXFLAGSP="$CXXFLAGSP -shared"

            CFLAGSO="$CFLAGSO -shared"
            CFLAGSG="$CFLAGSG -shared"
            CFLAGSP="$CFLAGSP -shared"

            LDFLAGS="$LDFLAGS -shared"
          fi
          ;;
  
      sun_studio | sun_forte)
          CXXFLAGSG="-DDEBUG -library=stlport4 -g"
          CXXFLAGSO="-DNDEBUG -library=stlport4 -xO4"
          CFLAGSG="-DDEBUG -g"
          CFLAGSO="-DNDEBUG -xO4"

          CXXSHAREDFLAG="-G"
          CSHAREDFLAG="-G"

          dnl Linker flags & librpcsvc for XDR
          RPATHFLAG="-R"
          LIBS="-lrpcsvc $LIBS"

          dnl Position-independent code for shared libraries
          if test "$enableshared" = yes ; then
            CXXFLAGSO="$CXXFLAGSO -KPIC"
            CXXFLAGSG="$CXXFLAGSG -KPIC"
            CXXFLAGSP="$CXXFLAGSP -KPIC"

            CFLAGSO="$CFLAGSO -KPIC"
            CFLAGSG="$CFLAGSG -KPIC"
            CFLAGSP="$CFLAGSP -KPIC"

            LDFLAGS="$LDFLAGS -KPIC"
          fi
          ;;
  
      portland_group)
	  CXXFLAGSG="-g --no_using_std --instantiate=used -DDEBUG"
          CXXFLAGSO="-O2 --no_using_std --instantiate=used -DNDEBUG"
	  CFLAGSG="-g -DDEBUG"
          CFLAGSO="-O2 -DNDEBUG"

	  LDFLAGS="$LDFLAGS -fpic"
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
          CXXFLAGSG="+DA2.0W -AA +z -ext -g"
          CXXFLAGSO="+DA2.0W -AA +z -ext -O +Onolimit"
	  CFLAGSG="+DA2.0W -Aa +z -Ae -g"
          CFLAGSO="+DA2.0W -Aa +z -Ae -O +Onolimit"
          LDFLAGS="$LDFLAGS -I/usr/lib/pa20_64"
          LIBS="$LIBS -lrpcsvc"
          FLIBS="$FLIBS -lF90 -lcl -I/opt/fortran90/lib/pa20_64"
          ;;

      *)
          AC_MSG_ERROR(No compiler options for this C++ compiler
                       specified at present)
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

  else
    if (test -r $PETSC_DIR/include/petsc.h) ; then
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
      AC_SUBST(petscversion)
      AC_SUBST(petscmajorminor)
      AC_SUBST(MPI_IMPL)
  
      else
  
      dnl PETSc config failed.  Try MPI.
      enablepetsc=no
  
      dnl -------------------------------------------------------------
      if (test "$enablempi" != no) ; then
        ACX_MPI
      fi
      dnl -------------------------------------------------------------
  
    fi
  fi
  AC_SUBST(enablepetsc)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Mpi
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_MPI, 
[
  AC_CHECK_FILE($MPIHOME/include/mpi.h,
                MPI_INCLUDE_PATH=-I$MPIHOME/include)

  dnl check for libmpich
  AC_CHECK_FILE($MPIHOME/lib/libmpich.a,
                MPI_LIBRARY_PATH=$MPIHOME/lib/libmpich.a,
                MPI_LIBRARY_PATH=/mpich_bar_not_there)

  if (test -r $MPIHOME/include/mpi.h -a -r $MPI_LIBRARY_PATH) ; then
    AC_SUBST(MPI_INCLUDE_PATH)
    dnl Here is a little hack.  If MPI_LIBRARY_PATH is valid 
    dnl for libmpich, we assume it will also be available for
    dnl libpmpich, which we will also require.  Therefore we
    dnl deftly change MPI_LIBRARY_PATH to link with both
    dnl mpich and pmpich.  This has no hope of working if
    dnl you are using some sort of specialized mpi instead
    dnl of mpich.
    MPI_LIBRARY_PATH="-L $MPIHOME/lib -lmpich -lpmpich"
    AC_SUBST(MPI_LIBRARY_PATH)
    AC_DEFINE(HAVE_MPI, 1,
	      [Flag indicating whether or not MPI is available])
    AC_MSG_RESULT(<<< Configuring library with MPI support >>>)
  else
    enablempi=no
  fi

  AC_SUBST(enablempi)	
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Space Filling Curves
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_SFC, 
[
  AC_CHECK_FILE(./contrib/sfcurves/sfcurves.h,
                [
                  SFC_INCLUDE=-I$PWD/contrib/sfcurves
                  SFC_LIB="\$(EXTERNAL_LIBDIR)/libsfcurves\$(EXTERNAL_LIBEXT)"
                  AC_SUBST(SFC_INCLUDE)
                  AC_SUBST(SFC_LIB)
                  AC_DEFINE(HAVE_SFCURVES, 1,
                             [Flag indicating whether or not Space filling curves are available])
                  AC_MSG_RESULT(<<< Configuring library with SFC support >>>)
		  enablesfc=yes
                ],
                [enablesfc=no])
])
dnl -------------------------------------------------------------




dnl -------------------------------------------------------------
dnl Read/Write Compressed Streams with gzstream
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_GZ, 
[
  AC_CHECK_HEADERS(zlib.h, have_zlib_h=yes)
  AC_CHECK_LIB(z, gzopen, have_libz=yes)
  if (test "$have_zlib_h" = yes \
        -a "$have_libz"   = yes) ; then
     AC_CHECK_FILE(./contrib/gzstream/gzstream.h,
                   [
                     GZSTREAM_INCLUDE=-I$PWD/contrib/gzstream
                     GZSTREAM_LIB="\$(EXTERNAL_LIBDIR)/libgzstream\$(EXTERNAL_LIBEXT) -lz"
                     AC_SUBST(GZSTREAM_INCLUDE)
                     AC_SUBST(GZSTREAM_LIB)
                     AC_DEFINE(HAVE_GZSTREAM, 1,
                                [Flag indicating whether or not gzstreams are available])
                     AC_MSG_RESULT(<<< Configuring library with gzstreams support >>>)
		     enablegz=yes                     
                   ],
                   [enablegz=no])
  fi
])
dnl -------------------------------------------------------------




dnl -------------------------------------------------------------
dnl LASPACK Iterative Solvers
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_LASPACK, 
[
  AC_CHECK_FILE(./contrib/laspack/lastypes.h,
		[
                  LASPACK_INCLUDE_PATH=$PWD/contrib/laspack
                  LASPACK_INCLUDE=-I$LASPACK_INCLUDE_PATH
                  LASPACK_LIB="\$(EXTERNAL_LIBDIR)/liblaspack\$(EXTERNAL_LIBEXT)"
                  AC_SUBST(LASPACK_INCLUDE)
                  AC_SUBST(LASPACK_LIB)
                  AC_DEFINE(HAVE_LASPACK, 1,
                            [Flag indicating whether or not LASPACK iterative solvers are available])
                  laspack_version=`grep "define LASPACK_VERSION " $LASPACK_INCLUDE_PATH/version.h | sed -e "s/[[^0-9.]]*//g"`
                  AC_MSG_RESULT(<<< Configuring library with LASPACK version $laspack_version support >>>)
                ],
                [])
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
                  METIS_LIB="\$(EXTERNAL_LIBDIR)/libmetis\$(EXTERNAL_LIBEXT)"
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
                      PARMETIS_LIB="\$(EXTERNAL_LIBDIR)/libparmetis\$(EXTERNAL_LIBEXT)"
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
              AC_HELP_STRING([--with-tecplot=PATH],[Specify the path where the Tecplot is installed]),
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
dnl if TetGen is used we need the header path and the lib
dnl all necessary files are believed to be in environment 
dnl variable $TETGEN_DIR
  if test $TETGEN_DIR; then
     TETGEN_VERSION="-DTETGEN_13"
     TETGEN_INCLUDE="-I$TETGEN_DIR"
     TETGEN_LIBRARY="$TETGEN_DIR/libtet.a"
  else
     TETGEN_VERSION=""
     TETGEN_INCLUDE=""
     TETGEN_LIBRARY=""
     enabletetgen=no
  fi
dnl TetGen version 1.3:
	AC_SUBST(TETGEN_VERSION)
	AC_SUBST(TETGEN_INCLUDE)
	AC_SUBST(TETGEN_LIBRARY)	
	AC_SUBST(enabletetgen)
	AC_DEFINE(HAVE_TETGEN, 1, [Flag indicating whether the library will be compiled with TetGen support])
        AC_MSG_RESULT(<<< Configuring library with TetGen support >>>)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl netCDF
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_NETCDF,
[
  AC_CHECK_FILE(./contrib/netcdf/lib/$host/libnetcdf.a,
		NETCDF_LIB=$PWD/contrib/netcdf/lib/$host/libnetcdf.a)
  AC_CHECK_FILE(./contrib/netcdf/include/netcdf.h,
		NETCDF_INCLUDE_PATH=$PWD/contrib/netcdf/include)

  if (test -r $NETCDF_INCLUDE_PATH/netcdf.h -a "x$NETCDF_LIB" != x) ; then
    NETCDF_INCLUDE=-I$NETCDF_INCLUDE_PATH
    AC_SUBST(NETCDF_LIB)
    AC_SUBST(NETCDF_INCLUDE)
    AC_DEFINE(HAVE_NETCDF, 1,
              [Flag indicating whether the library shall be compiled to support netcdf files])
    AC_MSG_RESULT(<<< Configuring library with netCDF support >>>)
    have_netcdf=yes
  fi
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl ExodusII 
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_EXODUS,
[
  AC_CHECK_FILE(./contrib/exodus/lib/$host/libexoIIv2c.a,
		EXODUS_LIB=$PWD/contrib/exodus/lib/$host/libexoIIv2c.a)
  AC_CHECK_FILE(./contrib/exodus/include/exodusII.h,
		EXODUS_INCLUDE_PATH=$PWD/contrib/exodus/include)

  if (test -r $EXODUS_INCLUDE_PATH/exodusII.h -a "x$EXODUS_LIB" != x) ; then
    EXODUS_INCLUDE=-I$EXODUS_INCLUDE_PATH
    AC_SUBST(EXODUS_LIB)
    AC_SUBST(EXODUS_INCLUDE)
    AC_DEFINE(HAVE_EXODUS_API, 1,
	      [Flag indicating whether the library shall be compiled to use the Exodus interface])
    AC_MSG_RESULT(<<< Configuring library with Exodus API support >>>)
  fi
])
dnl -------------------------------------------------------------






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
                       AC_MSG_RESULT([Found valid MPICH installlaion...])
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
