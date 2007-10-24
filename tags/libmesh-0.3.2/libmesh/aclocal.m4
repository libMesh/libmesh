
dnl -------------------------------------------------------------
dnl $Id: aclocal.m4,v 1.15 2003-02-24 14:35:52 benkirk Exp $
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
      *"egcs-1.1"*)
  	AC_MSG_RESULT(<<< C++ compiler is egcs-1.1 >>>)
  	GXX_VERSION=egcs1.1
  	;;
      *2.95*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-2.95 >>>)
  	GXX_VERSION=gcc2.95
	AC_DEFINE(BROKEN_IOSTREAM, 1,
             [This compiler is known not to support some iostream
              functionality])
         AC_MSG_RESULT(<<< Configuring library for broken iostream >>>)
  	;;
      *2.96*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-2.96 >>>)
  	GXX_VERSION=gcc2.96
	AC_DEFINE(BROKEN_IOSTREAM, 1,
             [This compiler is known not to support some iostream
              functionality])
         AC_MSG_RESULT(<<< Configuring library for broken iostream >>>)
  	;;
      *2.97*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-2.97 >>>)
  	GXX_VERSION=gcc2.97
  	;;
      *3.0*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.0 >>>)
  	GXX_VERSION=gcc3.0
  	;;
      *3.1*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.1 >>>)
  	GXX_VERSION=gcc3.1
  	;;
      *3.2*)
  	AC_MSG_RESULT(<<< C++ compiler is gcc-3.2 >>>)
  	GXX_VERSION=gcc3.2
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
    is_ibm_xlc="`($CXX 2>&1) | egrep 'VisualAge C++|C Set ++|C for AIX Compiler'`"
    if test "x$is_ibm_xlc" != "x"  ; then
      dnl Ah, this is IBM's C++ compiler. Unfortunately, we don't presently
      dnl know how to check the version number, so assume that is sufficiently
      dnl high...
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
  	
          dnl Intel's ICC C++ compiler for Itanium?
          is_intel_ecc="`($CXX -V 2>&1) | grep 'Intel(R) C++ Itanium(TM) Compiler'`"
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
  
  	      dnl Sun Workshop?
              is_sun_cc="`($CXX -V 2>&1) | grep 'Sun WorkShop'`"
              if test "x$is_sun_cc" != "x" ; then
                AC_MSG_RESULT(<<< C++ compiler is Sun Workshop compiler >>>)
                GXX_VERSION=sun_workshop
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
dnl CXXDEPFLAG : flag for creating dependency info
dnl
dnl Usage: SET_CXX_FLAGS
dnl
dnl -------------------------------------------------------------
AC_DEFUN(SET_CXX_FLAGS, dnl
[
  dnl Flag for creating dependencies; can be modified at a later stage
  CXXDEPFLAG="-MM"

  dnl Flag for creating shared objects; can be modified at a later stage
  CXXSHAREDFLAG="-shared"

  dnl First the flags for gcc compilers
  if test "$GXX" = yes ; then
    CXXFLAGSO="-O2 -felide-constructors -DNDEBUG"
    CXXFLAGSG="-g -ansi -pedantic -W -Wall -Wunused -Wpointer-arith -Wimplicit -Wformat -Wparentheses -O -Wuninitialized -DDEBUG"
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

      LDFLAGS="$LDFLAGS -fPIC"
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
  
  
  else
    dnl Non-gcc compilers
  
    case "$GXX_VERSION" in
      ibm_xlc)
          CXXFLAGSG="-DDEBUG -check=bounds -info=all -qrtti=all"
          CXXFLAGSO="-DNDEBUG -O3 -qmaxmem=-1 -w -qansialias -qrtti=all -Q"
          CXXFLAGSP="$CXXFLAGSO -g -pg"
          CFLAGSG="-DDEBUG -check=bounds -info=all -qrtti=all"
          CFLAGSO="-DNDEBUG -O2 -w -qansialias"
          CFLAGSP="$CFLAGSO -g -pg"
	  CXXSHAREDFLAG="-qmkshrobj"
	  CXXDEPFLAG="-qmakedep"
          ;;
  
      MIPSpro)
          CXXFLAGSG="-DDEBUG -LANG:std -no_auto_include -ansi -g -woff 1460"
          CXXFLAGSO="-DNDEBUG -LANG:std -no_auto_include -ansi -O2 -w"
          CFLAGSG="-DDEBUG"
          CFLAGSO="-DNDEBUG -O2 -w"
          CXXDEPFLAG="-M"

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
          CXXFLAGSG="-Kc++eh -Krtti -w1 -DDEBUG -inline_debug_info -g -wd504"
          CXXFLAGSO="-Kc++eh -Krtti -O2 -Ob2 -DNDEBUG -tpp6 -axiMK -unroll -w0 -vec"
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
          CFLAGSO="-O2 -DNDEBUG -axiMK -unroll -w0"
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
          LDFLAGS="$LDFLAGS -model ansi"
  
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
  
      sun_workshop | sun_forte)
          CXXFLAGSG="-DDEBUG -w"
          CXXFLAGSO="-DNDEBUG -w"
          CFLAGSG="-DDEBUG -w"
          CFLAGSO="-DNDEBUG -w"

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
	  CXXFLAGSG="-A -Xa -DDEBUG"
          CXXFLAGSO="-A -Xa -DNDEBUG"
	  CFLAGSG="-A -Xa -DDEBUG"
          CFLAGSO="-A -Xa -DNDEBUG"
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
          dnl override value for make dep's
          CXXDEPFLAG="-E +m"
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

  if (test -r $PETSC_DIR/include/petsc.h) ; then
    AC_PROG_F77            dnl Petsc requires linking with FORTRAN libraries 
    AC_F77_LIBRARY_LDFLAGS
    AC_SUBST(PETSC_DIR)
    AC_DEFINE(HAVE_PETSC, 1,
	      [Flag indicating whether or not Petsc is available])
    AC_DEFINE(HAVE_MPI, 1,
	      [Flag indicating whether or not MPI is available])
    petscmajor=`grep "define PETSC_VERSION_MAJOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MAJOR[ ]*//g"`
    petscminor=`grep "define PETSC_VERSION_MINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_MINOR[ ]*//g"`
    petscsubminor=`grep "define PETSC_VERSION_SUBMINOR" $PETSC_DIR/include/petscversion.h | sed -e "s/#define PETSC_VERSION_SUBMINOR[ ]*//g"`
    petscversion=$petscmajor.$petscminor.$petscsubminor
    AC_MSG_RESULT(<<< Configuring library with PETSc version $petscversion support >>>)
    AC_SUBST(petscversion)
  else
    enablepetsc=no  
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
                MPI_INCLUDE_PATH=$MPIHOME/include)

  AC_CHECK_FILE($MPIHOME/lib/libmpich.a,
                MPI_LIBRARY_PATH=$MPIHOME/lib/libmpich.a,
                MPI_LIBRARY_PATH=/mpich_bar_not_there)

  if (test -r $MPI_INCLUDE_PATH/mpi.h -a -r $MPI_LIBRARY_PATH) ; then
    AC_SUBST(MPI_INCLUDE_PATH)
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
                SFC_INCLUDE_PATH=$PWD/contrib/sfcurves)

  if (test -r $SFC_INCLUDE_PATH/sfcurves.h) ; then
    AC_SUBST(SFC_INCLUDE_PATH)
    AC_DEFINE(HAVE_SFCURVES, 1,
              [Flag indicating whether or not Space filling curves are available])
    AC_MSG_RESULT(<<< Configuring library with SFC support >>>)
  else
    enablesfc=no
  fi

  AC_SUBST(enablesfc)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl LASPACK Iterative Solvers
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_LASPACK, 
[
  AC_CHECK_FILE(./contrib/laspack/lastypes.h,
                LASPACK_INCLUDE_PATH=$PWD/contrib/laspack)

  if (test -r $LASPACK_INCLUDE_PATH/lastypes.h) ; then
    AC_SUBST(LASPACK_INCLUDE_PATH)
    AC_DEFINE(HAVE_LASPACK, 1,
              [Flag indicating whether or not LASPACK iterative solvers are available])
    laspack_version=`grep "define LASPACK_VERSION " $LASPACK_INCLUDE_PATH/version.h | sed -e "s/[[^0-9.]]*//g"`
    AC_MSG_RESULT(<<< Configuring library with LASPACK version $laspack_version support >>>)
  else
    enablelaspack=no
  fi

  AC_SUBST(enablelaspack)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Metis
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_METIS, 
[
  AC_CHECK_FILE(./contrib/metis/Lib/metis.h,
                METIS_INCLUDE_PATH=$PWD/contrib/metis/Lib)

  if (test -r $METIS_INCLUDE_PATH/metis.h) ; then
    AC_SUBST(METIS_INCLUDE_PATH)
    AC_DEFINE(HAVE_METIS, 1,
	      [Flag indicating whether or not Metis is available])
    AC_MSG_RESULT(<<< Configuring library with Metis support >>>)
  else
    enablemetis=no
  fi

  AC_SUBST(enablemetis)	
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl Tecplot
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_TECPLOT,
[
  AC_ARG_WITH(tecplot,
      [  --with-tecplot=[PATH] Specify the path where the Tecplot is installed],
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
    AC_SUBST(TECPLOT_LIBRARY)
    AC_SUBST(TECPLOT_INCLUDE_PATH)
    AC_DEFINE(HAVE_TECPLOT_API, 1,
              [Flag indicating whether the library shall be compiled to use the Tecplot interface])
    AC_MSG_RESULT(<<< Configuring library with Tecplot API support >>>)
  else
    enabletecplot=no
  fi

  AC_SUBST(enabletecplot)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl netCDF
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_NETCDF,
[
  AC_CHECK_FILE(./contrib/netcdf/lib/$host/libnetcdf.a,
		LIBNETCDF_PATH=$PWD/contrib/netcdf/lib/$host)
  AC_CHECK_FILE(./contrib/netcdf/include/netcdf.h,
		NETCDF_INCLUDE_PATH=$PWD/contrib/netcdf/include)

  if (test -r $LIBNETCDF_PATH/libnetcdf.a -a -r $NETCDF_INCLUDE_PATH/netcdf.h) ; then
    LIBNETCDF=$LIBNETCDF_PATH/libnetcdf.a
    AC_SUBST(LIBNETCDF)
    AC_SUBST(NETCDF_INCLUDE_PATH)
    AC_DEFINE(HAVE_NETCDF, 1,
              [Flag indicating whether the library shall be compiled to support netcdf files])
    AC_MSG_RESULT(<<< Configuring library with netCDF support >>>)
  else
    enablenetcdf=no
  fi

  AC_SUBST(enablenetcdf)
])
dnl -------------------------------------------------------------





dnl -------------------------------------------------------------
dnl HDF4
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_HDF4,
[
  AC_CHECK_FILE(./contrib/hdf4/include/mfhdf.h,
		HDF4_INCLUDE_PATH=$PWD/contrib/hdf4/include)

  dnl This isn't the world's best test, it
  dnl assumes that if the mfhdf library is there,
  dnl they are all there...
  AC_CHECK_FILE(./contrib/hdf4/lib/$host/libmfhdf.a,
		HDF4_LIB_PATH=$PWD/contrib/hdf4/lib/$host)
  HDF4_LIB_NAMES="-lmfhdf -ldf -ljpeg -lz"

  if (test -r $HDF4_LIB_PATH/libmfhdf.a -a -r $HDF4_INCLUDE_PATH/mfhdf.h) ; then
    AC_SUBST(HDF4_INCLUDE_PATH)
    AC_SUBST(HDF4_LIB_PATH)
    AC_SUBST(HDF4_LIB_NAMES)
    AC_DEFINE(HAVE_HDF4, 1,
	      [Flag indicating whether the library shall be compiled to support hdf4 files])
    AC_MSG_RESULT(<<< Configuring library with HDF support >>>)
  else
    enablehdf4=no
  fi

  AC_SUBST(enablehdf4)
])
dnl -------------------------------------------------------------



dnl -------------------------------------------------------------
dnl ExodusII 
dnl -------------------------------------------------------------
AC_DEFUN(CONFIGURE_EXODUS,
[
  AC_CHECK_FILE(./contrib/exodus/lib/$host/libexoIIv2c.a,
		EXODUS_EXOII_PATH=$PWD/contrib/exodus/lib/$host)
  AC_CHECK_FILE(./contrib/exodus/include/exodusII.h,
		EXODUS_INCLUDE_PATH=$PWD/contrib/exodus/include)

  if (test -r $EXODUS_EXOII_PATH/libexoIIv2c.a -a -r $EXODUS_INCLUDE_PATH/exodusII.h) ; then
    LIBEXOII=$EXODUS_EXOII_PATH/libexoIIv2c.a
    AC_SUBST(LIBEXOII)
    AC_SUBST(EXODUS_INCLUDE_PATH)
    AC_DEFINE(HAVE_EXODUS_API, 1,
	      [Flag indicating whether the library shall be compiled to use the Exodus interface])
    AC_MSG_RESULT(<<< Configuring library with Exodus API support >>>)
  else
    enableexodus=no
  fi

  AC_SUBST(enableexodus)
])
dnl -------------------------------------------------------------
