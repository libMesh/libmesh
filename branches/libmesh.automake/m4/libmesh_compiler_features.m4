# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_COMPILER_FEATURES],
[
AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(------- Configuring compiler features -------)
AC_MSG_RESULT(---------------------------------------------)

# --------------------------------------------------------------
# Real precision - double by default
# --------------------------------------------------------------
ACX_CHOOSE_PRECISION

# --------------------------------------------------------------
# Determine support for the C99 "restrict" keyword in the C++
# compiler.  This should place a pound-define statement in the
# config file along the lines of:
#
# #define _libmesh_restrict __restrict
#
# The restrict keyword allows the compiler to make
# certain optimizations which require non-overlapping pointer
# addresses, i.e. in a call to:
#
# foo(int* restrict a, int* restrict b)
#
# the compiler may assume that the memory addresses pointed to
# by a and b do not overlap, and therefore certain low-level
# optimizations can be made.  From the Autoconf documentation:
#
# If the C compiler recognizes a variant spelling for the restrict
# keyword (__restrict, __restrict__, or _Restrict), then define restrict
# to that; this is more likely to do the right thing with compilers that
# support language variants where plain restrict is not a
# keyword. Otherwise, if the C compiler recognizes the restrict keyword,
# don't do anything. Otherwise, define restrict to be empty. Thus,
# programs may simply use restrict as if every C compiler supported it;
# for those that do not, the makefile or configuration header defines it
# away.
# 
# Although support in C++ for the restrict keyword is not required,
# several C++ compilers do accept the keyword. This macro works for
# them, too.
# 
# This macro caches 'no' in the ac_cv_c_restrict variable if restrict is
# not supported, and a supported spelling otherwise.
AC_C_RESTRICT


# --------------------------------------------------------------
# getpwuid - enabled by default
# Some systems, for example Crays, actually have getpwuid on the head-node
# but (if I understand correctly) a dynamically-linked glibc is not available
# on the backend, which is running a reduced operating system like Compute
# Node Linux.  Thus functions like getpwuid cannot be called.  This makes
# automatically testing for the existence of getpwuid on the login node
# difficult.  If you know that you are on such a system, configure with
# --disable-getpwuid.
# --------------------------------------------------------------
AC_ARG_ENABLE(getpwuid,
              AC_HELP_STRING([--enable-getpwuid],
                             [allow calls to getpwuid]),
              enablegetpwuid=$enableval,
              enablegetpwuid=yes)

if test "$enablegetpwuid" != no ; then
  AC_DEFINE(HAVE_GETPWUID, 1,
           [Flag indicating if the library should be built with calls to getpwuid()])
  AC_MSG_RESULT(<<< Configuring library with getpwuid >>>)
fi
# --------------------------------------------------------------



# --------------------------------------------------------------
# C++ exceptions - enabled by default
# --------------------------------------------------------------
AC_ARG_ENABLE(exceptions,
              AC_HELP_STRING([--enable-exceptions],
                             [throw exceptions on unexpected errors]),
              enableexceptions=$enableval,
              enableexceptions=yes)

if test "$enableexceptions" != no ; then
  AC_DEFINE(ENABLE_EXCEPTIONS, 1,
           [Flag indicating if the library should be built to throw C++ exceptions on unexpected errors])
  AC_MSG_RESULT(<<< Configuring library with exception throwing support >>>)
fi
# --------------------------------------------------------------



# --------------------------------------------------------------
# Check for important type sizes
# --------------------------------------------------------------
AC_CHECK_SIZEOF(short int)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(long int)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(void *)
# Check the size of a function pointer.  This is the same
# as a void* on most systems which matter (POSIX).
# It turns out that AC_CHECK_SIZEOF can't do this
# automatically... so we use a kind of hack with the 3rd
# argument to do it.
AC_CHECK_SIZEOF([function_pointer], [], [typedef void (*function_pointer)();])
# AC_CHECK_SIZEOF([void(*)(void)]) <-- Does not work! 

# --------------------------------------------------------------
# Check for Run Time Type Identification support
# --------------------------------------------------------------
AC_CXX_RTTI



# --------------------------------------------------------------
# Check for headers
# --------------------------------------------------------------
AC_CHECK_HEADERS(getopt.h)
AC_CHECK_HEADERS(csignal)
AC_CHECK_HEADERS(sys/resource.h)
AC_CXX_HAVE_LOCALE
AC_CXX_HAVE_SSTREAM

AC_ARG_ENABLE(unordered-containers,
              AC_HELP_STRING([--enable-unordered-containers],
                             [Use unordered_map/unordered_set if available]),
              enableunorderedcontainers=$enableval,
              enableunorderedcontainers=yes)

  if test "$enableunorderedcontainers" != no ; then
    # The following routines, defined in unordered.m4, check to see if the compiler can compile programs using
    # various quasi-standard hash containers.
    ACX_BEST_UNORDERED_MULTIMAP
    ACX_BEST_UNORDERED_MAP
    ACX_BEST_UNORDERED_SET
  else
    ACX_STD_MAP
    ACX_STD_MULTIMAP
    ACX_STD_SET
  fi

AC_CHECK_HEADERS(dlfcn.h)
AX_CXX_GCC_ABI_DEMANGLE
AX_CXX_GLIBC_BACKTRACE





# -------------------------------------------------------------
# OpenMP Support  -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(openmp,
             AC_HELP_STRING([--enable-openmp],
                            [Build with OpenMP Support]),
             enableopenmp=$enableval,
             enableopenmp=yes)
if (test "$enableopenmp" != no) ; then
   AX_OPENMP
   #The above call only sets the flag for C++
   if (test "x$OPENMP_CXXFLAGS" != x) ; then
     AC_MSG_RESULT(<<< Configuring library with OpenMP support >>>)
   fi
   OPENMP_CFLAGS=$OPENMP_CXXFLAGS
   OPENMP_FFLAGS=$OPENMP_CXXFLAGS
   AC_SUBST(OPENMP_CXXFLAGS)
   AC_SUBST(OPENMP_CFLAGS)
   AC_SUBST(OPENMP_FFLAGS)
   AC_SUBST(enableopenmp)
fi
# -------------------------------------------------------------

])



# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_CONFIGURE_FEATURES],
[
# -------------------------------------------------------------
# size of subdomain_id -- default 2
# -------------------------------------------------------------
AC_ARG_WITH([subdomain_id_bytes],
	    AC_HELP_STRING([--with-subdomain-id-bytes=<1|2|4>],
                           [bytes of storage per element used to store the subdomain_id]),
	    [subdomain_bytes="$withval"],
	    [subdomain_bytes=2])

if test "$subdomain_bytes" == 1 ; then
  AC_DEFINE(SUBDOMAIN_ID_BYTES, 1,
           [size of subdomain_id])
  AC_MSG_RESULT(configuring size of subdomain_id... 1)
elif test "$subdomain_bytes" == 2 ; then
  AC_DEFINE(SUBDOMAIN_ID_BYTES, 2,
           [size of subdomain_id])
  AC_MSG_RESULT(configuring size of subdomain_id... 2)
elif test "$subdomain_bytes" == 4 ; then
  AC_DEFINE(SUBDOMAIN_ID_BYTES, 4,
           [size of subdomain_id])
  AC_MSG_RESULT(configuring size of subdomain_id... 4)
else
  AC_MSG_RESULT(unrecognized subdomain_id size - configuring size...2)
fi

# -------------------------------------------------------------



# -------------------------------------------------------------
# Allow user to specify --enable-everything
#
# This flag will cause all (non-conflicting) options to be
# enabled for the purposes of configuration.  For example, by
# performance logging is off by default, however
# --enable-everything will change it to be on by default.
#
# Note specific flags will override --enable-everything for
# that particular package, i.e.
#  ./configure --enable-everything --disable-perflog
#
# -------------------------------------------------------------
AC_ARG_ENABLE(everything,
              AC_HELP_STRING([--enable-everything],
                             [treat all applicable options as enabled]),
              enableeverything=$enableval,
              enableeverything=no)


# --------------------------------------------------------------
# Write stack trace output files on error() - disabled by default
# --------------------------------------------------------------
AC_ARG_ENABLE(tracefiles,
              AC_HELP_STRING([--enable-tracefiles],
                             [write stack trace files on unexpected errors]),
              enabletracefiles=$enableval,
              enabletracefiles=$enableeverything)

if test "$enabletracefiles" != no ; then
  AC_DEFINE(ENABLE_TRACEFILES, 1,
           [Flag indicating if the library should be built to write stack trace files on unexpected errors])
  AC_MSG_RESULT(<<< Configuring library with stack trace file support >>>)
fi
# --------------------------------------------------------------


# -------------------------------------------------------------
# AMR -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(amr,
              AC_HELP_STRING([--enable-amr],
                             [build with adaptive mesh refinement (AMR) suppport]),
              enableamr=$enableval,
              enableamr=yes)

if test "$enableamr" != no ; then
  AC_DEFINE(ENABLE_AMR, 1,
           [Flag indicating if the library should be built with AMR support])
  AC_MSG_RESULT(<<< Configuring library with AMR support >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Variational smoother -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(vsmoother,
              AC_HELP_STRING([--enable-vsmoother],
                             [build with variational smoother suppport]),
              enablevsmoother=$enableval,
              enablevsmoother=yes)

if test "$enablevsmoother" != no ; then
  AC_DEFINE(ENABLE_VSMOOTHER, 1,
           [Flag indicating if the library should be built with variational smoother support])
  AC_MSG_RESULT(<<< Configuring library with variational smoother support >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Periodic BCs -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(periodic,
              AC_HELP_STRING([--enable-periodic],
                             [build with periodic boundary condition suppport]),
              enableperiodic=$enableval,
              enableperiodic=yes)

if test "$enableperiodic" != no ; then
  AC_DEFINE(ENABLE_PERIODIC, 1,
           [Flag indicating if the library should be built with periodic boundary condition support])
  AC_MSG_RESULT(<<< Configuring library with periodic BC support >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# NodeConstraints -- disabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(nodeconstraint,
              AC_HELP_STRING([--enable-nodeconstraint],
                             [build with node constraints suppport]),
              enablenodeconstraint=$enableval,
              enablenodeconstraint=$enableeverything)

if test "$enablenodeconstraint" != no ; then
  AC_DEFINE(ENABLE_NODE_CONSTRAINTS, 1,
           [Flag indicating if the library should be built with node constraints support])
  AC_MSG_RESULT(<<< Configuring library with nodeconstraints support >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Mesh == ParallelMesh -- disabled until it's debugged/finished
# -------------------------------------------------------------
AC_ARG_ENABLE(parmesh,
              AC_HELP_STRING([--enable-parmesh],
                             [Use experimental ParallelMesh as Mesh]),
              enableparmesh=$enableval,
              enableparmesh=no)

if test "$enableparmesh" != no ; then
  AC_DEFINE(ENABLE_PARMESH, 1,
	   [Flag indicating if the library should use the experimental
ParallelMesh as its default Mesh type])
  AC_MSG_RESULT(<<< Configuring library to use ParallelMesh >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Ghosted instead of Serial local vectors -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(ghosted,
              AC_HELP_STRING([--enable-ghosted],
                             [Use ghosted local vectors when available]),
              enableghosted=$enableval,
              enableghosted=yes)

if test "$enableghosted" != no ; then
  AC_DEFINE(ENABLE_GHOSTED, 1,
	   [Flag indicating if the library should use ghosted local vectors])
  AC_MSG_RESULT(<<< Configuring library to use ghosted local vectors >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# 1D or 1D/2D only -- disabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(1D-only,
              AC_HELP_STRING([--enable-1D-only],
                             [build with support for 1D meshes only]),
              enable1D=$enableval,
              enable1D=no)

AC_ARG_ENABLE(2D-only,
              AC_HELP_STRING([--enable-2D-only],
                             [build with support for 1D and 2D meshes only]),
              enable2D=$enableval,
              enable2D=no)

if test "$enable2D" != no ; then
  AC_DEFINE(DIM, 2,
           [Integer indicating the highest spatial dimensionality supported by libMesh])
  AC_MSG_RESULT(<<< Configuring library for 1D/2D meshes only >>>)
elif test "$enable1D" != no ; then
  AC_DEFINE(DIM, 1,
           [Integer indicating the highest spatial dimensionality supported by libMesh])
  AC_MSG_RESULT(<<< Configuring library for 1D meshes only >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# higher order shapes -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(pfem,
              AC_HELP_STRING([--enable-pfem],
                             [build with support for higher order p-FEM shapes]),
              enablepfem=$enableval,
              enablepfem=yes)

if test "$enablepfem" != no ; then
  AC_DEFINE(ENABLE_HIGHER_ORDER_SHAPES, 1,
           [Flag indicating if the library should offer higher order p-FEM shapes])
  AC_MSG_RESULT(<<< Configuring library with higher order p-FEM shapes >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Infinite Elements  -- disabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(ifem,
              AC_HELP_STRING([--enable-ifem],
                             [build with infinite elements]),
              enableifem=$enableval,
              enableifem=$enableeverything)

if test "$enableifem" != no ; then
  AC_DEFINE(ENABLE_INFINITE_ELEMENTS, 1,
           [Flag indicating if the library should be built with infinite elements])
  AC_MSG_RESULT(<<< Configuring library with infinite elements >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Second Derivative Calculations -- disabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(second,
              AC_HELP_STRING([--enable-second],
                             [build with second derivatives]),
              enablesecond=$enableval,
              enablesecond=yes)

if test "$enablesecond" != no ; then
  AC_DEFINE(ENABLE_SECOND_DERIVATIVES, 1,
           [Flag indicating if the library should be built with second derivatives])
  AC_MSG_RESULT(<<< Configuring library with second derivatives >>>)
fi
# -------------------------------------------------------------



# --------------------------------------------------------------
# XDR binary IO support - enabled by default
# --------------------------------------------------------------
AC_ARG_ENABLE(xdr,
              AC_HELP_STRING([--enable-xdr],
                             [enable XDR platform-independent binary I/O]),
              enablexdr=$enableval,
              enablexdr=yes)

if test "$enablexdr" != no ; then
   AC_CHECK_HEADERS(rpc/rpc.h,
                    [
                     AC_CHECK_FUNC(xdrstdio_create,
                                   [
                                     AC_DEFINE(HAVE_XDR, 1,
                                               [Flag indicating headers and libraries for XDR IO are available])
                                     echo "<<< Configuring library with XDR support >>>"
                                   ],
                                   [enablexdr=no])
                    ],
                    [enablexdr=no])
fi
AC_SUBST(enablexdr)	
# -------------------------------------------------------------



# -------------------------------------------------------------
# complex numbers -- disabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(complex,
              AC_HELP_STRING([--enable-complex],
                             [build with complex number support]),
              enablecomplex=$enableval,
              enablecomplex=no)

if test "$enablecomplex" != no ; then
  AC_DEFINE(USE_COMPLEX_NUMBERS, 1,
     [Flag indicating if the library should be built using complex numbers])
  AC_MSG_RESULT(<<< Configuring library with complex number support >>>)
  AC_SUBST(enablecomplex)

else
  AC_DEFINE(USE_REAL_NUMBERS, 1,
     [Flag indicating if the library should be built using real numbers])
  AC_MSG_RESULT(<<< Configuring library with real number support >>>)
  AC_SUBST(enablecomplex)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# Reference Counting -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(reference-counting,
              AC_HELP_STRING([--enable-reference-counting],
                             [build with reference counting suppport]),
              enablerefct=$enableval,
              enablerefct=yes)

if test "$enablerefct" != no ; then
  if test "$ac_cv_cxx_rtti" = yes; then
    AC_DEFINE(ENABLE_REFERENCE_COUNTING, 1,
             [Flag indicating if the library should be built with reference counting support])
    AC_MSG_RESULT(<<< Configuring library with reference counting support >>>)
  else
    AC_MSG_RESULT(<<< No RTTI: disabling reference counting support >>>)
  fi
fi
# -------------------------------------------------------------


# -------------------------------------------------------------
# Performance Logging -- disabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(perflog,
              AC_HELP_STRING([--enable-perflog],
                             [build with performance logging turned on]),
              enableperflog=$enableval,
              enableperflog=$enableeverything)

if test "$enableperflog" != no ; then
  AC_DEFINE(ENABLE_PERFORMANCE_LOGGING, 1,
           [Flag indicating if the library should be built with performance logging support])
  AC_MSG_RESULT(<<< Configuring library with performance logging support >>>)
fi
# ------------------------------------------------------------



AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----- Done configuring compiler features ----)
AC_MSG_RESULT(---------------------------------------------)
])
