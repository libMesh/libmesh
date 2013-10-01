# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([LIBMESH_CORE_FEATURES],
[
AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(----- Configuring core library features -----)
AC_MSG_RESULT(---------------------------------------------)


# --------------------------------------------------------------
# blocked matrix/vector storage - disabled by default.
#   See http://sourceforge.net/mailarchive/forum.php?thread_name=B4613A7D-0033-43C7-A9DF-5A801217A097%40nasa.gov&forum_name=libmesh-devel
# --------------------------------------------------------------
AC_ARG_ENABLE(blocked-storage,
              [AC_HELP_STRING([--enable-blocked-storage],[Support for blocked matrix/vector storage])],
              enableblockedstorage=$enableval,
              enableblockedstorage=no)

if test "$enableblockedstorage" != no ; then
  AC_MSG_RESULT([<<< Configuring library to use blocked storage data structures >>>])
  AC_DEFINE(ENABLE_BLOCKED_STORAGE, 1,
           [Flag indicating if the library should use blocked matrix/vector storage])
fi
# --------------------------------------------------------------


# --------------------------------------------------------------
# default comm_world - enabled by default
# --------------------------------------------------------------
AC_ARG_ENABLE(default-comm-world,
              [AC_HELP_STRING([--enable-default-comm-world],[Provide global libMesh::CommWorld])],
              enabledefaultcommworld=$enableval,
              enabledefaultcommworld=yes)

AC_SUBST(enabledefaultcommworld)
if test "$enabledefaultcommworld" != no ; then
  AC_MSG_RESULT([>>> WARNING: using a legacy option <<<])
  AC_MSG_RESULT([>>> Configuring library to enable a default libMesh::CommWorld <<<])
else
  AC_MSG_RESULT([<<< Configuring library to disable libMesh::CommWorld >>>])
  AC_DEFINE(DISABLE_COMMWORLD, 1,
           [Flag indicating if the library should disable libMesh::CommWorld])
fi
# --------------------------------------------------------------


# --------------------------------------------------------------
# legacy include paths - disabled by default
# --------------------------------------------------------------
AC_ARG_ENABLE(legacy-include-paths,
              [AC_HELP_STRING([--enable-legacy-include-paths],[allow for e.g. #include "header.h" instead of #include "libmesh/header.h"])],
              enablelegacyincludepaths=$enableval,
              enablelegacyincludepaths=no)

AC_SUBST(enablelegacyincludepaths)
if test "$enablelegacyincludepaths" != no ; then
  AC_MSG_RESULT([>>> WARNING: using a legacy option <<<])
  AC_MSG_RESULT([>>> Configuring library to dump old header paths into include args <<<])
else
  AC_MSG_RESULT([<<< Configuring library to require ``include "libmesh/etc.h"'' style >>>])
fi
# --------------------------------------------------------------


# --------------------------------------------------------------
# legacy "using namespace libMesh" - disabled by default
# --------------------------------------------------------------
AC_ARG_ENABLE(legacy-using-namespace,
              [AC_HELP_STRING([--enable-legacy-using-namespace],[add "using namespace libMesh" to libMesh headers])],
              enablelegacyusingnamespace=$enableval,
              enablelegacyusingnamespace=no)

if test "$enablelegacyusingnamespace" != no ; then
  AC_MSG_RESULT([>>> WARNING: using a legacy option <<<])
  AC_MSG_RESULT([>>> Configuring library to dump names into global namespace <<<])
else
  AC_MSG_RESULT([<<< Configuring library to keep names in libMesh namespace >>>])
  AC_DEFINE(REQUIRE_SEPARATE_NAMESPACE, 1,
           [Flag indicating if the library should keep names in libMesh namespace])
fi
# --------------------------------------------------------------



# -------------------------------------------------------------
# size of boundary_id_type -- default 2 bytes
# -------------------------------------------------------------
AC_ARG_WITH([boundary_id_bytes],
	    AC_HELP_STRING([--with-boundary-id-bytes=<1|2|4|8>],
                           [bytes used per boundary side per boundary_id]),
	    [boundary_bytes="$withval"],
	    [boundary_bytes=2])

case "$boundary_bytes" in
    1)
	AC_DEFINE(BOUNDARY_ID_BYTES, 1, [size of boundary_id])
	;;
    2)
	AC_DEFINE(BOUNDARY_ID_BYTES, 2, [size of boundary_id])
	;;
    4)
	AC_DEFINE(BOUNDARY_ID_BYTES, 4, [size of boundary_id])
	;;
    8)
	AC_DEFINE(BOUNDARY_ID_BYTES, 8, [size of boundary_id])
	;;
    *)
	AC_MSG_RESULT([>>> unrecognized boundary_id size: $boundary_bytes - configuring size...2])
	AC_DEFINE(BOUNDARY_ID_BYTES, 2, [size of boundary_id])
	boundary_bytes=2
	;;
esac
AC_MSG_RESULT([configuring size of boundary_id... $boundary_bytes])
# -------------------------------------------------------------



# -------------------------------------------------------------
# size of dof_id_type -- default 4 bytes
# -------------------------------------------------------------
AC_ARG_WITH([dof_id_bytes],
	    AC_HELP_STRING([--with-dof-id-bytes=<1|2|4|8>],
                           [bytes used per dof object id, dof index]),
	    [dof_bytes="$withval"],
	    [dof_bytes=4])

case "$dof_bytes" in
    1)
	AC_DEFINE(DOF_ID_BYTES, 1, [size of dof_id])
	;;
    2)
	AC_DEFINE(DOF_ID_BYTES, 2, [size of dof_id])
	;;
    4)
	AC_DEFINE(DOF_ID_BYTES, 4, [size of dof_id])
	;;
    8)
	AC_DEFINE(DOF_ID_BYTES, 8, [size of dof_id])
	;;
    *)
	AC_MSG_RESULT([>>> unrecognized dof_id size: $dof_bytes - configuring size...4])
	AC_DEFINE(DOF_ID_BYTES, 4, [size of dof_id])
	dof_bytes=4
	;;
esac
AC_MSG_RESULT([configuring size of dof_id... $dof_bytes])
# -------------------------------------------------------------



# -------------------------------------------------------------
# size of processor_id_type -- default 4 bytes
# -------------------------------------------------------------
AC_ARG_WITH([processor_id_bytes],
	    AC_HELP_STRING([--with-processor-id-bytes=<1|2|4|8>],
                           [bytes used for processor id]),
	    [processor_bytes="$withval"],
	    [processor_bytes=2])

case "$processor_bytes" in
    1)
	AC_DEFINE(PROCESSOR_ID_BYTES, 1, [size of processor_id])
	;;
    2)
	AC_DEFINE(PROCESSOR_ID_BYTES, 2, [size of processor_id])
	;;
    4)
	AC_DEFINE(PROCESSOR_ID_BYTES, 4, [size of processor_id])
	;;
    8)
	AC_DEFINE(PROCESSOR_ID_BYTES, 8, [size of processor_id])
	;;
    *)
	AC_MSG_RESULT([>>> unrecognized processor_id size: $processor_bytes - configuring size...2])
	AC_DEFINE(PROCESSOR_ID_BYTES, 2, [size of processor_id])
	processor_bytes=2
	;;
esac
AC_MSG_RESULT([configuring size of processor_id... $processor_bytes])
# -------------------------------------------------------------



# -------------------------------------------------------------
# size of subdomain_id_type -- default 2 bytes
# -------------------------------------------------------------
AC_ARG_WITH([subdomain_id_bytes],
	    AC_HELP_STRING([--with-subdomain-id-bytes=<1|2|4|8>],
                           [bytes of storage per element used to store the subdomain_id]),
	    [subdomain_bytes="$withval"],
	    [subdomain_bytes=2])

case "$subdomain_bytes" in
    1)
	AC_DEFINE(SUBDOMAIN_ID_BYTES, 1, [size of subdomain_id])
	;;
    2)
	AC_DEFINE(SUBDOMAIN_ID_BYTES, 2, [size of subdomain_id])
	;;
    4)
	AC_DEFINE(SUBDOMAIN_ID_BYTES, 4, [size of subdomain_id])
	;;
    8)
	AC_DEFINE(SUBDOMAIN_ID_BYTES, 8, [size of subdomain_id])
	;;
    *)
	AC_MSG_RESULT([>>> unrecognized subdomain_id size: $subdomain_bytes - configuring size...2])
	AC_DEFINE(SUBDOMAIN_ID_BYTES, 2, [size of subdomain_id])
	subdomain_bytes=2
	;;
esac
AC_MSG_RESULT([configuring size of subdomain_id... $subdomain_bytes])
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



# -------------------------------------------------------------
# unique_id -- disable by default
# -------------------------------------------------------------
AC_ARG_ENABLE(unique-id,
              AC_HELP_STRING([--enable-unique-id],
                             [build with unique id suppport]),
	      [case "${enableval}" in
	          yes)  enableuniqueid=yes ;;
		   no)	enableuniqueid=no ;;
		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-unique-id) ;;
	       esac],
	       [enableuniqueid=$enableeverything])

if test "$enableuniqueid" = yes ; then
  AC_DEFINE(ENABLE_UNIQUE_ID, 1,
           [Flag indicating if the library should be built with unique id support])
  AC_MSG_RESULT(<<< Configuring library with unique id support >>>)
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# size of unique_id_type -- default 8 bytes
# -------------------------------------------------------------
AC_ARG_WITH([unique_id_bytes],
	    AC_HELP_STRING([--with-unique-id-bytes=<1|2|4|8>],
                           [bytes used per unique id]),
	    [unique_bytes="$withval"],
	    [unique_bytes=8])

case "$unique_bytes" in
    1)
	AC_DEFINE(UNIQUE_ID_BYTES, 1, [size of unique_id])
	;;
    2)
	AC_DEFINE(UNIQUE_ID_BYTES, 2, [size of unique_id])
	;;
    4)
	AC_DEFINE(UNIQUE_ID_BYTES, 4, [size of unique_id])
	;;
    8)
	AC_DEFINE(UNIQUE_ID_BYTES, 8, [size of unique_id])
	;;
    *)
	AC_MSG_RESULT([>>> unrecognized unique_id size: $unique_bytes - configuring size...8])
	AC_DEFINE(UNIQUE_ID_BYTES, 8, [size of unique_id])
	unique_bytes=8
	;;
esac

if test "$enableuniqueid" = yes ; then
   AC_MSG_RESULT([configuring size of unique_id... $unique_bytes])
fi
# -------------------------------------------------------------



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
# Dirichlet BC constraints -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(dirichlet,
              AC_HELP_STRING([--enable-dirichlet],
                             [build with Dirichlet boundary constraint support]),
              enabledirichlet=$enableval,
              enabledirichlet=yes)

if test "$enabledirichlet" != no ; then
  AC_DEFINE(ENABLE_DIRICHLET, 1,
           [Flag indicating if the library should be built with Dirichlet boundary constraint support])
  AC_MSG_RESULT(<<< Configuring library with Dirichlet constraint support >>>)
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
  AC_MSG_RESULT(<<< Configuring library with node constraints support >>>)
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

AM_CONDITIONAL(LIBMESH_ENABLE_INFINITE_ELEMENTS, test x$enableifem != no )

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
                    [
                      AC_CHECK_HEADERS(rpc/xdr.h,
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
                     ])
fi
# -------------------------------------------------------------



# -------------------------------------------------------------
# complex numbers -- disabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(complex,
              AC_HELP_STRING([--enable-complex],
                             [build with complex number support]),
 	      [case "${enableval}" in
	          yes)  enablecomplex=yes ;;
		   no)  enablecomplex=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-complex) ;;
	       esac],
	       [enablecomplex=no])

if test "$enablecomplex" != no ; then
  AC_DEFINE(USE_COMPLEX_NUMBERS, 1,
     [Flag indicating if the library should be built using complxex numbers])
  AC_MSG_RESULT(<<< Configuring library with complex number support >>>)

else
  AC_DEFINE(USE_REAL_NUMBERS, 1,
     [Flag indicating if the library should be built using real numbers])
  AC_MSG_RESULT(<<< Configuring library with real number support >>>)
fi

AM_CONDITIONAL(LIBMESH_ENABLE_COMPLEX, test x$enablecomplex = xyes)
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



# -------------------------------------------------------------
# Examples - enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(examples,
              AC_HELP_STRING([--enable-examples],
                             [support compilation, installation, & running example suite]),
 	      [case "${enableval}" in
	          yes)  enableexamples=yes ;;
		   no)  enableexamples=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-examples) ;;
	       esac],
              [enableexamples=yes])

if test "$enableexamples" = yes ; then
  AC_MSG_RESULT(<<< Configuring library example suite support >>>)
fi
AM_CONDITIONAL(LIBMESH_ENABLE_EXAMPLES, test x$enableexamples = xyes)
# ------------------------------------------------------------

AC_MSG_RESULT(---------------------------------------------)
AC_MSG_RESULT(-- Done configuring core library features ---)
AC_MSG_RESULT(---------------------------------------------)
])
