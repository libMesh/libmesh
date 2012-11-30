# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#
# LAST MODIFICATION
#
#   $Id: config_summary.m4 27667 2012-02-05 03:36:06Z benkirk $
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

######################################################################################
echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version.................... : $PACKAGE-$VERSION
echo
echo C++ compiler type.................. : $GXX_VERSION
echo C++ compiler....................... : $CXX
echo C compiler......................... : $CC
echo Fortran compiler................... : $FC
echo Build Methods...................... : $METHODS
echo " "
for method in ${METHODS}; do
     case "${method}" in
         opt)   
echo CPPFLAGS...\(opt\)................... : $CPPFLAGS_OPT
echo CXXFLAGS...\(opt\)................... : $CXXFLAGS_OPT
echo CFLAGS.....\(opt\)................... : $CFLAGS_OPT
     ;;					 
         devel)   
echo CPPFLAGS...\(devel\)................. : $CPPFLAGS_DEVEL
echo CXXFLAGS...\(devel\)................. : $CXXFLAGS_DEVEL
echo CFLAGS.....\(devel\)................. : $CFLAGS_DEVEL
     ;;					 
         dbg)   
echo CPPFLAGS...\(dbg\)................... : $CPPFLAGS_DBG
echo CXXFLAGS...\(dbg\)................... : $CXXFLAGS_DBG
echo CFLAGS.....\(dbg\)................... : $CFLAGS_DBG
     ;;					 
         prof)   
echo CPPFLAGS...\(prof\).................. : $CPPFLAGS_PROF
echo CXXFLAGS...\(prof\).................. : $CXXFLAGS_PROF
echo CFLAGS.....\(prof\).................. : $CFLAGS_PROF
     ;;					 
         oprof)   
echo CPPFLAGS...\(oprof\)................. : $CPPFLAGS_OPROF
echo CXXFLAGS...\(oprof\)................. : $CXXFLAGS_OPROF
echo CFLAGS.....\(oprof\)................. : $CFLAGS_OPROF
     esac
     echo " "
done
echo Install dir........................ : $prefix 
echo Build user......................... : $USER
echo Build host......................... : $BUILD_HOST
#echo Configure date..................... : $BUILD_DATE
echo Build architecture................. : $BUILD_ARCH
echo SVN revision number................ : $BUILD_VERSION

######################################################################################
echo
echo Library Features:
echo '  'adaptive mesh refinement......... : $enableamr
echo '  'complex variables................ : $enablecomplex
echo '  'example suite.................... : $enableexamples
echo '  'ghosted vectors.................. : $enableghosted
echo '  'high-order shape functions....... : $enablepfem
echo '  'infinite elements................ : $enableifem
echo '  'Dirichlet constraints............ : $enabledirichlet
echo '  'node constraints................. : $enablenodeconstraint
echo '  'parallel mesh.................... : $enableparmesh
echo '  'performance logging.............. : $enableperflog
echo '  'periodic boundary conditions..... : $enableperiodic
echo '  'reference counting............... : $enablerefct
echo '  'shape function 2nd derivatives... : $enablesecond
echo '  'stack trace files................ : $enabletracefiles
echo '  'subdomain id size................ : $subdomain_bytes bytes
echo '  'variational smoother............. : $enablevsmoother
echo '  'xdr binary I/O................... : $enablexdr
if (test "x$enablelegacyincludepaths" = "xyes"); then
echo '  'non-prefixed include paths....... : $enablelegacyincludepaths ***LEGACY FEATURE***
fi
if (test "x$enablelegacyusingnamespace" = "xyes"); then
echo '  'adding using namespace libMesh... : $enablelegacyusingnamespace ***LEGACY FEATURE***
fi


		   
######################################################################################
if (test "x$enableoptional" = "xyes"); then
  echo
  echo Optional Packages:
  echo '  'boost............................ : $enableboost
  echo '  'cppunit.......................... : $enablecppunit
#  if (test "x$enablecppunit" = "xyes"); then
#  echo '     'CPPUNIT_CFLAGS................ : $CPPUNIT_CFLAGS
#  echo '     'CPPUNIT_LIBS.................. : $CPPUNIT_LIBS
#  fi
  echo '  'eigen............................ : $enableeigen
  echo '  'exodus........................... : $enableexodus
  echo '  'fparser.......................... : $enablefparser
  if (test "x$enablefparser" = "xyes" -a "x$enablefparserdevel" = "xno"); then
  echo '     'build from version............ : release
  fi
  if (test "x$enablefparser" = "xyes" -a "x$enablefparserdevel" = "xyes"); then
  echo '     'build from version............ : devel
  fi
  echo '  'glpk............................. : $enableglpk
  echo '  'gmv.............................. : $enablegmv
  echo '  'gzstream......................... : $enablegz
  echo '  'laspack.......................... : $enablelaspack
  echo '  'libhilbert....................... : $enablelibhilbert
  echo '  'metis............................ : $enablemetis
  echo '  'mpi.............................. : $enablempi
  echo '  'nemesis.......................... : $enablenemesis
  echo '  'netcdf........................... : $enablenetcdf
  echo '  'openmp........................... : $enableopenmp
  echo '  'parmetis......................... : $enableparmetis
  echo '  'petsc............................ : $enablepetsc
  if (test "x$enablepetsc" = "xyes"); then
  echo '     'version....................... : $petscversion
  fi
  echo '  'sfcurves......................... : $enablesfc
  echo '  'slepc............................ : $enableslepc
  echo '  'tbb.............................. : $enabletbb
  echo '  'c++ threads...................... : $enablecppthreads
  if (test "x$enablecppthreads" = "xyes"); then
  echo '     'flavor........................ : $cppthreadflavor
  fi
  echo '  'tecio............................ : $enabletecio
  echo '  'tecplot...\(vendor binaries\)...... : $enabletecplot
  echo '  'tetgen........................... : $enabletetgen
  echo '  'triangle......................... : $enabletriangle
  echo '  'trilinos......................... : $enabletrilinos
  if (test "x$enabletrilinos" = "xyes"); then
  echo '     'AztecOO....................... : $enableaztecoo
  echo '     'NOX........................... : $enablenox
  echo '     'ML............................ : $enableml
  fi
  echo '  'vtk.............................. : $enablevtk
  echo
  if (test "x$libmesh_optional_INCLUDES" != "x"); then
  echo '  'libmesh_optional_INCLUDES........ : $libmesh_optional_INCLUDES
  echo
  fi
  if (test "x$libmesh_optional_LIBS" != "x"); then
  echo '  'libmesh_optional_LIBS............ : $libmesh_optional_LIBS
  echo
  fi	
fi		   
echo '-------------------------------------------------------------------------------'

echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
