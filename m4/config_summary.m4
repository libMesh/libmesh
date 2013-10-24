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
echo CPPFLAGS...\(opt\)................... : $CPPFLAGS_OPT $CPFLAGS
echo CXXFLAGS...\(opt\)................... : $CXXFLAGS_OPT $CXXFLAGS
echo CFLAGS.....\(opt\)................... : $CFLAGS_OPT $CFLAGS
     ;;
         devel)
echo CPPFLAGS...\(devel\)................. : $CPPFLAGS_DEVEL $CPFLAGS
echo CXXFLAGS...\(devel\)................. : $CXXFLAGS_DEVEL $CXXFLAGS
echo CFLAGS.....\(devel\)................. : $CFLAGS_DEVEL $CFLAGS
     ;;
         dbg)
echo CPPFLAGS...\(dbg\)................... : $CPPFLAGS_DBG $CPFLAGS
echo CXXFLAGS...\(dbg\)................... : $CXXFLAGS_DBG $CXXFLAGS
echo CFLAGS.....\(dbg\)................... : $CFLAGS_DBG $CFLAGS
     ;;
         prof)
echo CPPFLAGS...\(prof\).................. : $CPPFLAGS_PROF $CPFLAGS
echo CXXFLAGS...\(prof\).................. : $CXXFLAGS_PROF $CXXFLAGS
echo CFLAGS.....\(prof\).................. : $CFLAGS_PROF $CFLAGS
     ;;
         oprof)
echo CPPFLAGS...\(oprof\)................. : $CPPFLAGS_OPROF $CPFLAGS
echo CXXFLAGS...\(oprof\)................. : $CXXFLAGS_OPROF $CXXFLAGS
echo CFLAGS.....\(oprof\)................. : $CFLAGS_OPROF $CFLAGS
     esac
     echo " "
done
echo Install dir........................ : $prefix
echo Build user......................... : $USER
echo Build host......................... : $BUILD_HOST
echo Build architecture................. : $BUILD_ARCH
echo Git revision....................... : $BUILD_VERSION

######################################################################################
echo
echo Library Features:
echo '  library warnings................. :' $enablewarnings
echo '  adaptive mesh refinement......... :' $enableamr
echo '  blocked matrix/vector storage.... :' $enableblockedstorage
echo '  complex variables................ :' $enablecomplex
echo '  example suite.................... :' $enableexamples
echo '  ghosted vectors.................. :' $enableghosted
echo '  high-order shape functions....... :' $enablepfem
echo '  unique-id support................ :' $enableuniqueid
echo '  id size (boundaries)............. :' $boundary_bytes bytes
echo '  id size (dofs)................... :' $dof_bytes bytes
if (test "x$enableuniqueid" = "xyes"); then
echo '  id size (unique)................. :' $unique_bytes bytes
fi
echo '  id size (processors)............. :' $processor_bytes bytes
echo '  id size (subdomains)............. :' $subdomain_bytes bytes
echo '  infinite elements................ :' $enableifem
echo '  Dirichlet constraints............ :' $enabledirichlet
echo '  node constraints................. :' $enablenodeconstraint
echo '  parallel mesh.................... :' $enableparmesh
echo '  performance logging.............. :' $enableperflog
echo '  periodic boundary conditions..... :' $enableperiodic
echo '  reference counting............... :' $enablerefct
echo '  shape function 2nd derivatives... :' $enablesecond
echo '  stack trace files................ :' $enabletracefiles
echo '  variational smoother............. :' $enablevsmoother
echo '  xdr binary I/O................... :' $enablexdr
if (test "x$enablelegacyincludepaths" = "xyes"); then
echo '  non-prefixed include paths....... :' $enablelegacyincludepaths ***LEGACY FEATURE***
fi
if (test "x$enablelegacyusingnamespace" = "xyes"); then
echo '  adding using namespace libMesh... :' $enablelegacyusingnamespace ***LEGACY FEATURE***
fi
if (test "x$enabledefaultcommworld" = "xyes"); then
echo '  providing libMesh::CommWorld..... :' $enabledefaultcommworld ***LEGACY FEATURE***
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
  if (test "x$exodusversion" != "xno"); then
  echo '     'version....................... : $exodusversion
  fi
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
  echo '  'hdf5............................. : $enablehdf5
  echo '  'laspack.......................... : $enablelaspack
  echo '  'libhilbert....................... : $enablelibhilbert
  echo '  'metis............................ : $enablemetis
  echo '  'mpi.............................. : $enablempi
  echo '  'nanoflann........................ : $enablenanoflann
  echo '  'nemesis.......................... : $enablenemesis
  if (test "x$nemesisversion" != "xno"); then
  echo '     'version....................... : $nemesisversion
  fi
  echo '  'netcdf........................... : $enablenetcdf
  if (test "x$netcdfversion" != "xno"); then
  echo '     'version....................... : $netcdfversion
  fi
  echo '  'openmp........................... : $enableopenmp
  echo '  'parmetis......................... : $enableparmetis
  echo '  'petsc............................ : $enablepetsc
  if (test "x$enablepetsc" = "xyes"); then
  echo '     'version....................... : $petscversion
  fi
  echo '  'pthreads......................... : $enablepthreads
  echo '  'sfcurves......................... : $enablesfc
  echo '  'slepc............................ : $enableslepc
  echo '  'tbb.............................. : $enabletbb
  echo '  'c++ threads...................... : $enablecppthreads
  if (test "x$enablecppthreads" = "xyes"); then
  echo '     'flavor........................ : $cppthreadflavor
  fi
  echo '  'c++ rtti ........................ : $ac_cv_cxx_rtti
  echo '  'tecio............................ : $enabletecio
  echo '  'tecplot...\(vendor binaries\)...... : $enabletecplot
  echo '  'tetgen........................... : $enabletetgen
  echo '  'triangle......................... : $enabletriangle
  echo '  'trilinos......................... : $enabletrilinos
  if (test "x$enabletrilinos" = "xyes"); then
  echo '     'AztecOO....................... : $enableaztecoo
  echo '     'NOX........................... : $enablenox
  echo '     'ML............................ : $enableml
  echo '     'Tpetra........................ : $enabletpetra
  echo '     'DTK........................... : $enabledtk
  fi
  echo '  'vtk.............................. : $enablevtk
  if (test "x$enablevtk" = "xyes"); then
  echo '     'version....................... : $vtkversion
  fi
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
