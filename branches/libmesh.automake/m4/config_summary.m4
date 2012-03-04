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
echo Package version............... : $PACKAGE-$VERSION
echo
echo C++ compiler type............. : $GXX_VERSION
echo C++ compiler.................. : $CXX
echo C compiler.................... : $CC
echo Fortran compiler.............. : $FC
echo Build Method.................. : $METHOD
echo CPPFLAGS...................... : $CPPFLAGS_METHOD
echo CXXFLAGS...................... : $CXXFLAGS_METHOD
echo CFLAGS........................ : $CFLAGS_METHOD
echo Install dir................... : $prefix 
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo SVN revision number........... : $BUILD_VERSION

######################################################################################
# echo
# echo Solver Configuration:

#  if test "x$ENABLE_COMP_NS" = "x1"; then
#    echo '   'Compressible Navier-Stokes.......... : yes
#  else
#    echo '   'Compressible Navier-Stokes.......... : no
#  fi

# # if test "x$ENABLE_PRIM_NS" = "x1"; then
# #   echo '   'Primitive Navier-Stokes.... : yes
# # else
# #   echo '   'Primitive Navier-Stokes.... : no
# # fi

#  if test "x$ENABLE_HEAT_TRANSFER" = "x1"; then
#    echo '   'Heat transfer support............... : yes
#  else
#    echo '   'Heat transfer support............... : no
#  fi

#  if test "x$ENABLE_FASTMATH" = "x1"; then
#    echo '   'Fast math function approximations... : yes
#  else
#    echo '   'Fast math function approximations... : no
#  fi


######################################################################################
echo
echo Library Features:
echo '  'adaptive mesh refinement.... : $enableamr
echo '  'complex variables........... : $enablecomplex
echo '  'ghosted vectors............. : $enableghosted
echo '  'high-order shape functions.. : $enablepfem
echo '  'infinite elements........... : $enableifem
echo '  'node constraints............ : $enablenodeconstraint
echo '  'parallel mesh............... : $enableparmesh
echo '  'performance logging......... : $enableperflog
echo '  'periodic boundary conditions : $enableperiodic
echo '  'reference counting.......... : $enablerefct
echo '  'shape function 2nd derivs... : $enablesecond
echo '  'stack trace files........... : $enabletracefiles
echo '  'subdomain id size........... : $subdomain_bytes bytes
echo '  'variational smoother........ : $enablevsmoother
echo '  'xdr binary I/O.............. : $enablexdr


		   
######################################################################################
if test "x$enableoptional" = "xyes"; then
  echo
  echo Optional Packages:
  echo '  'eigen....................... : $enableeigen
  echo '  'exodus...................... : $enableexodus
  echo '  'fparser..................... : $enablefparser
  echo '  'glpk........................ : $enableglpk
  echo '  'gmv......................... : $enablegmv
  echo '  'gzstream.................... : $enablegz
  echo '  'laspack..................... : $enablelaspack
  echo '  'libhilbert.................. : $enablelibhilbert
  echo '  'metis....................... : $enablemetis
  echo '  'mpi......................... : $enablempi
  echo '  'nemesis..................... : $enablenemesis
  echo '  'netcdf...................... : $enablenetcdf
  echo '  'openmp...................... : $enableopenmp
  echo '  'parmetis.................... : $enableparmetis
  echo '  'petsc....................... : $enablepetsc
  echo '  'sfcurves.................... : $enablesfc
  echo '  'slepc....................... : $enableslepc
  echo '  'tbb......................... : $enabletbb
  echo '  'tetgen...................... : $enabletetgen
  echo '  'triangle.................... : $enabletriangle
  echo '  'trilinos.................... : $enabletrilinos
  echo '  'vtk......................... : $enablevtk
  echo
  echo '  'libmesh_optional_INCLUDES... : $libmesh_optional_INCLUDES
  echo '  'libmesh_optional_LIBS....... : $libmesh_optional_LIBS
fi
		   
echo
echo '-------------------------------------------------------------------------------'

echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
