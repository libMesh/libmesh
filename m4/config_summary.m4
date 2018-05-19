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
#   git log -n1 m4/config_summary.m4
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

######################################################################################
AS_ECHO([])
AS_ECHO(["----------------------------------- SUMMARY -----------------------------------"])
AS_ECHO([])
AS_ECHO(["Package version.................... : $PACKAGE-$VERSION"])
AS_ECHO([])
AS_ECHO(["C++ compiler type.................. : $GXX_VERSION"])
AS_ECHO(["C++ compiler....................... : $CXX"])
AS_ECHO(["C compiler......................... : $CC"])
AS_ECHO(["Fortran compiler................... : $FC"])
AS_ECHO(["Build Methods...................... : $METHODS"])
AS_ECHO([])
for method in ${METHODS}; do
     AS_CASE("${method}",
             [opt],   [AS_ECHO(["CPPFLAGS...(opt)................... : $CPPFLAGS_OPT $CPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(opt)................... : $CXXFLAGS_OPT $CXXFLAGS"])
                       AS_ECHO(["CFLAGS.....(opt)................... : $CFLAGS_OPT $CFLAGS"])],
             [devel], [AS_ECHO(["CPPFLAGS...(devel)................. : $CPPFLAGS_DEVEL $CPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(devel)................. : $CXXFLAGS_DEVEL $CXXFLAGS"])
                       AS_ECHO(["CFLAGS.....(devel)................. : $CFLAGS_DEVEL $CFLAGS"])],
             [dbg],   [AS_ECHO(["CPPFLAGS...(dbg)................... : $CPPFLAGS_DBG $CPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(dbg)................... : $CXXFLAGS_DBG $CXXFLAGS"])
                       AS_ECHO(["CFLAGS.....(dbg)................... : $CFLAGS_DBG $CFLAGS"])],
             [prof],  [AS_ECHO(["CPPFLAGS...(prof).................. : $CPPFLAGS_PROF $CPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(prof).................. : $CXXFLAGS_PROF $CXXFLAGS"])
                       AS_ECHO(["CFLAGS.....(prof).................. : $CFLAGS_PROF $CFLAGS"])],
             [oprof], [AS_ECHO(["CPPFLAGS...(oprof)................. : $CPPFLAGS_OPROF $CPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(oprof)................. : $CXXFLAGS_OPROF $CXXFLAGS"])
                       AS_ECHO(["CFLAGS.....(oprof)................. : $CFLAGS_OPROF $CFLAGS"])])

     dnl blank line
     AS_ECHO([])
done
AS_ECHO(["Install dir........................ : $prefix"])
AS_ECHO(["Build user......................... : $USER"])
AS_ECHO(["Build host......................... : $BUILD_HOST"])
AS_ECHO(["Build architecture................. : $BUILD_ARCH"])
AS_ECHO(["Git revision....................... : $BUILD_VERSION"])

dnl Print out library features
AS_ECHO([])
AS_ECHO(["Library Features:"])
AS_ECHO(["  library warnings................. : $enablewarnings"])
AS_ECHO(["  library deprecated code support.. : $enabledeprecated"])
AS_ECHO(["  adaptive mesh refinement......... : $enableamr"])
AS_ECHO(["  blocked matrix/vector storage.... : $enableblockedstorage"])
AS_ECHO(["  complex variables................ : $enablecomplex"])
AS_ECHO(["  example suite.................... : $enableexamples"])
AS_ECHO(["  ghosted vectors.................. : $enableghosted"])
AS_ECHO(["  high-order shape functions....... : $enablepfem"])
AS_ECHO(["  unique-id support................ : $enableuniqueid"])
AS_ECHO(["  id size (boundaries)............. : $boundary_bytes bytes"])
AS_ECHO(["  id size (dofs)................... : $dof_bytes bytes"])
AS_IF([test "x$enableuniqueid" = "xyes"],
      [AS_ECHO(["  id size (unique)................. : $unique_bytes bytes"])])
AS_ECHO(["  id size (processors)............. : $processor_bytes bytes"])
AS_ECHO(["  id size (subdomains)............. : $subdomain_bytes bytes"])
AS_ECHO(["  infinite elements................ : $enableifem"])
AS_ECHO(["  Dirichlet constraints............ : $enabledirichlet"])
AS_ECHO(["  node constraints................. : $enablenodeconstraint"])
AS_ECHO(["  parallel mesh.................... : $enableparmesh"])
AS_ECHO(["  performance logging.............. : $enableperflog"])
AS_ECHO(["  periodic boundary conditions..... : $enableperiodic"])
AS_ECHO(["  reference counting............... : $enablerefct"])
AS_ECHO(["  shape function 2nd derivatives... : $enablesecond"])
AS_ECHO(["  stack trace files................ : $enabletracefiles"])
AS_ECHO(["  track node valence............... : $enablenodevalence"])
AS_ECHO(["  variational smoother............. : $enablevsmoother"])
AS_ECHO(["  xdr binary I/O................... : $enablexdr"])
AS_IF([test "x$enablelegacyincludepaths" = "xyes"],
      [AS_ECHO(["  non-prefixed include paths....... : $enablelegacyincludepaths ***LEGACY FEATURE***"])])
AS_IF([test "x$enablelegacyusingnamespace" = "xyes"],
      [AS_ECHO(["  adding using namespace libMesh... : $enablelegacyusingnamespace ***LEGACY FEATURE***"])])
AS_IF([test "x$enabledefaultcommworld" = "xyes"],
      [AS_ECHO(["  providing libMesh::CommWorld..... : $enabledefaultcommworld ***LEGACY FEATURE***"])])


dnl Print out which optional libraries have been configured.
AS_IF([test "x$enableoptional" = "xyes"],
      [
        AS_ECHO([])
        AS_ECHO(["Optional Packages:"])
        AS_ECHO(["  boost............................ : $enableboost"])
        AS_ECHO(["  capnproto........................ : $enablecapnproto"])
        AS_ECHO(["  cppunit.......................... : $enablecppunit"])
        AS_ECHO(["  curl............................. : $enablecurl"])
        AS_ECHO(["  eigen............................ : $enableeigen"])
        AS_ECHO(["  exodus........................... : $enableexodus"])
        AS_IF([test "x$exodusversion" != "xno"],
              [AS_ECHO(["     version....................... : $exodusversion"])])
        AS_ECHO(["  fparser.......................... : $enablefparser"])
        AS_IF([test "x$enablefparser" = "xyes" && test "x$enablefparserdevel" = "xno"],
              [AS_ECHO(["     build from version............ : release"])])
        AS_IF([test "x$enablefparser" = "xyes" && test "x$enablefparserdevel" = "xyes"],
              [AS_ECHO(["     build from version............ : devel"])])
        AS_IF([test "x$enablefparser" = "xyes" && test "x$enablefparserdebugging" = "xyes"],
              [AS_ECHO(["     fparser debugging............. : enabled"])])
        AS_ECHO(["  glpk............................. : $enableglpk"])
        AS_ECHO(["  gmv.............................. : $enablegmv"])
        AS_ECHO(["  gzstream......................... : $enablegz"])
        AS_ECHO(["  hdf5............................. : $enablehdf5"])
        AS_ECHO(["  laspack.......................... : $enablelaspack"])
        AS_ECHO(["  libhilbert....................... : $enablelibhilbert"])
        AS_ECHO(["  metis............................ : $enablemetis"])
        AS_ECHO(["  mpi.............................. : $enablempi"])
        AS_ECHO(["  nanoflann........................ : $enablenanoflann"])
        AS_ECHO(["  nemesis.......................... : $enablenemesis"])
        AS_IF([test "x$nemesisversion" != "xno"],
              [AS_ECHO(["     version....................... : $nemesisversion"])])
        AS_ECHO(["  netcdf........................... : $enablenetcdf"])
        AS_IF([test "x$netcdfversion" != "xno"],
              [AS_ECHO(["     version....................... : $netcdfversion"])])
        AS_ECHO(["  nlopt............................ : $enablenlopt"])
        AS_ECHO(["  parmetis......................... : $enableparmetis"])
        AS_ECHO(["  petsc............................ : $enablepetsc"])
        AS_IF([test "x$enablepetsc" = "xyes"],
              [AS_ECHO(["     version....................... : $petscversion"])])
        AS_ECHO(["  qhull............................ : $enableqhull"])
        AS_ECHO(["  sfcurves......................... : $enablesfc"])
        AS_ECHO(["  slepc............................ : $enableslepc"])
        AS_IF([test "x$enableslepc" = "xyes"],
              [AS_ECHO(["     version....................... : $slepcversion"])])
        AS_ECHO(["  thread model..................... : $found_thread_model"])
        AS_ECHO(["  c++ rtti ........................ : $ac_cv_cxx_rtti"])
        AS_ECHO(["  tecio............................ : $enabletecio"])
        AS_ECHO(["  tecplot...(vendor binaries)...... : $enabletecplot"])
        AS_ECHO(["  tetgen........................... : $enabletetgen"])
        AS_ECHO(["  triangle......................... : $enabletriangle"])
        AS_ECHO(["  trilinos......................... : $enabletrilinos"])
        AS_IF([test "x$enabletrilinos" = "xyes"],
              [
                AS_ECHO(["     AztecOO....................... : $enableaztecoo"])
                AS_ECHO(["     NOX........................... : $enablenox"])
                AS_ECHO(["     ML............................ : $enableml"])
                AS_ECHO(["     Tpetra........................ : $enabletpetra"])
                AS_ECHO(["     DTK........................... : $enabledtk"])
                AS_ECHO(["     Ifpack........................ : $enableifpack"])
                AS_ECHO(["     Epetra........................ : $enableepetra"])
                AS_ECHO(["     EpetraExt..................... : $enableepetraext"])
              ])
        AS_ECHO(["  vtk.............................. : $enablevtk"])
        AS_IF([test "x$enablevtk" = "xyes"],
              [AS_ECHO(["     version....................... : $vtkversion"])])
        AS_ECHO([])
        AS_IF([test "x$libmesh_optional_INCLUDES" != "x"],
              [
                AS_ECHO(["  libmesh_optional_INCLUDES........ : $libmesh_optional_INCLUDES"])
                AS_ECHO([])
              ])
        AS_IF([test "x$libmesh_optional_LIBS" != "x"],
              [
                AS_ECHO(["  libmesh_optional_LIBS............ : $libmesh_optional_LIBS"])
                AS_ECHO([])
              ])
      ])
AS_ECHO(["-------------------------------------------------------------------------------"])
AS_ECHO(["Configure complete, now type 'make' and then 'make install'."])
AS_ECHO([])

])
