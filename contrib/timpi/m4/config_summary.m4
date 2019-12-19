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
#   2019-10-31
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

AS_ECHO([])
AS_ECHO(["----------------------------------- SUMMARY -----------------------------------"])
AS_ECHO([])
AS_ECHO(["Package version............... : $PACKAGE-$VERSION"])
AS_ECHO([])
AS_ECHO(["C++ compiler.................. : $CXX"])
AS_ECHO(["Build Methods...................... : $METHODS"])
for method in ${METHODS}; do
     AS_CASE("${method}",
             [opt],   [AS_ECHO(["CPPFLAGS...(opt)................... : $CPPFLAGS_OPT $CPPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(opt)................... : $CXXFLAGS_OPT $CXXFLAGS"])],
             [devel], [AS_ECHO(["CPPFLAGS...(devel)................. : $CPPFLAGS_DEVEL $CPPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(devel)................. : $CXXFLAGS_DEVEL $CXXFLAGS"])],
             [dbg],   [AS_ECHO(["CPPFLAGS...(dbg)................... : $CPPFLAGS_DBG $CPPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(dbg)................... : $CXXFLAGS_DBG $CXXFLAGS"])],
             [prof],  [AS_ECHO(["CPPFLAGS...(prof).................. : $CPPFLAGS_PROF $CPPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(prof).................. : $CXXFLAGS_PROF $CXXFLAGS"])],
             [oprof], [AS_ECHO(["CPPFLAGS...(oprof)................. : $CPPFLAGS_OPROF $CPPFLAGS"])
                       AS_ECHO(["CXXFLAGS...(oprof)................. : $CXXFLAGS_OPROF $CXXFLAGS"])])
     AS_ECHO([])
done

AS_ECHO(["Install dir................... : $prefix"])
AS_ECHO(["Build user.................... : $USER"])
AS_ECHO(["Build host.................... : $BUILD_HOST"])
AS_ECHO(["Configure date................ : $BUILD_DATE"])
AS_ECHO(["Build architecture............ : $BUILD_ARCH"])
AS_ECHO(["Git revision number........... : $BUILD_VERSION"])
AS_ECHO([])
AS_ECHO(["-------------------------------------------------------------------------------"])

AS_ECHO(["Optional Packages for Testing:"])
AS_IF([test "x$enablempi" = "xyes"],
      [
        AS_ECHO(["  MPI......................... : yes"])
      ],
      [
        AS_ECHO(["  MPI......................... : no"])
      ])
AS_IF([test "x$timpi_optional_INCLUDES" != "x"],
      [
        AS_ECHO(["  timpi_optional_INCLUDES..... : $timpi_optional_INCLUDES"])
      ])
AS_IF([test "x$timpi_optional_LIBS" != "x"],
      [
        AS_ECHO(["  timpi_optional_LIBS......... : $timpi_optional_LIBS"])
      ])


AS_ECHO([])
AS_ECHO(["Configure complete, now type \'make\' and then \'make install\'."])
AS_ECHO([])

])
