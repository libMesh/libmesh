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

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo C++ compiler.................. : $CXX
echo C++ compiler flags............ : $CXXFLAGS
echo Install dir................... : $prefix 
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo Git revision number........... : $BUILD_VERSION
echo
echo '-------------------------------------------------------------------------------'

echo Optional Packages for Testing:
if test "x$enable_mpi" = "x1"; then
  echo '  'MPI......................... : yes
  echo '  'MPI_INCLUDES_PATH........... : $MPI_INCLUDES_PATH
  echo '  'MPI_LIBS_PATH............... : $MPI_LIBS_PATH
  echo '  'MPI_LIBS.................... : $MPI_LIBS
else
  echo '  'MPI......................... : no
fi


echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
