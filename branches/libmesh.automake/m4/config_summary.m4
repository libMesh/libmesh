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
echo C++ compiler.................. : $CXX
echo C++ compiler type..............: $GXX_VERSION
echo CXXFLAGS-opt...................: $CXXFLAGS_OPT
echo CXXFLAGS-devel.................: $CXXFLAGS_DVL
echo CXXFLAGS-dbg...................: $CXXFLAGS_DBG
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
# echo
# echo External Packages:

# if test "x$enable_cantera" = "xyes"; then
#   echo '   'Link with Cantera.......... : yes
#   echo '      'CANTERA_INCLUDE......... : $CANTERA_INCLUDE
#   echo '      'CANTERA_LIB............. : $CANTERA_LIB
#   echo 
# else
#   echo '   'Link with Cantera.......... : no
# fi

# if test "x$HAVE_CHEMAY" = "x1"; then
#   echo '   'Link with Chemay........... : yes
#   echo '      'CHEMAY_CPPFLAGS......... : $CHEMAY_CPPFLAGS
#   echo '      'CHEMAY_LIBS............. : $CHEMAY_LIBS
#   echo 
# else
#   echo '   'Link with Chemay........... : no
# fi

# if test "x$ENABLE_ABLATION" = "x1"; then
#   echo '   'PECOS Ablating surfaces.... : yes
#   echo '      'ABLATION_CFLAGS......... : $ABLATION_CFLAGS
#   echo '      'ABLATION_LIBS........... : $ABLATION_LIBS
#   echo 
# else
#   echo '   'PECOS Ablating surfaces.... : no
# fi

# if test "x$enable_boost" = "xyes"; then
#   echo '   'Link with Boost............ : yes
#   echo '      'BOOST_CPPFLAGS.......... : $BOOST_CPPFLAGS
#   echo '      'Boost LD Flags.......... : $BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS
#   echo '      'Boost Unit Test......... : yes
#   echo
# else
#   echo '   'Link with Boost............ : no

# fi	

# if test "x$HAVE_MASA" = "x1"; then
#   echo '   'Link with MASA............. : yes
#   echo '      'MASA_CXXFLAGS........... : $MASA_CXXFLAGS
#   echo '      'MASA_LIBS............... : $MASA_LIBS
#   echo
# else
#   echo '   'Link with MASA............. : no
# fi

# if test "x$HAVE_GCOV_TOOLS" = "x1"; then
#   echo '   'Enable gcov code coverage.. : yes
# else
#   echo '   'Enable gcov code coverage.. : no
# fi

		   
echo
echo '-------------------------------------------------------------------------------'

echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
