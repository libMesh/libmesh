cmake_minimum_required(VERSION 2.6)
PROJECT(ExodusII C CXX Fortran)

SET(EXODUSII_VERSION_MAJOR "5")
SET(EXODUSII_VERSION_MINOR "14")
SET(EXODUSII_VERSION_PATCH "0")
SET(EXODUSII_VERSION "${EXODUSII_VERSION_MAJOR}.${EXODUSII_VERSION_MINOR}")
SET(EXODUSII_VERSION_FULL "${EXODUSII_VERSION}.${EXODUSII_VERSION_PATCH}")

SET(EXODUSII_BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}")

SUBDIRS(cbind forbind)

find_path( NETCDF_INCLUDE_DIR netcdf.h
   $ENV{ACCESS}/inc
   $ENV{NETCDF_DIR}/inc
   $ENV{NETCDF_DIR}/include
   $ENV{NETCDF_DIR}/libsrc
)	

find_library( NETCDF_LIBRARY
   NAMES netcdf
   PATHS "$ENV{ACCESS}/lib"
         "$ENV{NETCDF_DIR}/lib"
	 "$ENV{NETCDF_DIR}/libsrc/.libs"
)

find_program (NETCDF_NCDUMP ncdump
   $ENV{ACCESS}/bin
   $ENV{NETCDF_DIR}/bin
   $ENV{NETCDF_DIR}/ncdump
)

# Hack for HDF5
find_library( HDF5_LIBRARY hdf5
  $ENV{ACCESS}/lib/shared
  $ENV{ACCESS}/lib
  $ENV{NETCDF_DIR}/lib
  $ENV{NETCDF_DIR}/libsrc/.libs
  /usr/local/hdf5/lib
  )
find_library( HDF5HL_LIBRARY hdf5_hl
  $ENV{ACCESS}/lib/shared
  $ENV{ACCESS}/lib
  $ENV{NETCDF_DIR}/lib
  $ENV{NETCDF_DIR}/libsrc/.libs
  /usr/local/hdf5/lib
  )


