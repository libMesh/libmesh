# Set Mangling by hand 
# Assume _name mangling see src/exo_jack.c code for how ADDC_ will be used in
# Fortran mangling
add_definitions (-DADDC_)


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

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_BINARY_DIR}/include"
  "${NETCDF_INCLUDE_DIR}"
)

SET(forbind_SRCS
   src/exo_jack.c
)

ADD_LIBRARY(exoIIv2for ${forbind_SRCS})

TARGET_LINK_LIBRARIES(exoIIv2for exoIIv2c ${NETCDF_LIBRARY} ${HDF5_LIBRARY} ${HDF5HL_LIBRARY})

IF(EXODUS_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(exoIIv2for PROPERTIES ${EXODUS_LIBRARY_PROPERTIES})
ENDIF(EXODUS_LIBRARY_PROPERTIES)

IF(NOT EXODUS_INSTALL_NO_LIBRARIES)
  INSTALL(TARGETS exoIIv2for
    RUNTIME DESTINATION lib COMPONENT Runtime
    LIBRARY DESTINATION lib COMPONENT Runtime
    ARCHIVE DESTINATION lib COMPONENT Development)
  INSTALL(FILES ${forbind_HEADERS} DESTINATION include COMPONENT Development)
ENDIF(NOT EXODUS_INSTALL_NO_LIBRARIES)


INCLUDE (${CMAKE_ROOT}/Modules/Dart.cmake)

IF (BUILD_TESTING)
  ENABLE_TESTING()
  SUBDIRS(test)
ENDIF (BUILD_TESTING)

