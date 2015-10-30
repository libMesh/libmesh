cmake_minimum_required(VERSION 2.6)
SET(NEMESIS_VERSION_MAJOR "5")
SET(NEMESIS_VERSION_MINOR "15")
SET(NEMESIS_VERSION_PATCH "0")
SET(NEMESIS_VERSION "${NEMESIS_VERSION_MAJOR}.${NEMESIS_VERSION_MINOR}")
SET(NEMESIS_VERSION_FULL "${NEMESIS_VERSION}.${NEMESIS_VERSION_PATCH}")

SET(NEMESIS_BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}")
SET(HEADERS "ne_nemesisI.h")


find_path( NETCDF_INCLUDE_DIR netcdf.h
   $ENV{ACCESS}/inc
   $ENV{NETCDF_DIR}/inc
   $ENV{NETCDF_DIR}/include
   $ENV{NETCDF_DIR}/libsrc
   )   

find_library( NETCDF_LIBRARY netcdf
  $ENV{ACCESS}/lib/shared
  $ENV{ACCESS}/lib
  $ENV{NETCDF_DIR}/lib
  $ENV{NETCDF_DIR}/libsrc/.libs
  )

find_program (NETCDF_NCDUMP ncdump
   $ENV{ACCESS}/bin
   $ENV{NETCDF_DIR}/bin
   $ENV{NETCDF_DIR}/ncdump
   )

find_library( EXODUS_LIBRARY exodus
  $ENV{ACCESS}/lib/shared
  $ENV{ACCESS}/lib
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

find_library(Z_LIBRARY z
  /usr/local/zlib/lib
  )

find_library(MATH_LIBRARY
             NAMES m
             PATHS /usr/lib64
                   /usr/lib
)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_BINARY_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/../exodusii/cbind/include/"
  "${CMAKE_CURRENT_BINARY_DIR}/../exodusii/cbind/include/"
  "${NETCDF_INCLUDE_DIR}")

SET(SOURCES
nemesis_wrapper.c
)

ADD_LIBRARY(nemesis ${SOURCES})
TARGET_LINK_LIBRARIES(nemesis ${EXODUS_LIBRARY} ${NETCDF_LIBRARY} ${HDF5HL_LIBRARY} ${HDF5_LIBRARY} ${Z_LIBRARY} ${MATH_LIBRARY})

ADD_EXECUTABLE(ne_test ne_ctest_wrap.c)
TARGET_LINK_LIBRARIES( ne_test nemesis ${EXODUS_LIBRARY} ${NETCDF_LIBRARY} ${HDF5HL_LIBRARY} ${HDF5_LIBRARY} ${Z_LIBRARY} ${MATH_LIBRARY})
 
INSTALL(TARGETS nemesis
  RUNTIME DESTINATION lib COMPONENT Runtime
  LIBRARY DESTINATION lib COMPONENT Runtime
  ARCHIVE DESTINATION lib COMPONENT Development)
INSTALL(FILES ${HEADERS} DESTINATION include COMPONENT Development)

