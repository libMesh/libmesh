## This is a CMake file, part of Unidata's netCDF package.
# Copyright 2012-2013, see the COPYRIGHT file for more information.

SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")


# Process these files with m4.
SET(m4_SOURCES attr ncx putget t_ncxx)
foreach (f ${m4_SOURCES})      
  GEN_m4(${f})
endforeach(f)

SET(libsrc_SOURCES v1hpg.c putget.c attr.c nc3dispatch.c nc3internal.c var.c dim.c ncx.c lookup3.c ncio.c)
SET(libsrc_SOURCES ${libsrc_SOURCES} pstdint.h)

IF (BUILD_DISKLESS)
  SET(libsrc_SOURCES ${libsrc_SOURCES} memio.c)

  IF (BUILD_MMAP)
    SET(libsrc_SOURCES ${libsrc_SOURCES} mmapio.c)
  ENDIF( BUILD_MMAP)
ENDIF (BUILD_DISKLESS)

IF (USE_FFIO)
  SET(libsrc_SOURCES ${libsrc_SOURCES} ffio.c)
ELSE (USE_FFIO)
  SET(libsrc_SOURCES ${libsrc_SOURCES} posixio.c)
ENDIF (USE_FFIO)

add_library(netcdf3 OBJECT ${libsrc_SOURCES}) 

# The C API man page.

SET(MAN_FILES netcdf.3)

set(ARGS_MANPAGE -DAPI=C)

IF (USE_NETCDF4)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DNETCDF4=TRUE)
ENDIF(USE_NETCDF4)

IF (BUILD_DAP)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DDAP=TRUE)
ENDIF(BUILD_DAP)

IF (BUILD_PARALLEL)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DPARALLEL_IO=TRUE)
ENDIF(BUILD_PARALLEL)

# Generate the man page.
IF(NOT MSVC)

  FILE(GLOB libsrc_MANPAGE ${CMAKE_SOURCE_DIR}/man4/netcdf.m4)
  FILE(COPY ${libsrc_MANPAGE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  ADD_CUSTOM_TARGET(manpage ALL
    COMMAND m4 ${ARGS_MANPAGE} ${CMAKE_CURRENT_BINARY_DIR}/netcdf.m4 > ${CMAKE_CURRENT_BINARY_DIR}/netcdf.3
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/netcdf.3 DESTINATION "share/man/man3" COMPONENT documentation)
ENDIF()

# Add files to the distribution
FILE(GLOB libsrc_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
SET(CUR_EXTRA_DIST attr.m4 ncx.m4 putget.m4 XGetopt.c memio.c mmapio.c ffio.c posixio.c CMakeLists.txt Makefile.am)
ADD_EXTRA_DIST("${libsrc_SOURCES}")
ADD_EXTRA_DIST("${libsrc_HEADERS}") 
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")

