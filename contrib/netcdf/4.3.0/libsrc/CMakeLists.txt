SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")


IF(NOT MSVC)
    # Process these files with m4.
    SET(m4_SOURCES attr ncx putget t_ncxx)
    foreach (f ${m4_SOURCES})      
    	       GEN_m4(${f})
    endforeach(f)
ENDIF()

SET(libsrc_SOURCES v1hpg.c putget.c attr.c nc3dispatch.c nc3internal.c var.c dim.c ncx.c lookup3.c ncio.c)

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
set(man_MANS netcdf.3)

SET(MAN_FILES netcdf.3)
IF(EXISTS ${MAN_FILES} AND NOT MSVC)
	INSTALL(FILES ${MAN_FILES} DESTINATION "share/man/man3" COMPONENT documentation)
ENDIF()

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

