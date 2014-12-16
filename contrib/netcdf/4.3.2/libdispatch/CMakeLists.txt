SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

SET(libdispatch_SOURCES dparallel.c dcopy.c dfile.c ddim.c datt.c dattinq.c dattput.c dattget.c derror.c dvar.c dvarget.c dvarput.c dvarinq.c ddispatch.c nclog.c dstring.c dutf8proc.c ncuri.c nclist.c ncbytes.c nchashmap.c nctime.c dsubstrate.c nc.c nclistmgr.c)

IF(USE_NETCDF4)
  SET(libdispatch_SOURCES ${libdispatch_SOURCES} dgroup.c dvlen.c dcompound.c dtype.c denum.c dopaque.c ncaux.c)
ENDIF(USE_NETCDF4)

IF(BUILD_V2)
  SET(libdispatch_SOURCES ${libdispatch_SOURCES} dv2i.c)
ENDIF(BUILD_V2)	

add_library(dispatch OBJECT ${libdispatch_SOURCES}) 

FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} CMakeLists.txt Makefile.am)
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")