SET(dap4_SOURCES d4crc32.c d4curlfunctions.c d4fix.c d4data.c d4file.c d4parser.c d4meta.c d4varx.c d4dump.c d4swap.c d4chunk.c d4printer.c d4read.c d4http.c d4util.c d4odom.c d4cvt.c d4debug.c ncd4dispatch.c ezxml_extra.c ezxml.c)

add_library(dap4 OBJECT ${dap4_SOURCES})

###
# Options related to the man page generation.
###

###
# The C API man page.
###

set(man_MANS netcdf.3)

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

FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} dce.y CMakeLists.txt Makefile.am)
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")

