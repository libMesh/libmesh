SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

SET(nctest_SRC add.c atttests.c cdftests.c dimtests.c driver.c emalloc.c error.c misctest.c rec.c slabs.c val.c vardef.c varget.c vargetg.c varput.c varputg.c vartests.c vputget.c vputgetg.c)

FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.nc)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

ADD_EXECUTABLE(nctest ${nctest_SRC})
ADD_EXECUTABLE(tst_rename tst_rename.c)

TARGET_LINK_LIBRARIES(nctest netcdf)
TARGET_LINK_LIBRARIES(tst_rename netcdf)

ADD_TEST(nctest ${EXECUTABLE_OUTPUT_PATH}/nctest)
ADD_TEST(tst_rename ${EXECUTABLE_OUTPUT_PATH}/tst_rename)

## Specify files to be distributed by 'make dist'
FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} CMakeLists.txt Makefile.am)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} ref_nctest_classic.nc ref_nctest_64bit_offset.nc)
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")
