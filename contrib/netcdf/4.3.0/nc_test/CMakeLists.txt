SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")
SET (nc_test_m4_SOURCES test_get test_put)
FOREACH (F ${nc_test_m4_SOURCES})
	GEN_m4(${F})
ENDFOREACH()


SET (nc_test_SRC nc_test.c error.c test_get.c test_put.c test_read.c test_write.c util.c)
ADD_EXECUTABLE(nc_test ${nc_test_SRC})
TARGET_LINK_LIBRARIES(nc_test netcdf)

# Some extra tests
SET(TESTS t_nc tst_small tst_misc tst_norm tst_names tst_nofill tst_nofill2 tst_nofill3)

IF(NOT MSVC)
	SET(TESTS ${TESTS} tst_atts3)
ENDIF()

IF(USE_NETCDF4)
  SET(TESTS ${TESTS} tst_atts)
ENDIF()

IF(LARGE_FILE_TESTS)
  SET(TESTS ${TESTS} quick_large_files tst_big_var6 tst_big_var2 tst_big_rvar tst_big_var tst_large large_files)
ENDIF()

IF(BUILD_DISKLESS)
  SET(TESTFILES ${TESTFILES} tst_diskless tst_diskless3 tst_diskless4)
  IF(USE_NETCDF4)
    SET(TESTFILES ${TESTFILES} tst_diskless2)
  ENDIF()
ENDIF()

# Build executables required for the shell scripts.
FOREACH(BINFILE ${TESTFILES})
	build_bin_test(${BINFILE})
ENDFOREACH()

# Process single-file test files.
FOREACH(CTEST ${TESTS})
	add_bin_test(nc_test ${CTEST})
ENDFOREACH()

ADD_TEST(nc_test ${EXECUTABLE_OUTPUT_PATH}/nc_test)

IF(BUILD_DISKLESS)
	add_sh_test(nc_test run_diskless)  	
	IF(BUILD_MMAP)
		add_sh_test(nc_test run_mmap)
	ENDIF()
	IF(LARGE_FILE_TESTS)
		add_sh_test(nc_test run_diskless2.sh)
	ENDIF()

ENDIF()

# Copy some test files from current source dir to out-of-tree build dir.
FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.nc ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

