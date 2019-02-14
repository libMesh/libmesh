SET(abs_top_srcdir ${CMAKE_CURRENT_SOURCE_DIR})

remove_definitions(-DDLL_EXPORT)

ADD_SUBDIRECTORY(baseline)
ADD_SUBDIRECTORY(baselineraw)
ADD_SUBDIRECTORY(baselineremote)
ADD_SUBDIRECTORY(cdltestfiles)
ADD_SUBDIRECTORY(daptestfiles)
ADD_SUBDIRECTORY(dmrtestfiles)
ADD_SUBDIRECTORY(nctestfiles)
ADD_SUBDIRECTORY(misctestfiles)

FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)

FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ FILE_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE)

IF(ENABLE_TESTS)

  # Base tests
  # The tests are set up as a combination of shell scripts and executables that
  # must be run in a particular order. It is painful but will use macros to help
  # keep it from being too bad.

  IF(BUILD_UTILITIES)
    add_sh_test(dap4_test test_raw)
  ENDIF(BUILD_UTILITIES)

  IF(ENABLE_DAP_REMOTE_TESTS)

    # Change name (add '4') to avoid cmake
    # complaint about duplicate targets.
    BUILD_BIN_TEST(findtestserver4)

    IF(BUILD_UTILITIES)
      add_sh_test(dap4_test test_remote)
    ENDIF(BUILD_UTILITIES)
  ENDIF(ENABLE_DAP_REMOTE_TESTS)
ENDIF(ENABLE_TESTS)

FILE(COPY ./baselineraw DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
#FILE(COPY ./baseline DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
#FILE(COPY ./cdltestfiles DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
#FILE(COPY ./daptestfiles DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
#FILE(COPY ./dmrtestfiles DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

## Specify files to be distributed by 'make dist'
FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.sh
#${CMAKE_CURRENT_SOURCE_DIR}/daptestfiles
#${CMAKE_CURRENT_SOURCE_DIR}/dmrtestfiles
#${CMAKE_CURRENT_SOURCE_DIR}/cdltestfiles
#${CMAKE_CURRENT_SOURCE_DIR}/baseline
${CMAKE_CURRENT_SOURCE_DIR}/baselineraw
)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} CMakeLists.txt Makefile.am)
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")
