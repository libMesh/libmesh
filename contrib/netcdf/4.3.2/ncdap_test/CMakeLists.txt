SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")


SET(abs_top_srcdir ${CMAKE_CURRENT_SOURCE_DIR})

remove_definitions(-DDLL_EXPORT)
add_definitions(-D"TOPSRCDIR=${CMAKE_SOURCE_DIR}")

FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)

FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ FILE_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE) 

IF(ENABLE_DAP_REMOTE_TESTS)
  ADD_EXECUTABLE(nctestserver nctestserver.c)
  TARGET_LINK_LIBRARIES(nctestserver netcdf ${ALL_TLL_LIBS})
ENDIF()

IF(ENABLE_TESTS)
  # Base tests
  # The tests are set up as a combination of shell scripts and executables that
  # must be run in a particular order. It is painful but will use macros to help
  # keep it from being too bad.
  # Binary Test Macro
  
  add_sh_test(ncdap tst_ncdap3)
  
  IF(NOT MSVC)
    add_bin_env_test(ncdap t_dap3a)
    add_bin_env_test(ncdap test_cvt)
    add_bin_env_test(ncdap test_vara)
  ENDIF()

  IF(ENABLE_DAP_REMOTE_TESTS)
    add_sh_test(ncdap tst_remote3)
    add_sh_test(ncdap tst_formatx)
    add_bin_test(ncdap test_partvar)
    IF(ENABLE_DAP_LONG_TESTS)
      add_sh_test(ncdap tst_longremote3)
    ENDIF()
    
    add_bin_test(ncdap test_nstride_cached)
    add_bin_test(ncdap t_misc)
    
    IF(ENABLE_DAP_AUTH_TESTS)
      add_bin_test(ncdap t_auth)
    ENDIF()
    
    add_bin_test(ncdap test_varm3)
    
  ENDIF()

ENDIF()
ADD_SUBDIRECTORY(testdata3)
ADD_SUBDIRECTORY(expected3)
ADD_SUBDIRECTORY(expected4)
ADD_SUBDIRECTORY(expectremote3)
ADD_SUBDIRECTORY(expectremote4)


## Specify files to be distributed by 'make dist'
FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} CMakeLists.txt Makefile.am)
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")
