IF(BUILD_SHARED_LIBS AND WIN32)
  remove_definitions(-DDLL_EXPORT)
  remove_definitions(-DDLL_NETCDF)
ENDIF()

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ctest.c
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/ref_ctest.c"
  "${CMAKE_CURRENT_BINARY_DIR}/ctest.c"
  )
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ctest64.c
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/ref_ctest64.c"
  "${CMAKE_CURRENT_BINARY_DIR}/ctest64.c"
  )

SET(ncdump_FILES ncdump.c vardata.c dumplib.c indent.c nctime0.c utils.c nciter.c)
SET(nccopy_FILES nccopy.c nciter.c chunkspec.c utils.c dimmap.c)

IF(USE_X_GETOPT)
  SET(ncdump_FILES ${ncdump_FILES} XGetopt.c)
  SET(nccopy_FILES ${nccopy_FILES} XGetopt.c)
ENDIF()

ADD_EXECUTABLE(ncdump ${ncdump_FILES})
ADD_EXECUTABLE(nccopy ${nccopy_FILES})

TARGET_LINK_LIBRARIES(ncdump netcdf ${ALL_TLL_LIBS})
TARGET_LINK_LIBRARIES(nccopy netcdf ${ALL_TLL_LIBS})

####
# We have to do a little tweaking
# to remove the Release/ and Debug/ directories
# in MSVC builds.  This is required to get
# test scripts to work.
####
IF(MSVC)
  SET_TARGET_PROPERTIES(ncdump PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(ncdump PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(ncdump PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
    ${CMAKE_CURRENT_BINARY_DIR})

  SET_TARGET_PROPERTIES(nccopy PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(nccopy PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(nccopy PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
    ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()


IF(ENABLE_TESTS)
  ADD_EXECUTABLE(rewrite-scalar rewrite-scalar.c)
  ADD_EXECUTABLE(bom bom.c)
  ADD_EXECUTABLE(tst_dimsizes tst_dimsizes.c)
  TARGET_LINK_LIBRARIES(rewrite-scalar netcdf)
  TARGET_LINK_LIBRARIES(bom netcdf)
  TARGET_LINK_LIBRARIES(tst_dimsizes netcdf)

IF(USE_NETCDF4)
  ADD_EXECUTABLE(tst_fileinfo tst_fileinfo.c)
  TARGET_LINK_LIBRARIES(tst_fileinfo netcdf)
  add_sh_test(ncdump tst_fileinfo)
ENDIF()

IF(MSVC)
  SET_TARGET_PROPERTIES(rewrite-scalar PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(rewrite-scalar PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(rewrite-scalar PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
    ${CMAKE_CURRENT_BINARY_DIR})

  SET_TARGET_PROPERTIES(bom PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(bom PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(bom PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
    ${CMAKE_CURRENT_BINARY_DIR})

  SET_TARGET_PROPERTIES(tst_dimsizes PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(tst_dimsizes PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(tst_dimsizes PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(tst_fileinfo PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
    ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()

  # Base tests
  # The tests are set up as a combination of shell scripts and executables that
  # must be run in a particular order. It is painful but will use macros to help
  # keep it from being too bad.

  ## Start adding tests in the appropriate order
  add_sh_test(ncdump run_tests)
  add_sh_test(ncdump tst_64bit)
  add_bin_test(ncdump ctest)
  add_bin_test(ncdump ctest64)
  add_sh_test(ncdump tst_output)
  add_sh_test(ncdump tst_lengths)
  add_sh_test(ncdump tst_calendars)
  add_bin_test(ncdump tst_utf8)
  add_sh_test(ncdump run_utf8_tests)
  IF(USE_NETCDF4)
    add_sh_test(ncdump run_utf8_nc4_tests)
  ENDIF(USE_NETCDF4)

  add_sh_test(ncdump tst_nccopy3)
  add_sh_test(ncdump tst_charfill)

  add_sh_test(ncdump tst_formatx3)
  add_sh_test(ncdump tst_bom)
  add_sh_test(ncdump tst_dimsizes)

  # The following test script invokes
  # gcc directly.
  IF(CMAKE_COMPILER_IS_GNUCC OR APPLE)
    IF(ENABLE_LARGE_FILE_TESTS)
      add_sh_test(ncdump tst_iter)
    ENDIF(ENABLE_LARGE_FILE_TESTS)
  ENDIF(CMAKE_COMPILER_IS_GNUCC OR APPLE)

  IF(BUILD_DISKLESS)
    add_sh_test(ncdump tst_inmemory_nc3)
    IF(USE_NETCDF4)
      add_sh_test(ncdump tst_inmemory_nc4)
    ENDIF(USE_NETCDF4)
  ENDIF(BUILD_DISKLESS)

  IF(USE_NETCDF4)
    add_bin_test(ncdump tst_create_files)
    add_bin_test(ncdump tst_group_data)
    add_bin_test(ncdump tst_enum_data)
    add_bin_test(ncdump tst_opaque_data)
    add_bin_test(ncdump tst_string_data)
    add_bin_test(ncdump tst_vlen_data)
    add_bin_test(ncdump tst_comp)
    add_bin_test(ncdump tst_comp2)
    add_bin_test(ncdump tst_nans)
    add_bin_test(ncdump tst_h_scalar)
    add_bin_test(ncdump tst_bug324)
    add_sh_test(ncdump tst_formatx4)
    # Add this test by hand, as it is also called from a script.
    # Editing the script would break autotools compatibility.
    add_bin_test_no_prefix(tst_special_atts)
    add_bin_test_no_prefix(tst_compress)
    add_bin_test_no_prefix(tst_chunking)


    ###
    # This test fails on Visual Studio builds with bash.
    # It passes, but technically fails because the scientific
    # formatting omits a 0.
    ###
    IF(EXTRA_TESTS AND NOT MSVC)
      add_sh_test(ncdump run_back_comp_tests)
    ENDIF()

    # Known failure on MSVC; the number of 0's padding
    # is different, but the result is actually correct.
    #IF(NOT MSVC)
    add_sh_test(ncdump tst_netcdf4)
    add_bin_test(ncdump tst_h_rdc0)
    #ENDIF()

    add_bin_test(ncdump tst_unicode)

    add_bin_test(ncdump tst_fillbug)
    add_sh_test(ncdump_sh tst_fillbug)

    IF(NOT MSVC)
      add_sh_test(ncdump tst_netcdf4_4)
    ENDIF()



    IF(NOT MSVC)
      add_sh_test(ncdump tst_nccopy4)
    ENDIF()

    add_sh_test(ncdump tst_grp_spec)
    add_sh_test(ncdump tst_mud)
    add_sh_test(ncdump_shell tst_h_scalar)

  ENDIF()

  add_sh_test(ncdump tst_ncgen4_classic)
  IF(USE_NETCDF4)
    add_sh_test(ncdump tst_ncgen4)
  ENDIF()

  add_sh_test(ncdump tst_inttags)
  IF(USE_NETCDF4)
    add_sh_test(ncdump tst_inttags4)
  ENDIF()


ENDIF()

IF(MSVC)
  SET_TARGET_PROPERTIES(ncdump
    PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
    )
  SET_TARGET_PROPERTIES(nccopy
    PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
    )

ENDIF()


INSTALL(TARGETS ncdump RUNTIME DESTINATION bin COMPONENT utilities)
INSTALL(TARGETS nccopy RUNTIME DESTINATION bin COMPONENT utilities)
SET(MAN_FILES nccopy.1 ncdump.1)


FILE(GLOB COPY_FILES ${CMAKE_BINARY_DIR}/ncgen/*.nc ${CMAKE_BINARY_DIR}/nc_test4/*.nc ${CMAKE_CURRENT_SOURCE_DIR}/*.ncml ${CMAKE_CURRENT_SOURCE_DIR}/*.nc ${CMAKE_CURRENT_SOURCE_DIR}/*.cdl ${CMAKE_CURRENT_SOURCE_DIR}/*.sh ${CMAKE_CURRENT_SOURCE_DIR}/*.1)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ FILE_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE)

ADD_SUBDIRECTORY(cdl)
ADD_SUBDIRECTORY(expected)

SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")

IF(NOT MSVC)
  INSTALL(FILES ${MAN_FILES} DESTINATION "share/man/man1" COMPONENT documentation)
ENDIF()
