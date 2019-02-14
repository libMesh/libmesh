IF(BUILD_SHARED_LIBS AND WIN32)
  remove_definitions(-DDLL_EXPORT)
  remove_definitions(-DDLL_NETCDF)
ENDIF()

SET(ncdump_FILES ncdump.c vardata.c dumplib.c indent.c nctime0.c utils.c nciter.c)
SET(nccopy_FILES nccopy.c nciter.c chunkspec.c utils.c dimmap.c list.c)
SET(ocprint_FILES ocprint.c)

IF(USE_X_GETOPT)
  SET(ncdump_FILES ${ncdump_FILES} XGetopt.c)
  SET(nccopy_FILES ${nccopy_FILES} XGetopt.c)
  SET(ocprint_FILES ${ocprint_FILES} XGetopt.c)

ENDIF()

ADD_EXECUTABLE(ncdump ${ncdump_FILES})
ADD_EXECUTABLE(nccopy ${nccopy_FILES})

IF(ENABLE_DAP)
  ADD_EXECUTABLE(ocprint ${ocprint_FILES})
ENDIF(ENABLE_DAP)

TARGET_LINK_LIBRARIES(ncdump netcdf ${ALL_TLL_LIBS})
TARGET_LINK_LIBRARIES(nccopy netcdf ${ALL_TLL_LIBS})

IF(ENABLE_DAP)
  TARGET_LINK_LIBRARIES(ocprint netcdf ${ALL_TLL_LIBS})
ENDIF(ENABLE_DAP)

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

  IF(ENABLE_DAP)
    SET_TARGET_PROPERTIES(ocprint PROPERTIES RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(ocprint PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(ocprint PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
      ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF(ENABLE_DAP)

ENDIF()


IF(ENABLE_TESTS)
  ADD_EXECUTABLE(rewrite-scalar rewrite-scalar.c)
  ADD_EXECUTABLE(bom bom.c)
  ADD_EXECUTABLE(tst_dimsizes tst_dimsizes.c)
  ADD_EXECUTABLE(nctrunc nctrunc.c)
  TARGET_LINK_LIBRARIES(rewrite-scalar netcdf)
  TARGET_LINK_LIBRARIES(bom netcdf)
  TARGET_LINK_LIBRARIES(tst_dimsizes netcdf)

  IF(USE_NETCDF4)
    ADD_EXECUTABLE(tst_fileinfo tst_fileinfo.c)
    TARGET_LINK_LIBRARIES(tst_fileinfo netcdf)
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

    SET_TARGET_PROPERTIES(nctrunc PROPERTIES RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(nctrunc PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(nctrunc PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
      ${CMAKE_CURRENT_BINARY_DIR})

    IF(USE_NETCDF4)
      SET_TARGET_PROPERTIES(tst_fileinfo PROPERTIES RUNTIME_OUTPUT_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR})
      SET_TARGET_PROPERTIES(tst_fileinfo PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
        ${CMAKE_CURRENT_BINARY_DIR})
      SET_TARGET_PROPERTIES(tst_fileinfo PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
        ${CMAKE_CURRENT_BINARY_DIR})

  ENDIF(USE_NETCDF4)
ENDIF(MSVC)

  # Base tests
  # The tests are set up as a combination of shell scripts and executables that
  # must be run in a particular order. It is painful but will use macros to help
  # keep it from being too bad.

  ## Start adding tests in the appropriate order
  add_sh_test(ncdump run_tests)
  add_sh_test(ncdump tst_64bit)
  add_bin_test_no_prefix(ref_ctest)
  add_bin_test_no_prefix(ref_ctest64)
  add_sh_test(ncdump tst_output)
  add_sh_test(ncdump tst_lengths)
  add_sh_test(ncdump tst_calendars)
  build_bin_test_no_prefix(tst_utf8)
  add_sh_test(ncdump run_utf8_tests)
  IF(USE_NETCDF4)
    add_sh_test(ncdump run_utf8_nc4_tests)
    add_sh_test(ncdump tst_fileinfo)
    add_sh_test(ncdump tst_hdf5_offset)
  ENDIF(USE_NETCDF4)

  add_sh_test(ncdump tst_null_byte_padding)
  IF(USE_STRICT_NULL_BYTE_HEADER_PADDING)
    SET_TESTS_PROPERTIES(ncdump_tst_null_byte_padding PROPERTIES WILL_FAIL TRUE)
  ENDIF(USE_STRICT_NULL_BYTE_HEADER_PADDING)

  add_sh_test(ncdump tst_nccopy3)
  IF(HAVE_BASH)
    SET_TESTS_PROPERTIES(ncdump_tst_nccopy3 PROPERTIES RUN_SERIAL TRUE)
  ENDIF(HAVE_BASH)

  add_sh_test(ncdump tst_nccopy3_subset)
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

  add_sh_test(ncdump tst_inmemory_nc3)
  IF(USE_NETCDF4)
      add_sh_test(ncdump tst_inmemory_nc4)
  ENDIF(USE_NETCDF4)

  IF(USE_NETCDF4)
    build_bin_test_no_prefix(tst_create_files)
    build_bin_test_no_prefix(tst_group_data)
    build_bin_test_no_prefix(tst_enum_data)
    build_bin_test_no_prefix(tst_opaque_data)
    build_bin_test_no_prefix(tst_string_data)
    build_bin_test_no_prefix(tst_vlen_data)
    build_bin_test_no_prefix(tst_comp)
    build_bin_test_no_prefix(tst_comp2)
    build_bin_test_no_prefix(tst_nans)
    build_bin_test_no_prefix(tst_h_scalar)
    add_sh_test(ncdump tst_formatx4)
    # Add this test by hand, as it is also called from a script.
    # Editing the script would break autotools compatibility.
    build_bin_test_no_prefix(tst_special_atts)
    build_bin_test_no_prefix(tst_compress)
    build_bin_test_no_prefix(tst_chunking)

    ###
    # This test fails on Visual Studio builds with bash.
    # It passes, but technically fails because the scientific
    # formatting omits a 0.
    ###
    IF(EXTRA_TESTS)
      add_sh_test(ncdump run_back_comp_tests)
      IF(MSVC)
        SET_TESTS_PROPERTIES(ncdump_run_back_comp_tests PROPERTIES WILL_FAIL TRUE)
      ENDIF(MSVC)
    ENDIF(EXTRA_TESTS)

    # Known failure on MSVC; the number of 0's padding
    # is different, but the result is actually correct.
    IF(HAVE_BASH)
      add_sh_test(ncdump tst_netcdf4)
      build_bin_test_no_prefix(tst_h_rdc0)
    ENDIF()

    build_bin_test_no_prefix(tst_unicode)

    build_bin_test_no_prefix(tst_fillbug)
    add_sh_test(ncdump_sh tst_fillbug)

    add_sh_test(ncdump tst_netcdf4_4)
    IF(MSVC AND HAVE_BASH)
      SET_TESTS_PROPERTIES(ncdump_tst_netcdf4_4 PROPERTIES WILL_FAIL TRUE)
    ENDIF(MSVC AND HAVE_BASH)

    ###
    # Some test reordering was required to ensure these tests
    # only ran after ncdump was built.
    ###
    add_sh_test(ncdump run_ncgen_tests)
    IF(USE_NETCDF4)
      add_sh_test(ncdump run_ncgen_nc4_tests)
    ENDIF(USE_NETCDF4)

    add_sh_test(ncdump tst_grp_spec)
    add_sh_test(ncdump tst_mud)
    add_sh_test(ncdump_shell tst_h_scalar)

    add_sh_test(ncdump tst_nccopy4)
    IF(HAVE_BASH)
      SET_TESTS_PROPERTIES(ncdump_tst_nccopy4 PROPERTIES RUN_SERIAL TRUE)
    ENDIF(HAVE_BASH)

    add_sh_test(ncdump tst_nccopy5)
    IF(HAVE_BASH)
      SET_TESTS_PROPERTIES(ncdump_tst_nccopy5 PROPERTIES RUN_SERIAL TRUE)
    ENDIF(HAVE_BASH)

    IF(USE_NETCDF4)
      add_sh_test(ncdump tst_ncgen4)
    ENDIF(USE_NETCDF4)

  add_sh_test(ncdump tst_inttags)
  IF(USE_NETCDF4)
    add_sh_test(ncdump tst_inttags4)
  ENDIF(USE_NETCDF4)

ENDIF()

ENDIF()

IF(MSVC)
  SET_TARGET_PROPERTIES(ncdump
    PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
    )
  SET_TARGET_PROPERTIES(nccopy
    PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
    )

  IF(ENABLE_DAP)
    SET_TARGET_PROPERTIES(ocprint
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF(ENABLE_DAP)

ENDIF()


INSTALL(TARGETS ncdump RUNTIME DESTINATION bin COMPONENT utilities)
INSTALL(TARGETS nccopy RUNTIME DESTINATION bin COMPONENT utilities)

SET(MAN_FILES nccopy.1 ncdump.1)

# Note, the L512.bin file is file containing exactly 512 bytes each of value 0.
# It is used for creating hdf5 files with varying offsets for testing.

FILE(GLOB COPY_FILES ${CMAKE_BINARY_DIR}/ncgen/*.nc ${CMAKE_BINARY_DIR}/nc_test4/*.nc ${CMAKE_CURRENT_SOURCE_DIR}/*.ncml ${CMAKE_CURRENT_SOURCE_DIR}/*.nc ${CMAKE_CURRENT_SOURCE_DIR}/*.cdl ${CMAKE_CURRENT_SOURCE_DIR}/*.sh ${CMAKE_CURRENT_SOURCE_DIR}/*.1 ${CMAKE_CURRENT_SOURCE_DIR}/L512.bin)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ FILE_PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE)

ADD_SUBDIRECTORY(cdl)
ADD_SUBDIRECTORY(expected)

SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")

IF(NOT MSVC)
  INSTALL(FILES ${MAN_FILES} DESTINATION "share/man/man1" COMPONENT documentation)
ENDIF()
