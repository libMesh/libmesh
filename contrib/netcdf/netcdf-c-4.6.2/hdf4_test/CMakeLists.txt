# Copy some test files from current source dir to out-of-tree build dir.
FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.sh ${CMAKE_CURRENT_SOURCE_DIR}/*.hdf4)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
IF(MSVC)
  FILE(COPY ${COPY_FILES} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY}/)
ENDIF()

IF(USE_HDF4_FILE_TESTS AND NOT MSVC)
  build_bin_test_no_prefix(tst_interops2)
  build_bin_test_no_prefix(tst_interops3)
  add_bin_test(hdf4_test tst_chunk_hdf4)
  add_bin_test(hdf4_test tst_h4_lendian)
  add_bin_test(hdf4_test tst_hdf4_extra)
  add_sh_test(hdf4_test run_get_hdf4_files)
  add_sh_test(hdf4_test run_formatx_hdf4)
ENDIF()

