# Some extra tests
SET(NC4_TESTS tst_dims tst_dims2 tst_dims3 tst_files tst_files4 tst_vars
  tst_varms tst_unlim_vars tst_converts tst_converts2 tst_grps tst_grps2
  tst_compounds tst_compounds2 tst_compounds3 tst_opaques tst_strings
  tst_strings2 tst_interops tst_interops4 tst_interops6
  tst_enums tst_coords tst_coords2 tst_coords3 tst_vars3 tst_vars4
  tst_chunks tst_chunks2 tst_utf8 tst_fills tst_fills2 tst_fillbug
  tst_xplatform2 tst_h_atts2 tst_endian_fill tst_atts
  t_type cdm_sea_soundings tst_vl tst_atts1 tst_atts2
  tst_vars2 tst_files5 tst_files6 tst_sync tst_h_strbug tst_h_refs
  tst_h_scalar tst_rename tst_h5_endians tst_atts_string_rewrite
  tst_put_vars_two_unlim_dim tst_hdf5_file_compat tst_fill_attr_vanish
  tst_rehash)

# Note, renamegroup needs to be compiled before run_grp_rename

IF(BUILD_UTILITIES)
  SET(NC4_TESTS ${NC4_TESTS} tst_xplatform)
  build_bin_test(renamegroup)
  add_sh_test(nc_test4 run_grp_rename)
  ADD_SH_TEST(nc_test4 tst_misc)
ENDIF(BUILD_UTILITIES)

##
# The shell script, run_empty_vlen_test.sh,
# depends on the 'tst_empty_vlen_unlim' binary.
##
BUILD_BIN_TEST(tst_empty_vlen_unlim)
ADD_SH_TEST(nc_test4 run_empty_vlen_test)



IF(NOT MSVC)
  SET(NC4_TESTS ${NC4_TESTS} tst_interops5 tst_camrun)
ENDIF()

# If the v2 API was built, add the test program.
IF(ENABLE_V2_API)
  build_bin_test(tst_v2)
ENDIF()

IF(LARGE_FILE_TESTS)
  SET(NC4_TESTS ${NC4_TESTS} tst_large tst_large2)
ENDIF()


IF(USE_HDF4)
  SET(NC4_TESTS ${NC4_TESTS} tst_interops2)
ENDIF()

IF(BUILD_BENCHMARKS)
  add_sh_test(nc_test4 run_bm_test1)
  add_sh_test(nc_test4 run_bm_elena)
  add_sh_test(nc_test4 run_bm_test2)
  add_sh_test(nc_test4 run_tst_chunks)
  add_sh_test(nc_test4 run_bm_ar4)
  add_sh_test(nc_test4 run_get_knmi_files)

  SET(NC4_TESTS ${NC4_TESTS} tst_create_files bm_file tst_chunks3 tst_ar4 tst_ar4_3d tst_ar4_4d bm_many_objs tst_h_many_atts bm_many_atts tst_files2 tst_files3 tst_ar5 tst_h_files3 tst_mem tst_knmi bm_netcdf4_recs)
  IF(TEST_PARALLEL)
    add_sh_test(nc_test4 run_par_bm_test)
  ENDIF()
ENDIF()

# Copy some test files from current source dir to out-of-tree build dir.
FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.nc ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
IF(MSVC)
  FILE(COPY ${COPY_FILES} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY}/)
ENDIF()


FOREACH(CTEST ${NC4_TESTS})
  add_bin_test(nc_test4 ${CTEST})
ENDFOREACH()

# This must go after the 'foreach' stanza
# immediately above this comment.
IF(USE_HDF4_FILE_TESTS AND NOT MSVC)
    add_sh_test(nc_test4 run_get_hdf4_files)
    add_bin_test(nc_test4 tst_interops3)
    add_sh_test(nc_test4 tst_formatx_hdf4)
    build_bin_test(tst_chunk_hdf4)
    add_sh_test(nc_test4 run_chunk_hdf4)
    add_bin_test(nc_test4 tst_h4_lendian)
    add_sh_test(nc_test4 tst_hdf4_read_var)
ENDIF()


IF(TEST_PARALLEL4)
  build_bin_test(tst_mpi_parallel)
  build_bin_test(tst_parallel)
  build_bin_test(tst_parallel3)
  build_bin_test(tst_parallel4)
  build_bin_test(tst_nc4perf)
  build_bin_test(tst_mode)
  build_bin_test(tst_simplerw_coll_r)
  add_sh_test(nc_test4 run_par_test)
ENDIF()
