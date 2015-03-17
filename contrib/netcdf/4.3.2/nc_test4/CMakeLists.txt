SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")
# Some extra tests
SET(NC4_TESTS tst_dims tst_dims2 tst_dims3 tst_files tst_files4 tst_vars 
  tst_varms tst_unlim_vars tst_converts tst_converts2 tst_grps tst_grps2
  tst_compounds tst_compounds2 tst_compounds3 tst_opaques tst_strings
  tst_strings2 tst_interops tst_interops4 tst_interops6
  tst_enums tst_coords tst_coords2 tst_coords3 tst_vars3 tst_vars4
  tst_chunks tst_chunks2 tst_utf8 tst_fills tst_fills2 tst_fillbug
  tst_xplatform tst_xplatform2 tst_h_atts2 tst_endian_fill tst_atts 
  t_type cdm_sea_soundings tst_vl tst_atts1 tst_atts2
  tst_vars2 tst_files5 tst_files6 tst_sync tst_h_strbug tst_h_refs
  tst_h_scalar)

# Note, renamegroup needs to be compiled before run_grp_rename
build_bin_test(renamegroup)
add_sh_test(nc_test4 run_grp_rename)

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

IF(USE_HDF4_FILE_TESTS)
  add_sh_test(nc_test4 run_get_hdf4_files)
  add_bin_test(nc_test4 tst_interops3)
  add_sh_test(nc_test4 tst_formatx_hdf4)
ENDIF()

IF(TEST_PARALLEL)
  build_bin_test(tst_mpi_parallel)	
  build_bin_test(tst_parallel)
  build_bin_test(tst_parallel3)
  build_bin_test(tst_parallel4)
  build_bin_test(tst_nc4perf)
  add_sh_test(nc_test4 run_par_test)
ENDIF()	

## Specify files to be distributed by 'make dist'
FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} CMakeLists.txt Makefile.am)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} ref_chunks1.cdl ref_chunks2.cdl ref_tst_compounds.nc ref_tst_xplatform2_1.nc ref_tst_xplatform2_2.nc ref_tst_dims.nc ref_tst_interops4.nc ref_grp_rename.cdl ref_tst_nvars.nc)

ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")
