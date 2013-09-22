SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")
# Some extra tests
SET(NC4_TESTS tst_dims tst_dims2 tst_dims3 tst_files tst_files4 tst_vars tst_varms tst_unlim_vars tst_converts tst_converts2 tst_grps tst_grps2
tst_compounds tst_compounds2 tst_compounds3 tst_opaques tst_strings
tst_strings2 tst_interops tst_interops4 tst_interops6
tst_enums tst_coords tst_coords2 tst_coords3 tst_vars3 tst_vars4
tst_chunks tst_chunks2 tst_utf8 tst_fills tst_fills2 tst_fillbug
tst_xplatform tst_xplatform2 tst_h_atts2 tst_endian_fill tst_atts 
t_type cdm_sea_soundings tst_vl tst_atts1 tst_atts2
tst_vars2 tst_files5 tst_files6 tst_sync tst_h_strbug tst_h_refs)

IF(NOT MSVC)
	SET(NC4_TESTS ${NC4_TESTS} tst_interops5 tst_camrun)
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
ENDIF()

IF(TEST_PARALLEL)
  	add_bin_test(nc_test4 tst_mpi_parallel)	
	add_bin_test(nc_test4 tst_parallel)
	add_bin_test(nc_test4 tst_parallel3)
	add_bin_test(nc_test4 tst_parallel4)
	add_bin_test(nc_test4 tst_nc4perf)
ENDIF()	

IF(USE_PNETCDF)
	build_bin_test_no_prefix(tst_pnetcdf)
	build_bin_test_no_prefix(tst_parallel2)
	add_sh_test(nc_test4 run_pnetcdf_test)
	#add_sh_test(nc_test4 run_pnetcdf_test)
ENDIF()
