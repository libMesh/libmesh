SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

IF(BUILD_SHARED_LIBS AND WIN32)
  remove_definitions(-DDLL_EXPORT)
  remove_definitions(-DDLL_NETCDF)
ENDIF()



#IF(EXTRA_TESTS)
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/ctest.c
  COMMAND cmake -E copy "${CMAKE_CURRENT_SOURCE_DIR}/ref_ctest.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest.c"
  )
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/ctest64.c
  COMMAND cmake -E copy "${CMAKE_CURRENT_SOURCE_DIR}/ref_ctest64.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/ctest64.c"
  )
#ENDIF()

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

IF(ENABLE_TESTS)
  ADD_EXECUTABLE(rewrite-scalar rewrite-scalar.c)
  ADD_EXECUTABLE(bom bom.c)
  TARGET_LINK_LIBRARIES(rewrite-scalar netcdf)
  TARGET_LINK_LIBRARIES(bom netcdf)
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
  add_sh_test(ncdump tst_nccopy3)
  add_sh_test(ncdump tst_charfill)
  add_sh_test(ncdump tst_iter)
  add_sh_test(ncdump tst_formatx3)
  add_sh_test(ncdump tst_bom)
  
  IF(EXTRA_TESTS)
    add_sh_test(ncdump run_back_comp_tests)
  ENDIF()
  
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
    add_sh_test(ncdump tst_formatx4)
    #add_bin_test(ncdump tst_special_atts)
    # Add this test by hand, as it is also called from a script.
    # Editing the script would break autotools compatibility.
    set(MANUAL_TESTS tst_special_atts tst_compress tst_chunking)
    FOREACH (MTEST ${MANUAL_TESTS})
      ADD_EXECUTABLE(${MTEST} ${MTEST}.c)
      TARGET_LINK_LIBRARIES(${MTEST} netcdf)
      IF(MSVC)
	SET_TARGET_PROPERTIES(${MTEST}
	  PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
	  )
      ENDIF()
      #ADD_TEST(${MTEST} ${EXECUTABLE_OUTPUT_PATH}/${MTEST})
    ENDFOREACH()
    
    IF(NOT WIN32)
      # tst_netcdf4 depends on the output from other scripts
      # that don't currently run. So for time being, leave
      # it disabled.

      #add_sh_test(ncdump tst_netcdf4)
      #add_bin_test(ncdump tst_h_rdc0)
    ENDIF()
    add_bin_test(ncdump tst_unicode)
    if(NOT WIN32)
      #add_bin_test(ncdump tst_fillbug)
      #add_sh_test(ncdump tst_fillbug)
    ENDIF()
    IF(NOT WIN32)	
      add_sh_test(ncdump tst_netcdf4_4)
    ENDIF()	
    
    IF(NOT WIN32)
      #add_sh_test(ncdump tst_nccopy4)
    ENDIF()
    add_sh_test(ncdump tst_grp_spec)
    add_sh_test(ncdump tst_mud)
    add_sh_test(ncdump_shell tst_h_scalar)
    
  ENDIF()
  
  add_sh_test(ncdump tst_ncgen4_classic)
  IF(USE_NETCDF4)
    add_sh_test(ncdump tst_ncgen4)
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

ADD_SUBDIRECTORY(cdl4)
ADD_SUBDIRECTORY(expected4)

## Specify files to be distributed by 'make dist'
FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} test0.cdl ref_ctest1_nc4.cdl ref_ctest1_nc4c.cdl ref_tst_solar_1.cdl ref_tst_solar_2.cdl ref_tst_small.cdl tst_ncml.cdl ref1.ncml ref_tst_group_data.cdl ref_tst_enum_data.cdl ref_tst_opaque_data.cdl ref_tst_string_data.cdl ref_tst_vlen_data.cdl ref_tst_comp.cdl ref_tst_unicode.cdl ref_tst_nans.cdl small.cdl small2.cdl ref_tst_utf8.cdl ref_tst_fillbug.cdl tst_calendars.cdl ref_times.cdl ref_tst_special_atts.cdl ref_tst_noncoord.cdl ref_tst_compounds2.nc ref_tst_compounds2.cdl ref_tst_compounds3.nc ref_tst_compounds3.cdl ref_tst_compounds4.nc ref_tst_compounds4.cdl ref_tst_group_data_v23.cdl tst_mslp.cdl ref_tst_format_att.cdl ref_tst_format_att_64.cdl ref_nc_test_netcdf4_4_0.nc ref_nc_test_netcdf4.cdl ref_tst_special_atts3.cdl tst_brecs.cdl ref_tst_grp_spec0.cdl ref_tst_grp_spec.cdl ref_tst_grp_spec0.cdl ref_tst_grp_spec.cdl ref_tst_charfill.cdl tst_charfill.cdl ref_tst_mud4.cdl ref_tst_mud4-bc.cdl ref_tst_ncf213.cdl CMakeLists.txt Makefile.am XGetopt.c nccopy.1 ncdump.1) 
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")



SET(CLEANFILES test0.nc test1.cdl test1.nc test2.cdl ctest1.cdl ctest.c ctest64.c ctest0.nc  ctest0_64.nc c1.cdl c1_4.cdl ctest1_64.cdl c0.nc c0_4.nc small.nc small2.nc c0tmp.nc c1.ncml utf8.cdl utf8_64.cdl utf8.nc utf8_64.nc tmp.cdl tst_vlen_data.nc tst_utf8.nc tst_special_atts.nc tst_unicode.nc tst_solar_2.nc tst_string_data.nc tst_calendars.nc tst_nans.nc tst_opaque_data.nc tst_solar_cmp.nc tst_enum_data.nc tst_solar_1.nc tst_mslp_64.nc tst_mslp.nc tst_comp2.nc tst_ncml.nc tst_fillbug.nc tst_group_data.nc tst_small.nc tst_comp.nc tst_unicode.cdl tst_group_data.cdl tst_compounds2.cdl tst_comp.cdl tst_enum_data.cdl tst_small.cdl tst_times.cdl tst_solar_2.cdl tst_string_data.cdl tst_fillbug.cdl tst_opaque_data.cdl tst_compounds4.cdl tst_utf8.cdl tst_compounds3.cdl tst_special_atts.cdl tst_nans.cdl tst_format_att_64.cdl tst_vlen_data.cdl tst_solar_1.cdl tst_format_att.cdl tst_inflated.nc tmp_subset.cdl tst_inflated4.nc tst_deflated.nc tst_chunking.nc tmp*.nc tst_charfill.nc tmp_tst_charfill.cdl iter.* tst_nc_test_netcdf4_4_0.cdl tst_mud4.nc tst_mud4.cdl tst_mud4-bc.cdl)

SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")

IF(NOT MSVC)	
  INSTALL(FILES ${MAN_FILES} DESTINATION "share/man/man1" COMPONENT documentation)
ENDIF()

