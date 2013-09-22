SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

SET(H5TESTS tst_h_files tst_h_files2 tst_h_files4 tst_h_atts tst_h_atts3 tst_h_atts4 tst_h_vars tst_h_vars2 tst_h_vars3 tst_h_grps tst_h_compounds tst_h_compounds2 tst_h_wrt_cmp tst_h_rd_cmp tst_h_vl tst_h_opaques tst_h_strings tst_h_strings1 tst_h_strings2 tst_h_ints tst_h_dimscales tst_h_dimscales1 tst_h_dimscales2 tst_h_dimscales3 tst_h_enums tst_h_dimscales4)

FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h5 ${CMAKE_CURRENT_SOURCE_DIR}/*.nc)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

FOREACH(CTEST ${H5TESTS})
  ADD_EXECUTABLE(${CTEST} ${CTEST}.c)
  TARGET_LINK_LIBRARIES(${CTEST} netcdf)
  ADD_TEST(${CTEST} ${EXECUTABLE_OUTPUT_PATH}/${CTEST})
ENDFOREACH()






