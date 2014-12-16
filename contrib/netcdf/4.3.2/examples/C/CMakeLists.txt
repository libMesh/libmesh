SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

SET(exam_C_tests simple_xy_wr simple_xy_rd sfc_pres_temp_wr sfc_pres_temp_rd pres_temp_4D_wr pres_temp_4D_rd)
SET(exam_C_tests_source "")
FOREACH(F ${exam_C_tests}) 
  set(exam_C_tests_source ${exam_C_test_source} ${F}.c)
ENDFOREACH()

FOREACH(F ${exam_C_tests})
  add_bin_test(C_tests ${F})
ENDFOREACH()


## Specify files to be distributed by 'make dist'
FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} CMakeLists.txt Makefile.am)
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")


SET(CLEANFILES sfc_pres_temp.nc simple_xy.nc pres_temp_4D.nc simple_nc4.nc simple_xy_nc4.nc)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")

