SET(exam_C_tests simple_xy_wr simple_xy_rd sfc_pres_temp_wr sfc_pres_temp_rd pres_temp_4D_wr pres_temp_4D_rd)
SET(exam_C_tests_source "")
FOREACH(F ${exam_C_tests})
  set(exam_C_tests_source ${exam_C_test_source} ${F}.c)
ENDFOREACH()

FOREACH(F ${exam_C_tests})
  add_bin_test(C_tests ${F})
ENDFOREACH()

SET_TESTS_PROPERTIES(C_tests_simple_xy_rd PROPERTIES DEPENDS C_tests_simple_xy_wr)
SET_TESTS_PROPERTIES(C_tests_sfc_pres_temp_rd PROPERTIES DEPENDS C_tests_sfc_pres_temp_wr)
SET_TESTS_PROPERTIES(C_tests_pres_temp_4D_rd PROPERTIES DEPENDS C_tests_pres_temp_4D_wr)

SET(CLEANFILES sfc_pres_temp.nc simple_xy.nc pres_temp_4D.nc simple_nc4.nc simple_xy_nc4.nc)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")
