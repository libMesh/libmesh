ADD_DEFINITIONS(-DUSING_CMAKE)

SET( cbind_TEST_SRCS
  CreateEdgeFace.c
  ReadEdgeFace.c
)

# Always report error/warning messages during testing:
# The ex_opts call precipitates the EXTRA_INCLUDE below.
SET( CMAKE_TESTDRIVER_BEFORE_TESTMAIN
  "  ex_opts( EX_VERBOSE );\n"
)
CREATE_TEST_SOURCELIST( ExoIICTests
  ExoIICTests.cxx
  ${cbind_TEST_SRCS}
  EXTRA_INCLUDE exodusII_test.h
  )

ADD_EXECUTABLE( ExoIICTests ${ExoIICTests} )
TARGET_LINK_LIBRARIES( ExoIICTests exoIIv2c )

SET( TestsToRun ${ExoIICTests} )
REMOVE( TestsToRun ExoIICTests.cxx )

FOREACH( test ${TestsToRun} )
  GET_FILENAME_COMPONENT(TName ${test} NAME_WE)
  ADD_TEST( ${TName} ${EXECUTABLE_OUTPUT_PATH}/ExoIICTests ${TName} )
ENDFOREACH( test )
ADD_TEST( CreateEdgeFaceWithConcats ${EXECUTABLE_OUTPUT_PATH}/ExoIICTests CreateEdgeFace -pcab -pcset -pvpax )
ADD_TEST( ReadEdgeFaceWithConcats ${EXECUTABLE_OUTPUT_PATH}/ExoIICTests ReadEdgeFace )

# ===============================================

SET( cbind_OLDTESTS
  testwt
  testwtd
  testrd
  testrd1
  testrdd
  testwt1
  testwt2
  testwtm
  testwt_ss
  testwt_nossnsdf
  testrd_ss
  testrdwt
  testcp
  testcp_ln
  testcp_nl
  testcpd
  testwt_clb
  testwt_nc
  testrd_nc
  testwt-zeroe
  testwt-zeron
  testwt-one-attrib
  testwt-partial
  testrd-nsided
  testwt-nsided
  testrd-nfaced
  testwt-nfaced
  create_mesh
  rd_wt_mesh
  test_nemesis
  testwt-long-name
  testrd-long-name)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/testall.in
  ${CMAKE_CURRENT_BINARY_DIR}/testall
  @ONLY
)

FOREACH ( testName ${cbind_OLDTESTS} )
  ADD_EXECUTABLE( ${testName} ${testName}.c )
  TARGET_LINK_LIBRARIES( ${testName} exoIIv2c )
ENDFOREACH ( testName )

#SET( EXECUTABLE_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR} )
ADD_CUSTOM_TARGET( check
  COMMAND /bin/csh ${CMAKE_CURRENT_BINARY_DIR}/testall
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${cbind_OLDTESTS}
)

