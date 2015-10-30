ADD_DEFINITIONS(-DUSING_CMAKE)

SET( forbind_OLDTESTS
  testwt
  testwtd
  testrd
  testrd1
  testrdd
  testwt1
  testwt2
  testwtm
  testcp
  testcpnl
  testcpln
  testwt_nsid
  testrd_nsid
  test_nem
)

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/testall.in
    ${CMAKE_CURRENT_BINARY_DIR}/testall
    @ONLY
  )

  # NOTE: The "LINKER_LANGUAGE CXX" is used to workaround an issue with
  #       the trilinos cmake adding the pamgen and pamgen_extras libraries as
  #       a dependency.  These are C++ libraries and result in unresolved depedencies
  #       when built with shareable libraries on some platforms....

FOREACH ( testName ${forbind_OLDTESTS} )
  ADD_EXECUTABLE( f_${testName} ${testName}.f )
  TARGET_LINK_LIBRARIES( f_${testName} exoIIv2for )
  SET_TARGET_PROPERTIES(f_${testName} PROPERTIES LINKER_LANGUAGE Fortran)
ENDFOREACH ( testName )

ADD_CUSTOM_TARGET( f_check
  COMMAND /bin/csh ${CMAKE_CURRENT_BINARY_DIR}/testall
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
