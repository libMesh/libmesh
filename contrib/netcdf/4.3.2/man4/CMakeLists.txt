SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")


#####
# Build doxygen documentation, if need be.
#####
IF(BUILD_DOCS)
  FIND_PACKAGE(Doxygen)
  IF(DOXYGEN_FOUND)
    FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.html ${CMAKE_CURRENT_SOURCE_DIR}/images ${CMAKE_CURRENT_SOURCE_DIR}/*.doc ${CMAKE_CURRENT_SOURCE_DIR}/*.xml ${CMAKE_CURRENT_SOURCE_DIR}/*.m4 ${CMAKE_CURRENT_SOURCE_DIR}/*.texi)
    FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
    # Copy the image files used by markdown documentation 
    # manually.
    FILE(GLOB IMG_FILES ${CMAKE_CURRENT_BINARY_DIR}/images/*.jpg ${CMAKE_CURRENT_BINARY_DIR}/images/*.png)
    
    # Set abs_top_srcdir to work with the autotools
    # doxyfile template.
    SET(abs_top_srcdir ${CMAKE_SOURCE_DIR})

    # Create general and guide doxyfiles.
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.guide.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.guide @ONLY)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.tutorial.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.tutorial @ONLY)

    # We need to loop through twice. First to generate tags.
    ADD_CUSTOM_TARGET(all_doc_tags ALL
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API Tags" VERBATIM
      )
    ADD_CUSTOM_TARGET(all_doc_tags_guide ALL
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.guide
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating Guide Tags" VERBATIM
      )
    ADD_CUSTOM_TARGET(all_doc_tags_tutorial ALL
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.tutorial
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating Tutorial Tags" VERBATIM
      )

    # Second, generate documentation that uses the tags.
    ADD_CUSTOM_TARGET(doc_all ALL
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API Documentation" VERBATIM
      )
    ADD_CUSTOM_TARGET(doc_guide ALL
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.guide
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating NetCDF Users Guide" VERBATIM
      )
    ADD_CUSTOM_TARGET(doc_tutorial ALL
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.tutorial
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating Tutorial Tags" VERBATIM
      )

    FILE(COPY ${IMG_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/html)		
  ENDIF()
ENDIF()

SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} CMakeLists.txt Makefile.am netcdf.m4 DoxygenLayout.xml Doxyfile.in Doxyfile.guide.in footer.html mainpage.dox tutorial.dox install.dox dispatch.dox guide.dox types.dox notes.dox cdl.dox architecture.dox internal.dox install-fortran.dox Doxyfile.in.cmake windows-binaries.md building-with-cmake.md)

ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")