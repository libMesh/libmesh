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
		FILE(GLOB IMG_FILES ${CMAKE_CURRENT_BINARY_DIR}/images/*.jpg)


		CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in.cmake
		${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
		ADD_CUSTOM_TARGET(doc ALL
			${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			COMMENT "Generating API Documentation with Doxygen" VERBATIM
		)
		FILE(COPY ${IMG_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/html/)
	ENDIF()
ENDIF()

SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} CMakeLists.txt Makefile.am netcdf.m4 DoxygenLayout.xml Doxyfile.in footer.html mainpage.doc tutorial.doc install.doc dispatch.doc guide.doc types.doc notes.doc cdl.doc architecture.doc internal.doc Doxyfile.in.cmake windows-binaries.md Building-with-CMake.md)

ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")