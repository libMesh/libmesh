SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

SET(ncgen_FILES generate.c main.c cdata.c bindata.c genchar.c cvt.c data.c debug.c escapes.c genc.c genbin.c generr.c genlib.c getfill.c odom.c offsets.c semantics.c ncgentab.c dump.c util.c bytebuffer.c list.c genf77.c f77data.c genj.c jdata.c nc_iter.c ConvertUTF.c  )

IF(USE_X_GETOPT)
  SET(ncgen_FILES ${ncgen_FILES} XGetopt.c)
ENDIF()

ADD_EXECUTABLE(ncgen ${ncgen_FILES})
TARGET_LINK_LIBRARIES(ncgen netcdf ${ALL_TLL_LIBS})

INSTALL(TARGETS ncgen DESTINATION bin COMPONENT utilities)

# Copy the input for the flex file (ncgen.l)
FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.l ${CMAKE_CURRENT_SOURCE_DIR}/*.y)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# Stuff to build if tests are enabled.
IF(ENABLE_TESTS)
  FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cdl ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
  FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE) 
  add_sh_test(ncgen run_tests)
  IF(USE_NETCDF4)
    add_sh_test(ncgen run_nc4_tests)
  ENDIF()
ENDIF()

SET(MAN_FILES ncgen.1)
IF(NOT MSVC)
  INSTALL(FILES ${MAN_FILES} DESTINATION "share/man/man1"
    COMPONENT documentation)
ENDIF()
SET(CLEANFILES c0.nc c0_64.nc c0_4.nc c0_4c.nc ref_camrun.c)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")

ADD_CUSTOM_COMMAND(
  OUTPUT ncgentab.h
  COMMAND flex -Pncg -B ncgen.l
  COMMAND rm -f ncgenyy.c
  COMMAND mv lex.ncg.c ncgenyy.c
  COMMAND bison -pncg -t -d ncgen.y
  COMMAND rm -f ncgentab.c ncgentab.h
  COMMAND mv ncgen.tab.c ncgentab.c
  COMMAND mv ncgen.tab.h ncgentab.h
  COMMAND mv ncgentab.h ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND mv ncgentab.c ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND mv ncgenyy.c ${CMAKE_CURRENT_SOURCE_DIR}

  )

ADD_CUSTOM_TARGET(makeparser DEPENDS ncgentab.h)

## Specify files to be distributed by 'make dist'
FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} ncgen.y ncgenyy.c ncgen.l internals.html c0.cdl c0_4.cdl ref_camrun.cdl ncf199.cdl CMakeLists.txt Makefile.am ncgen.1)
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")
