IF(BUILD_SHARED_LIBS AND WIN32)
  remove_definitions(-DDLL_EXPORT)
  remove_definitions(-DDLL_NETCDF)
ENDIF()

SET(ncgen_FILES generate.c main.c cdata.c bindata.c genchar.c cvt.c data.c debug.c escapes.c genc.c genbin.c generr.c genlib.c getfill.c odom.c semantics.c ncgeny.c dump.c util.c bytebuffer.c list.c genf77.c f77data.c genj.c jdata.c nc_iter.c ncgen.h)

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)


IF(USE_X_GETOPT)
  SET(ncgen_FILES ${ncgen_FILES} XGetopt.c)
ENDIF()

ADD_EXECUTABLE(ncgen ${ncgen_FILES})
TARGET_LINK_LIBRARIES(ncgen netcdf ${ALL_TLL_LIBS})

####
# We have to do a little tweaking
# to remove the Release/ and Debug/ directories
# in MSVC builds.  This is required to get
# test scripts to work.
####
IF(MSVC)
  SET_TARGET_PROPERTIES(ncgen PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(ncgen PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
    ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(ncgen PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
    ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()

INSTALL(TARGETS ncgen DESTINATION bin COMPONENT utilities)

# Copy the input for the flex file (ncgen.l)
FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.l ${CMAKE_CURRENT_SOURCE_DIR}/*.y)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# Stuff to build if tests are enabled.
IF(ENABLE_TESTS)
  FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cdl ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
  FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)

  ##
  # The following tests were moved to ncdump/, to ensure
  # that they are run in the proper order.
  ##
  #add_sh_test(ncgen run_tests)
  #IF(USE_NETCDF4)
  #  add_sh_test(ncgen run_nc4_tests)
  #ENDIF()
ENDIF()

SET(MAN_FILES ncgen.1)
IF(NOT MSVC)
  INSTALL(FILES ${MAN_FILES} DESTINATION "share/man/man1"
    COMPONENT documentation)
ENDIF()
SET(CLEANFILES c0.nc c0_64.nc c0_4.nc c0_4c.nc ref_camrun.c)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEANFILES}")

# These rules are used if someone wants to rebuild ncgeny.c or ncgenl.c
# Otherwise never invoked, but records how to do it.
# BTW: note that renaming is essential because otherwise
# autoconf will forcibly delete files of the name *.tab.*

ADD_CUSTOM_COMMAND(
  OUTPUT ncgentab.h
  COMMAND flex -Pncg -8 ncgen.l
  COMMAND rm -f ncgenl.c
  COMMAND mv lex.ncg.c ncgenl.c

  COMMAND bison -pncg -t -d ncgen.y
  COMMAND rm -f ncgeny.c ncgeny.h
  COMMAND mv ncgen.tab.c ncgeny.c
  COMMAND mv ncgen.tab.h ncgeny.h
  COMMAND mv ncgeny.h ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND mv ncgeny.c ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND mv ncgenl.c ${CMAKE_CURRENT_SOURCE_DIR}

  )
ADD_CUSTOM_TARGET(makeparser DEPENDS ncgentab.h)

## Specify files to be distributed by 'make dist'
FILE(GLOB CUR_EXTRA_DIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.sh)
SET(CUR_EXTRA_DIST ${CUR_EXTRA_DIST} ncgen.y ncgenl.c ncgen.l internals.html c0.cdl c0_4.cdl ref_camrun.cdl ncf199.cdl tst_gattenum.cdl CMakeLists.txt Makefile.am ncgen.1)
ADD_EXTRA_DIST("${CUR_EXTRA_DIST}")
