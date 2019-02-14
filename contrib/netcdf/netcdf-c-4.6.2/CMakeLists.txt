## This is a CMake file, part of Unidata's netCDF package.
# Copyright 2012-2018, see the COPYRIGHT file for more information.
#

##################################
# Set Project Properties
##################################

#Minimum required CMake Version
cmake_minimum_required(VERSION 3.6.1)

#Project Name
project(netCDF C)
set(PACKAGE "netCDF" CACHE STRING "")

#####
# Version Info:
#
# Release Version
# Library Version
# SO Version
#
# SO Version is computed from library version. See:
# http://www.gnu.org/software/libtool/manual/libtool.html#Libtool-versioning
#####

SET(NC_VERSION_MAJOR 4)
SET(NC_VERSION_MINOR 6)
SET(NC_VERSION_PATCH 2)
SET(NC_VERSION_NOTE "")
SET(netCDF_VERSION ${NC_VERSION_MAJOR}.${NC_VERSION_MINOR}.${NC_VERSION_PATCH}${NC_VERSION_NOTE})
SET(VERSION ${netCDF_VERSION})
SET(NC_VERSION ${netCDF_VERSION})
SET(netCDF_LIB_VERSION 13)
SET(netCDF_SO_VERSION 13)
SET(PACKAGE_VERSION ${VERSION})

# Get system configuration, Use it to determine osname, os release, cpu. These
# will be used when committing to CDash.
find_program(UNAME NAMES uname)
IF(UNAME)
  macro(getuname name flag)
    exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
  endmacro(getuname)
  getuname(osname -s)
  getuname(osrel  -r)
  getuname(cpu    -m)
  set(TMP_BUILDNAME "${osname}-${osrel}-${cpu}")
ENDIF()

###
# Allow for some customization of the buildname.
# This will make it easier to identify different builds,
# based on values passed from command line/shell scripts.
#
# For ctest scripts, we can use CTEST_BUILD_NAME.
###

SET(BUILDNAME_PREFIX "" CACHE STRING "")
SET(BUILDNAME_SUFFIX "" CACHE STRING "")

IF(BUILDNAME_PREFIX)
  SET(TMP_BUILDNAME "${BUILDNAME_PREFIX}-${TMP_BUILDNAME}")
ENDIF()

IF(BUILDNAME_SUFFIX)
  SET(TMP_BUILDNAME "${TMP_BUILDNAME}-${BUILDNAME_SUFFIX}")
ENDIF()

IF(NOT BUILDNAME)
  SET(BUILDNAME "${TMP_BUILDNAME}" CACHE STRING "Build name variable for CDash")
ENDIF()
###
# End BUILDNAME customization.
###

# For CMAKE_INSTALL_LIBDIR
INCLUDE(GNUInstallDirs)

IF(MSVC)
  SET(GLOBAL PROPERTY USE_FOLDERS ON)
ENDIF()

#Add custom CMake Module
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/"
  CACHE INTERNAL "Location of our custom CMake modules.")

# auto-configure style checks, other CMake modules.
INCLUDE(CheckLibraryExists)
INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckTypeSize)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckCSourceCompiles)
INCLUDE(TestBigEndian)
INCLUDE(CheckSymbolExists)
INCLUDE(GetPrerequisites)

INCLUDE(CheckCCompilerFlag)
FIND_PACKAGE(PkgConfig QUIET)

# A check to see if the system is big endian
TEST_BIG_ENDIAN(BIGENDIAN)
IF(${BIGENDIAN})
  SET(WORDS_BIGENDIAN "1")
ENDIF(${BIGENDIAN})

# A macro to check if a C linker supports a particular flag.
MACRO(CHECK_C_LINKER_FLAG M_FLAG M_RESULT)
  SET(T_REQ_FLAG "${CMAKE_REQUIRED_FLAGS}")
  SET(CMAKE_REQUIRED_FLAGS "${M_FLAG}")
  CHECK_C_SOURCE_COMPILES("int main() {return 0;}" ${M_RESULT})
  SET(CMAKE_REQUIRED_FLAGS "${T_REQ_FLAG}")
ENDMACRO()

# Enable 'dist and distcheck'.
# File adapted from http://ensc.de/cmake/FindMakeDist.cmake
FIND_PACKAGE(MakeDist)
# End 'enable dist and distcheck'

# Set the build type.
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Choose the type of build, options are: None, Debug, Release."
    FORCE)
ENDIF()

# Set build type uppercase
STRING(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

# Determine the configure date.

IF(DEFINED ENV{SOURCE_DATE_EPOCH})
	EXECUTE_PROCESS(
	  COMMAND "date" "-u" "-d" "@$ENV{SOURCE_DATE_EPOCH}"
	  OUTPUT_VARIABLE CONFIG_DATE
	  )
ELSE()
	EXECUTE_PROCESS(
	  COMMAND date
	  OUTPUT_VARIABLE CONFIG_DATE
	  )
ENDIF()
IF(CONFIG_DATE)
	string(STRIP ${CONFIG_DATE} CONFIG_DATE)
ENDIF()
##
# Allow for extra dependencies.
##

SET(EXTRA_DEPS "")

################################
# End Project Properties
################################


################################
# Set CTest Properties
################################



ENABLE_TESTING()
INCLUDE(CTest)

# Copy the CTest customization file into binary directory, as required.
FILE(COPY ${CMAKE_SOURCE_DIR}/CTestCustom.cmake DESTINATION ${CMAKE_BINARY_DIR})

# Set Memory test program for non-MSVC based builds.
# Assume valgrind for now.
IF(NOT MSVC)
  SET(CTEST_MEMORYCHECK_COMMAND valgrind CACHE STRING "")
ENDIF()

# Set variable to define the build type.
INCLUDE(GenerateExportHeader)

################################
# End CTest Properties
################################


################################
# Compiler and Linker Configuration
################################

##
# Default building shared libraries.
# BUILD_SHARED_LIBS is provided by/used by
# CMake directly.
##
OPTION(BUILD_SHARED_LIBS "Configure netCDF as a shared library." ON)
IF(BUILD_SHARED_LIBS)
  SET(CMAKE_POSITION_INDEPENDENT_CODE ON)
ENDIF()

OPTION(NC_FIND_SHARED_LIBS "Find dynamically-built versions of dependent libraries" ${BUILD_SHARED_LIBS})

##
# We've had a request to allow for non-versioned shared libraries.
# This seems reasonable enough to accommodate.  See
# https://github.com/Unidata/netcdf-c/issues/228 for more info.
##
OPTION(ENABLE_SHARED_LIBRARY_VERSION "Encode the library SO version in the file name of the generated library file." ON)

# Set some default linux gcc & apple compiler options for
# debug builds.
IF(CMAKE_COMPILER_IS_GNUCC OR APPLE)
  OPTION(ENABLE_COVERAGE_TESTS "Enable compiler flags needed to perform coverage tests." OFF)
  OPTION(ENABLE_CONVERSION_WARNINGS "Enable warnings for implicit conversion from 64 to 32-bit datatypes." ON)
  OPTION(ENABLE_LARGE_FILE_TESTS "Enable large file tests." OFF)

  # Debugging flags
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")

  # Check to see if -Wl,--no-undefined is supported.
  CHECK_C_LINKER_FLAG("-Wl,--no-undefined" LIBTOOL_HAS_NO_UNDEFINED)

  IF(LIBTOOL_HAS_NO_UNDEFINED)
    SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -Wl,--no-undefined")
  ENDIF()
  SET(CMAKE_REQUIRED_FLAGS "${TMP_CMAKE_REQUIRED_FLAGS}")

  # Coverage tests need to have optimization turned off.
  IF(ENABLE_COVERAGE_TESTS)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
    MESSAGE(STATUS "Coverage Tests: On.")
  ENDIF()

    # Warnings for 64-to-32 bit conversions.
  IF(ENABLE_CONVERSION_WARNINGS)
    CHECK_C_COMPILER_FLAG(-Wconversion CC_HAS_WCONVERSION)
    CHECK_C_COMPILER_FLAG(-Wshorten-64-to-32 CC_HAS_SHORTEN_64_32)

    IF(CC_HAS_SHORTEN_64_32)
      SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wshorten-64-to-32")
    ENDIF()
    IF(CC_HAS_WCONVERSION)
      SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wconversion")
    ENDIF()

  ENDIF(ENABLE_CONVERSION_WARNINGS)

ENDIF(CMAKE_COMPILER_IS_GNUCC OR APPLE)

# End default linux gcc & apple compiler options.

ADD_DEFINITIONS()

# Suppress CRT Warnings.
# Only necessary for Windows
IF(MSVC)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF()

#####
# System inspection checks
#####
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/oc2)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/libsrc)
SET(CMAKE_REQUIRED_INCLUDES ${CMAKE_SOURCE_DIR}/libsrc)

################################
# End Compiler Configuration
################################



##
# Configuration for post-install RPath
# Adapted from http://www.cmake.org/Wiki/CMake_RPATH_handling
##
IF(NOT MSVC AND BUILD_SHARED_LIBS)
  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
  endif(APPLE)

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing,
  # but only if it's not a system directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  ENDIF("${isSystemDir}" STREQUAL "-1")

ENDIF()

##
# End configuration for post-install RPath
##

################################
# Option checks
################################

# HDF5 cache variables.
SET(DEFAULT_CHUNK_SIZE 4194304 CACHE STRING "Default Chunk Cache Size.")
SET(DEFAULT_CHUNKS_IN_CACHE 10 CACHE STRING "Default number of chunks in cache.")
SET(CHUNK_CACHE_SIZE 4194304 CACHE STRING "Default Chunk Cache Size.")
SET(CHUNK_CACHE_NELEMS 1009 CACHE STRING "Default maximum number of elements in cache.")
SET(CHUNK_CACHE_PREEMPTION 0.75 CACHE STRING "Default file chunk cache preemption policy for HDf5 files(a number between 0 and 1, inclusive.")
SET(MAX_DEFAULT_CACHE_SIZE 67108864 CACHE STRING "Default maximum cache size.")
SET(NETCDF_LIB_NAME "" CACHE STRING "Default name of the netcdf library.")
SET(TEMP_LARGE "." CACHE STRING "Where to put large temp files if large file tests are run.")
SET(NCPROPERTIES_EXTRA "" CACHE STRNG "Specify extra pairs for _NCProperties.")
SET(MPIEXEC "mpiexec" CACHE STRING "Command to run MPI programs if parallel tests are run.")

IF(NOT NETCDF_LIB_NAME STREQUAL "")
  SET(MOD_NETCDF_NAME ON)
ENDIF()

# Set the appropriate compiler/architecture for universal OSX binaries.
IF(${CMAKE_SYSTEM_NAME} EQUAL "Darwin")
  SET(CMAKE_OSX_ARCHITECTURES i386;x86_64)
ENDIF(${CMAKE_SYSTEM_NAME} EQUAL "Darwin")

# Macro for replacing '/MD' with '/MT'.
# Used only on Windows, /MD tells VS to use the shared
# CRT libs, MT tells VS to use the static CRT libs.
#
# Taken From:
#   http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F
#
MACRO(specify_static_crt_flag)
  SET(vars
    CMAKE_C_FLAGS
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_MINSIZEREL
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELWITHDEBINFO)

  FOREACH(flag_var ${vars})
    IF(${flag_var} MATCHES "/MD")
      STRING(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    ENDIF()
  ENDFOREACH()

  FOREACH(flag_var ${vars})
    MESSAGE(STATUS " '${flag_var}': ${${flag_var}}")
  ENDFOREACH()
  MESSAGE(STATUS "")
ENDMACRO()

# Option to use Static Runtimes in MSVC
IF(MSVC)
  OPTION(NC_USE_STATIC_CRT "Use static CRT Libraries ('\\MT')." OFF)
  IF(NC_USE_STATIC_CRT)
    SET(USE_STATIC_CRT ON)
    specify_static_crt_flag()
  ENDIF()
ENDIF()

# Option to build netCDF Version 2
OPTION (ENABLE_V2_API "Build netCDF Version 2." ON)
SET(BUILD_V2 ${ENABLE_V2_API})
IF(NOT ENABLE_V2_API)
  SET(NO_NETCDF_2 ON)
ELSE(NOT ENABLE_V2_API)
  SET(USE_NETCDF_2 TRUE)
ENDIF(NOT ENABLE_V2_API)

# Option to build utilities
OPTION(BUILD_UTILITIES "Build ncgen, ncgen3, ncdump." ON)

# Option to use MMAP
OPTION(ENABLE_MMAP "Use MMAP." ON)

# Option to use examples.
OPTION(ENABLE_EXAMPLES "Build Examples" ON)

# Option to automatically build netcdf-fortran.
IF(NOT MSVC)
  OPTION(ENABLE_REMOTE_FORTRAN_BOOTSTRAP "Download and build netcdf-fortran automatically (EXPERIMENTAL)." OFF)
  IF(ENABLE_REMOTE_FORTRAN_BOOTSTRAP)
    SET(BUILD_FORTRAN ON)
  ENDIF()
  IF(BUILD_FORTRAN)
    CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/postinstall.sh.in"
      "${CMAKE_BINARY_DIR}/postinstall.sh"
      @ONLY)

    ADD_CUSTOM_TARGET(build-netcdf-fortran
      COMMAND sh -c "${CMAKE_BINARY_DIR}/postinstall.sh -t cmake -a build"
      DEPENDS netcdf
      )

    ADD_CUSTOM_TARGET(install-netcdf-fortran
      COMMAND sh -c "${CMAKE_BINARY_DIR}/postinstall.sh -t cmake -a install"
      DEPENDS build-netcdf-fortran
      )

  ENDIF(BUILD_FORTRAN)
ENDIF()


###
# Allow the user to specify libraries
# to link against, similar to automakes 'LIBS' variable.
###
SET(NC_EXTRA_DEPS "" CACHE STRING "Additional libraries to link against.")
IF(NC_EXTRA_DEPS)
  STRING(REPLACE " " ";" DEPS_LIST ${NC_EXTRA_DEPS})
  FOREACH(_DEP ${DEPS_LIST})
    STRING(REGEX REPLACE "^-l" "" _LIB ${_DEP})
    FIND_LIBRARY("${_LIB}_DEP" NAMES "${_LIB}" "lib${_LIB}")
    MESSAGE(${${_LIB}_DEP})
    IF("${${_LIB}_DEP}" STREQUAL "${_LIB}_DEP-NOTFOUND")
      MESSAGE(FATAL_ERROR "Error finding ${_LIB}.")
    ELSE()
      MESSAGE(STATUS "Found ${_LIB}: ${${_LIB}_DEP}")
    ENDIF()
    SET(EXTRA_DEPS ${EXTRA_DEPS} "${${_LIB}_DEP}")
  ENDFOREACH()
  MESSAGE("Extra deps: ${EXTRA_DEPS}")
  LIST(REMOVE_DUPLICATES EXTRA_DEPS)
  SET(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${EXTRA_DEPS})
ENDIF()
###
# End user-specified dependent libraries.
###

# Option to use HDF4
OPTION(ENABLE_HDF4 "Build netCDF-4 with HDF4 read capability(HDF4, HDF5 and Zlib required)." OFF)
IF(ENABLE_HDF4)
  SET(USE_HDF4 ON)
  # Check for include files, libraries.

  FIND_PATH(MFHDF_H_INCLUDE_DIR mfhdf.h)
  IF(NOT MFHDF_H_INCLUDE_DIR)
    MESSAGE(FATAL_ERROR "HDF4 Support specified, cannot find file mfhdf.h")
  ELSE()
    INCLUDE_DIRECTORIES(${MFHDF_H_INCLUDE_DIR})
  ENDIF()

  FIND_LIBRARY(HDF4_DF_LIB NAMES df libdf hdf)
  IF(NOT HDF4_DF_LIB)
    MESSAGE(FATAL_ERROR "Can't find or link to the hdf4 df library.")
  ENDIF()

  FIND_LIBRARY(HDF4_MFHDF_LIB NAMES mfhdf libmfhdf)
  IF(NOT HDF4_MFHDF_LIB)
    MESSAGE(FATAL_ERROR "Can't find or link to the hdf4 mfhdf library.")
  ENDIF()

  SET(HAVE_LIBMFHDF TRUE)

  SET(HDF4_LIBRARIES ${HDF4_DF_LIB} ${HDF4_MFHDF_LIB})
  # End include files, libraries.
  MESSAGE(STATUS "HDF4 libraries: ${HDF4_DF_LIB}, ${HDF4_MFHDF_LIB}")

  MESSAGE(STATUS "Seeking HDF4 jpeg dependency.")

  # Look for the jpeglib.h header file.
  FIND_PATH(JPEGLIB_H_INCLUDE_DIR jpeglib.h)
  IF(NOT JPEGLIB_H_INCLUDE_DIR)
    MESSAGE(FATAL_ERROR "HDF4 Support enabled but cannot find jpeglib.h")
  ELSE()
    SET(HAVE_JPEGLIB_H ON CACHE BOOL "")
    SET(HAVE_LIBJPEG TRUE)
    INCLUDE_DIRECTORIES(${JPEGLIB_H_INCLUDE_DIR})
  ENDIF()

  FIND_LIBRARY(JPEG_LIB NAMES jpeg libjpeg)
  IF(NOT JPEG_LIB)
    MESSAGE(FATAL_ERROR "HDF4 Support enabled but cannot find libjpeg")
  ENDIF()
  SET(HDF4_LIBRARIES ${JPEG_LIB} ${HDF4_LIBRARIES})
  MESSAGE(STATUS "Found JPEG libraries: ${JPEG_LIB}")

  # Option to enable HDF4 file tests.
  OPTION(ENABLE_HDF4_FILE_TESTS "Run HDF4 file tests.  This fetches sample HDF4 files from the Unidata ftp site to test with (requires curl)." ON)
  IF(ENABLE_HDF4_FILE_TESTS)
    FIND_PROGRAM(PROG_CURL NAMES curl)
    IF(PROG_CURL)
      SET(USE_HDF4_FILE_TESTS ON)
    ELSE()
      MESSAGE(STATUS "Unable to locate 'curl'.  Disabling hdf4 file tests.")
      SET(USE_HDF4_FILE_TESTS OFF)
    ENDIF()
  ENDIF()
ENDIF()

# Option to Build DLL
IF(WIN32)
  OPTION(ENABLE_DLL "Build a Windows DLL." ${BUILD_SHARED_LIBS})
  IF(ENABLE_DLL)
    SET(BUILD_DLL ON CACHE BOOL "")
    ADD_DEFINITIONS(-DDLL_NETCDF)
    ADD_DEFINITIONS(-DDLL_EXPORT)
    ADD_DEFINITIONS(-DUTF8PROC_DLLEXPORT)
  ENDIF()
ENDIF()
# Did the user specify a default minimum blocksize for posixio?
SET(NCIO_MINBLOCKSIZE 256 CACHE STRING "Minimum I/O Blocksize for netCDF classic and 64-bit offset format files.")

# Build netCDF4
OPTION(ENABLE_NETCDF_4 "Enable netCDF-4" ON)
IF(ENABLE_NETCDF_4)
  SET(USE_NETCDF4 ON CACHE BOOL "")
  SET(ENABLE_NETCDF_4 ON CACHE BOOL "")
  SET(ENABLE_NETCDF4 ON CACHE BOOL "")
ELSE()
  SET(USE_HDF4_FILE_TESTS OFF)
  SET(USE_HDF4 OFF)
  SET(ENABLE_HDF4_FILE_TESTS OFF)
  SET(ENABLE_HDF4 OFF)
ENDIF()

# Option Logging, only valid for netcdf4.
OPTION(ENABLE_LOGGING "Enable Logging." OFF)
IF(NOT ENABLE_NETCDF_4)
SET(ENABLE_LOGGING OFF)
ENDIF()
IF(ENABLE_LOGGING)
  ADD_DEFINITIONS(-DLOGGING)
  ADD_DEFINITIONS(-DENABLE_SET_LOG_LEVEL)
  SET(LOGGING ON)
  SET(ENABLE_SET_LOG_LEVEL ON)
ENDIF()
OPTION(ENABLE_SET_LOG_LEVEL_FUNC "Enable definition of nc_set_log_level()." ON)
IF(ENABLE_NETCDF_4 AND NOT ENABLE_LOGGING AND ENABLE_SET_LOG_LEVEL_FUNC)
  ADD_DEFINITIONS(-DENABLE_SET_LOG_LEVEL)
  SET(ENABLE_SET_LOG_LEVEL ON)
ENDIF()

# Option to allow for strict null file padding.
# See https://github.com/Unidata/netcdf-c/issues/657 for more information
OPTION(ENABLE_STRICT_NULL_BYTE_HEADER_PADDING "Enable strict null byte header padding." OFF)

IF(ENABLE_STRICT_NULL_BYTE_HEADER_PADDING)
  SET(USE_STRICT_NULL_BYTE_HEADER_PADDING ON CACHE BOOL "")
ENDIF(ENABLE_STRICT_NULL_BYTE_HEADER_PADDING)

# Option for building RPC
OPTION(ENABLE_RPC "Enable RPC Client and Server." OFF)
IF(ENABLE_RPC)
  SET(BUILD_RPC ON CACHE BOOL "")
ENDIF()

##
# Option to Enable HDF5
#
# The HDF5 cmake variables differ between platform (linux/osx and Windows),
# as well as between HDF5 versions.  As a result, this section is a bit convoluted.
#
# Note that the behavior seems much more stable across HDF5 versions under linux,
# so we do not have to do as much version-based tweaking.
#
# At the end of it, we should have the following defined:
#
# * HDF5_C_LIBRARY
# * HDF5_HL_LIBRARY
# * HDF5_LIBRARIES
# * HDF5_INCLUDE_DIR
# *
##
OPTION(USE_HDF5 "Use HDF5." ${ENABLE_NETCDF_4})
IF(USE_HDF5 OR ENABLE_NETCDF_4)
  SET(USE_HDF5 ON)
  SET(USE_NETCDF4 ON)
  ##
  # Accommodate developers who have hdf5 libraries and
  # headers on their system, but do not have a the hdf
  # .cmake files.  If this is the case, they should
  # specify HDF5_HL_LIB, HDF5_LIB, HDF5_INCLUDE_DIR manually.
  ##
  IF(HDF5_C_LIBRARY AND HDF5_HL_LIBRARY AND HDF5_INCLUDE_DIR)
    SET(HDF5_LIBRARIES ${HDF5_C_LIBRARY} ${HDF5_HL_LIBRARY})
    SET(HDF5_C_LIBRARIES ${HDF5_C_LIBRARY})
    SET(HDF5_C_LIBRARY_hdf5 ${HDF5_C_LIBRARY})
    SET(HDF5_HL_LIBRARIES ${HDF5_HL_LIBRARY})
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
    MESSAGE(STATUS "Using HDF5 C Library: ${HDF5_C_LIBRARY}")
    MESSAGE(STATUS "Using HDF5 HL LIbrary: ${HDF5_HL_LIBRARY}")
  ELSE(HDF5_C_LIBRARY AND HDF5_HL_LIBRARY AND HDF5_INCLUDE_DIR) # We are seeking out HDF5 with Find Package.
    ###
    # For now we assume that if we are building netcdf
    # as a shared library, we will use hdf5 as a shared
    # library. If we are building netcdf statically,
    # we will use a static library.  This can be toggled
    # by explicitly modifying NC_FIND_SHARED_LIBS.
    ##
    IF(NC_FIND_SHARED_LIBS)
      SET(NC_HDF5_LINK_TYPE "shared")
      SET(NC_HDF5_LINK_TYPE_UPPER "SHARED")
      ADD_DEFINITIONS(-DH5_BUILT_AS_DYNAMIC_LIB)
    ELSE(NC_FIND_SHARED_LIBS)
      SET(NC_HDF5_LINK_TYPE "static")
      SET(NC_HDF5_LINK_TYPE_UPPER "STATIC")
      ADD_DEFINITIONS(-DH5_BUILT_AS_STATIC_LIB)
    ENDIF(NC_FIND_SHARED_LIBS)

    #####
    # First, find the C and HL libraries.
    #
    # This has been updated to reflect what is in the hdf5
    # examples, even though the previous version of what we
    # had worked.
    #####
    IF(MSVC)
      SET(SEARCH_PACKAGE_NAME ${HDF5_PACKAGE_NAME})
      FIND_PACKAGE(HDF5 NAMES ${SEARCH_PACKAGE_NAME} COMPONENTS C HL NO_MODULES REQUIRED ${NC_HDF5_LINK_TYPE})
    ELSE(MSVC)
      FIND_PACKAGE(HDF5 COMPONENTS C HL REQUIRED)
    ENDIF(MSVC)

    ##
    # Next, check the HDF5 version. This will inform which
    # HDF5 variables we need to munge.
    ##

    ##
    # Assert HDF5 version meets minimum required version.
    ##
    SET(HDF5_VERSION_REQUIRED 1.8.10)

    IF(HDF5_VERSION_STRING AND NOT HDF5_VERSION)
      SET(HDF5_VERSION ${HDF5_VERSION_STRING})
    ENDIF()

    IF("${HDF5_VERSION}" STREQUAL "")
      MESSAGE(STATUS "Unable to determine hdf5 version.  NetCDF requires at least version ${HDF5_VERSION_REQUIRED}")
    ELSE()
      IF(${HDF5_VERSION} VERSION_LESS ${HDF5_VERSION_REQUIRED})
        MESSAGE(FATAL_ERROR
	      "netCDF requires at least HDF5 ${HDF5_VERSION_REQUIRED}. Found ${HDF5_VERSION}.")
      ELSE()
        MESSAGE(STATUS "Found HDF5 libraries version ${HDF5_VERSION}")
      ENDIF()
    ENDIF()

    ##
    # Include the HDF5 include directory.
    ##
    IF(HDF5_INCLUDE_DIRS AND NOT HDF5_INCLUDE_DIR)
      SET(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})
    ENDIF()
    MESSAGE(STATUS "Using HDF5 include dir: ${HDF5_INCLUDE_DIR}")
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

    ###
    # This is the block where we figure out what the appropriate
    # variables are, and we ensure that we end up with
    # HDF5_C_LIBRARY, HDF5_HL_LIBRARY and HDF5_LIBRARIES.
    ###
    IF(MSVC)
      ##
      # HDF5 1.8.15 defined HDF5_LIBRARIES.
      ##
      IF(${HDF5_VERSION} VERSION_LESS "1.8.16")
        SET(HDF5_C_LIBRARY hdf5)
        SET(HDF5_C_LIBRARY_hdf5 hdf5)
      ENDIF(${HDF5_VERSION} VERSION_LESS "1.8.16")

      IF(${HDF5_VERSION} VERSION_GREATER "1.8.15")
        IF(NOT HDF5_LIBRARIES AND HDF5_C_${NC_HDF5_LINK_TYPE_UPPER}_LIBRARY AND HDF5_HL_${NC_HDF5_LINK_TYPE_UPPER}_LIBRARY)
          SET(HDF5_C_LIBRARY ${HDF5_C_${NC_HDF5_LINK_TYPE_UPPER}_LIBRARY})
          SET(HDF5_C_LIBRARY_hdf5 ${HDF5_C_${NC_HDF5_LINK_TYPE_UPPER}_LIBRARY})
          SET(HDF5_HL_LIBRARY ${HDF5_HL_${NC_HDF5_LINK_TYPE_UPPER}_LIBRARY})

      	  SET(HDF5_LIBRARIES ${HDF5_C_${NC_HDF5_LINK_TYPE_UPPER}_LIBRARY} ${HDF5_HL_${NC_HDF5_LINK_TYPE_UPPER}_LIBRARY})
        ENDIF()
      ENDIF(${HDF5_VERSION} VERSION_GREATER "1.8.15")

    ELSE(MSVC)

      # Depending on the install, either HDF5_hdf_library or
      # HDF5_C_LIBRARIES may be defined.  We must check for either.
      IF(HDF5_C_LIBRARIES AND NOT HDF5_hdf5_LIBRARY)
        SET(HDF5_hdf5_LIBRARY ${HDF5_C_LIBRARIES})
      ENDIF()

    ENDIF(MSVC)
    IF(NOT HDF5_C_LIBRARY)
      SET(HDF5_C_LIBRARY hdf5)
    ENDIF()
  ENDIF(HDF5_C_LIBRARY AND HDF5_HL_LIBRARY AND HDF5_INCLUDE_DIR)

  FIND_PACKAGE(Threads)

  # There is a missing case in the above code so default it
  IF(NOT HDF5_C_LIBRARY_HDF5 OR "${HDF5_C_LIBRARY_hdf5}" STREQUAL "" )
    SET(HDF5_C_LIBRARY_hdf5 "${HDF5_C_LIBRARY}")
  ENDIF()

  ###
  # The following options are not used in Windows currently.
  ###
  IF(NOT MSVC)
    # Find out if HDF5 was built with parallel support.
    # Do that by checking for the targets H5Pget_fapl_mpiposx and
    # H5Pget_fapl_mpio in ${HDF5_LIB}.

    # H5Pset_fapl_mpiposix and H5Pget_fapl_mpiposix have been removed since HDF5 1.8.12.
    # Use H5Pset_fapl_mpio and H5Pget_fapl_mpio, instead.
    # CHECK_LIBRARY_EXISTS(${HDF5_C_LIBRARY_hdf5} H5Pget_fapl_mpiposix "" HDF5_IS_PARALLEL_MPIPOSIX)

    CHECK_LIBRARY_EXISTS(${HDF5_C_LIBRARY_hdf5} H5Pget_fapl_mpio "" HDF5_IS_PARALLEL_MPIO)
    IF(HDF5_IS_PARALLEL_MPIO)
      SET(HDF5_PARALLEL ON)
    ELSE()
      SET(HDF5_PARALLEL OFF)
    ENDIF()

    #Check to see if HDF5 library has collective metadata APIs, (HDF5 >= 1.10.0)
    CHECK_LIBRARY_EXISTS(${HDF5_C_LIBRARY_hdf5} H5Pset_all_coll_metadata_ops "" HDF5_HAS_COLL_METADATA_OPS)
  ENDIF(NOT MSVC)

  #Check to see if H5Z_SZIP exists in HDF5_Libraries. If so, we must use szip.
  CHECK_LIBRARY_EXISTS(${HDF5_C_LIBRARY_hdf5} H5Z_SZIP "" USE_SZIP)
  IF(USE_SZIP)
    FIND_LIBRARY(SZIP NAMES szip sz)
    IF(SZIP)
      SET(HAVE_H5Z_SZIP 1)
      SET(SZIP_LIBRARY ${SZIP})
      SET(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${SZIP})
    ELSE()
      MESSAGE(FATAL_ERROR "HDF5 Requires SZIP, but cannot find libszip or libsz.")
    ENDIF()
  ENDIF()

  CHECK_LIBRARY_EXISTS(${HDF5_C_LIBRARY_hdf5} H5Pset_libver_bounds "" HDF5_HAS_LIBVER_BOUNDS)
  CHECK_LIBRARY_EXISTS(${HDF5_C_LIBRARY_hdf5} H5free_memory "" HDF5_HAS_H5FREE)

  IF(HDF5_PARALLEL)
	SET(HDF5_CC h5pcc)
  ELSE()
	SET(HDF5_CC h5cc)
  ENDIF()

  SET(H5_USE_16_API 1)
  OPTION(NC_ENABLE_HDF_16_API "Enable HDF5 1.6.x Compatibility(Required)" ON)
  IF(NOT NC_ENABLE_HDF_16_API)
    SET(H5_USE_16_API 0)
  ENDIF()

  FIND_PATH(HAVE_HDF5_H hdf5.h)
  IF(NOT HAVE_HDF5_H)
    MESSAGE(FATAL_ERROR "Compiling a test with hdf5 failed. Either hdf5.h cannot be found, or the log messages should be checked for another reason.")
  ELSE(NOT HAVE_HDF5_H)
    INCLUDE_DIRECTORIES(${HAVE_HDF5_H})
  ENDIF(NOT HAVE_HDF5_H)

  #option to include HDF5 High Level header file (hdf5_hl.h) in case we are not doing a make install
  INCLUDE_DIRECTORIES(${HDF5_HL_INCLUDE_DIR})


ENDIF(USE_HDF5 OR ENABLE_NETCDF_4)

# Option to Build DAP2+DAP4 Clients
OPTION(ENABLE_DAP "Enable DAP2 and DAP4 Client." ON)
IF(ENABLE_DAP)
  SET(USE_DAP ON CACHE BOOL "")
  SET(ENABLE_DAP2 ON CACHE BOOL "")


  IF(ENABLE_NETCDF_4)
    SET(ENABLE_DAP4 ON CACHE BOOL "")
  ELSE(ENABLE_NETCDF_4)
    SET(ENABLE_DAP4 OFF CACHE BOOL "")
  ENDIF(ENABLE_NETCDF_4)

  FIND_PACKAGE(CURL)
  IF(NOT CURL_LIBRARY)
    MESSAGE(FATAL_ERROR "DAP Support specified, CURL libraries are not found.")
  ENDIF()

  ADD_DEFINITIONS(-DCURL_STATICLIB=1)
  INCLUDE_DIRECTORIES(${CURL_INCLUDE_DIRS})

  # Check to see if CURLOPT_USERNAME is defined.
  # It is present starting version 7.19.1.
  CHECK_C_SOURCE_COMPILES("
  #include <curl/curl.h>
  int main() {int x = CURLOPT_USERNAME;}" HAVE_CURLOPT_USERNAME)

  # Check to see if CURLOPT_PASSWORD is defined.
  # It is present starting version 7.19.1.
  CHECK_C_SOURCE_COMPILES("
  #include <curl/curl.h>
  int main() {int x = CURLOPT_PASSWORD;}" HAVE_CURLOPT_PASSWORD)

  # Check to see if CURLOPT_KEYPASSWD is defined.
  # It is present starting version 7.16.4.
  CHECK_C_SOURCE_COMPILES("
  #include <curl/curl.h>
  int main() {int x = CURLOPT_KEYPASSWD;}" HAVE_CURLOPT_KEYPASSWD)

  # Check to see if CURLINFO_RESPONSE_CODE is defined.
  # It showed up in curl 7.10.7.
  CHECK_C_SOURCE_COMPILES("
  #include <curl/curl.h>
  int main() {int x = CURLINFO_RESPONSE_CODE;}" HAVE_CURLINFO_RESPONSE_CODE)

  # Check to see if CURLINFO_HTTP_CONNECTCODE is defined.
  # It showed up in curl 7.10.7.
  CHECK_C_SOURCE_COMPILES("
  #include <curl/curl.h>
  int main() {int x = CURLINFO_HTTP_CONNECTCODE;}" HAVE_CURLINFO_HTTP_CONNECTCODE)

  # Check to see if CURLOPT_BUFFERSIZE is defined.
  # It is present starting version 7.59
  CHECK_C_SOURCE_COMPILES("
  #include <curl/curl.h>
  int main() {int x = CURLOPT_BUFFERSIZE;}" HAVE_CURLOPT_BUFFERSIZE)

  # Check to see if CURLOPT_TCP_KEEPALIVE is defined.
  # It is present starting version 7.25
  CHECK_C_SOURCE_COMPILES("
  #include <curl/curl.h>
  int main() {int x = CURLOPT_TCP_KEEPALIVE;}" HAVE_CURLOPT_KEEPALIVE)

ELSE()
  SET(ENABLE_DAP2 OFF)
  SET(ENABLE_DAP4 OFF)
ENDIF()


# Check for the math library so it can be explicitly linked.
IF(NOT WIN32)
  FIND_LIBRARY(HAVE_LIBM NAMES math m libm)
  MESSAGE(STATUS "Found Math library: ${HAVE_LIBM}")
  IF(NOT HAVE_LIBM)
    MESSAGE(FATAL_ERROR "Unable to find the math library.")
  ENDIF()
ENDIF()

# Option to Enable DAP long tests, remote tests.
OPTION(ENABLE_DAP_LONG_TESTS "Enable DAP long tests." OFF)
OPTION(ENABLE_DAP_REMOTE_TESTS "Enable DAP remote tests." ON)
SET(REMOTETESTSERVERS "149.165.169.123:8080,remotetest.unidata.ucar.edu" CACHE STRING "test servers to use for remote test")

# If netCDF4 and DAP, Option for DAP groups.
IF(ENABLE_NETCDF_4 AND ENABLE_DAP2)
  OPTION(ENABLE_DAP_GROUPS "Whether netcdf4 group names should be enabled." ON)
ELSE()
  SET(ENABLE_DAP_GROUPS OFF CACHE BOOL "Whether netcdf4 group names should be enabled.")
ENDIF()

# Enable some developer-only tests
OPTION(ENABLE_EXTRA_TESTS "Enable Extra tests. Some may not work because of known issues. Developers only." OFF)
IF(ENABLE_EXTRA_TESTS)
  SET(EXTRA_TESTS ON)
ENDIF()

# Option to use bundled XGetopt in place of getopt(). This is mostly useful
# for MSVC builds. If not building utilities, getopt() isn't required at all.
IF(MSVC)
  OPTION(ENABLE_XGETOPT "Enable bundled XGetOpt instead of external getopt()." ON)
  IF(ENABLE_XGETOPT)
    SET(USE_X_GETOPT ON CACHE BOOL "")
  ENDIF()
ENDIF()

SET(MATH "")
IF(NOT WIN32)

  # STDIO instead of posixio.
  OPTION(ENABLE_STDIO "If true, use stdio instead of posixio (ex. on the Cray)" OFF)
  IF(ENABLE_STDIO)
    SET(USE_STDIO ON CACHE BOOL "")
  ENDIF()

  # FFIO insteaad of PosixIO
  OPTION(ENABLE_FFIO "If true, use ffio instead of posixio" OFF)
  IF(ENABLE_FFIO)
    SET(USE_FFIO ON CACHE BOOL "")
  ENDIF()
ENDIF()

##
# Enable Tests
##
OPTION(ENABLE_TESTS "Enable basic tests, run with 'make test'." ON)
IF(ENABLE_TESTS)
  SET(BUILD_TESTSETS ON CACHE BOOL "")

  # Options for CTest-based tests, dashboards.
  SET(NC_CTEST_PROJECT_NAME "netcdf-c" CACHE STRING "Project Name for CTest-based testing purposes.")
  SET(NC_CTEST_DROP_SITE "cdash.unidata.ucar.edu" CACHE STRING "Dashboard location for CTest-based testing purposes.")
  SET(NC_CTEST_DROP_LOC_PREFIX "" CACHE STRING "Prefix for Dashboard location on remote server when using CTest-based testing.")

  FIND_PROGRAM(HOSTNAME_CMD NAMES hostname)
  IF(NOT MSVC)
    SET(HOSTNAME_ARG "-s")
  ENDIF()
  IF(HOSTNAME_CMD)
    EXEC_PROGRAM(${HOSTNAME_CMD} ARGS "${HOSTNAME_ARG}" OUTPUT_VARIABLE HOSTNAME)
    SET(NC_CTEST_SITE "${HOSTNAME}" CACHE STRING "Hostname of test machine.")
  ENDIF()

  IF(NC_CTEST_SITE)
    SET(SITE "${NC_CTEST_SITE}" CACHE STRING "")
  ENDIF()

  # Create a CTestConfig file from the template.
  CONFIGURE_FILE("${netCDF_SOURCE_DIR}/CTestConfig.cmake.in"
    "${netCDF_SOURCE_DIR}/CTestConfig.cmake"
    @ONLY
    )

  ###
  # This option dictates whether or not to turn on
  # tests which are known to fail.  This is not the
  # same thing as an 'expected failure'. Rather, these
  # are tests that will need to be fixed eventually.
  #
  # By placing them here, we can occasionaly turn this
  # flag on and see if any known failures have been
  # fixed in the course of code improvement/other bug
  # fixes.
  #
  # To use this, simply add as a fencepost around tests
  # which are known to fail.
  ###

  OPTION(ENABLE_FAILING_TESTS "Run tests which are known to fail, check to see if any have been fixed." OFF)

  ###
  # End known-failures.
  ###
  MARK_AS_ADVANCED(ENABLE_FAILING_TESTS)
ENDIF()

###
# Option to enable extreme numbers during testing.
###
OPTION(ENABLE_EXTREME_NUMBERS "Enable extreme numbers during testing, such as MAX_INT-1" ON)
IF(ENABLE_EXTREME_NUMBERS)
  SET(USE_EXTREME_NUMBERS ON)
ENDIF()

# Enable Large file tests
IF(ENABLE_LARGE_FILE_TESTS)
  SET(LARGE_FILE_TESTS ON)
ENDIF()

OPTION(ENABLE_METADATA_PERF_TESTS "Enable test of metadata performance." OFF)
IF(ENABLE_METADATA_PERF_TESTS)
  SET(ENABLE_METADATA_PERF ON)
ENDIF()

# Location for large file tests.
SET(TEMP_LARGE "." CACHE STRING "Location to store large file tests.")

OPTION(ENABLE_FSYNC "Enable experimental fsync code." OFF)
IF(ENABLE_FSYNC)
  SET(USE_FSYNC ON)
ENDIF()

# Temporary
OPTION (ENABLE_JNA "Enable jna bug fix code." OFF)
IF(ENABLE_JNA)
  SET(JNA ON)
ENDIF()

# Linux specific large file support flags.
# Modelled after check in CMakeLists.txt for hdf5.
OPTION(ENABLE_LARGE_FILE_SUPPORT "Enable large file support." ON)
IF(ENABLE_LARGE_FILE_SUPPORT)
  IF(MSVC)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /LARGEADDRESSAWARE")
    SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /LARGEADDRESSAWARE")
  ELSE()
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64")
  ENDIF()
ENDIF()

OPTION(ENABLE_EXAMPLE_TESTS "Run extra example tests.  Requires GNU Sed. Ignored if netCDF-4 is not Enabled" OFF)
IF(NOT ENABLE_NETCDF_4 AND ENABLE_EXAMPLE_TESTS)
  SET(ENABLE_EXAMPLE_TESTS OFF)
ENDIF()

# Enable Parallel IO with netCDF-4/HDF5 files using HDF5 parallel I/O.
SET(STATUS_PARALLEL "OFF")
OPTION(ENABLE_PARALLEL4 "Build netCDF-4 with parallel IO" "${HDF5_PARALLEL}")
IF(ENABLE_PARALLEL4 AND ENABLE_NETCDF_4)
  IF(NOT HDF5_PARALLEL)
    SET(USE_PARALLEL OFF CACHE BOOL "")
    MESSAGE(STATUS "Cannot find HDF5 library built with parallel support. Disabling parallel build.")
  ELSE()
    SET(HDF5_PARALLEL ON CACHE BOOL "")
    SET(USE_PARALLEL ON CACHE BOOL "")
    SET(USE_PARALLEL4 ON CACHE BOOL "")
    SET(STATUS_PARALLEL "ON")
    configure_file("${netCDF_SOURCE_DIR}/nc_test4/run_par_test.sh.in"
      "${netCDF_BINARY_DIR}/tmp/run_par_test.sh" @ONLY NEWLINE_STYLE LF)
    FILE(COPY "${netCDF_BINARY_DIR}/tmp/run_par_test.sh"
      DESTINATION ${netCDF_BINARY_DIR}/nc_test4
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  ENDIF()
ENDIF()

# Options to enable parallel IO for classic formats with PnetCDF library.
SET(STATUS_PNETCDF "OFF")
OPTION(ENABLE_PNETCDF "Build with parallel I/O for CDF-1, 2, and 5 files using PnetCDF." OFF)

IF(ENABLE_PNETCDF)
  # Check for ncmpi_create in libpnetcdf, define USE_PNETCDF
  # Does the user want to turn on PnetCDF read ability?
  SET(USE_PNETCDF ON CACHE BOOL "")
  FIND_LIBRARY(PNETCDF NAMES pnetcdf)
  FIND_PATH(PNETCDF_INCLUDE_DIR pnetcdf.h)
  IF(NOT PNETCDF)
    MESSAGE(STATUS "Cannot find PnetCDF library. Disabling PnetCDF support.")
    SET(USE_PNETCDF OFF CACHE BOOL "")
  ELSE(NOT PNETCDF)
    SET(USE_PARALLEL ON CACHE BOOL "")

    # Check PnetCDF version. Must be >= 1.6.0
    set(pnetcdf_h "${PNETCDF_INCLUDE_DIR}/pnetcdf.h" )
    message(STATUS "PnetCDF include file ${pnetcdf_h} will be searched for version")
    file(STRINGS "${pnetcdf_h}" pnetcdf_major_string REGEX "^#define PNETCDF_VERSION_MAJOR")
    string(REGEX REPLACE "[^0-9]" "" pnetcdf_major "${pnetcdf_major_string}")
    file(STRINGS "${pnetcdf_h}" pnetcdf_minor_string REGEX "^#define PNETCDF_VERSION_MINOR")
    string(REGEX REPLACE "[^0-9]" "" pnetcdf_minor "${pnetcdf_minor_string}")
    file(STRINGS "${pnetcdf_h}" pnetcdf_sub_string REGEX "^#define PNETCDF_VERSION_SUB")
    string(REGEX REPLACE "[^0-9]" "" pnetcdf_sub "${pnetcdf_sub_string}")
    set(pnetcdf_version "${pnetcdf_major}.${pnetcdf_minor}.${pnetcdf_sub}")
    message(STATUS "Found PnetCDF version ${pnetcdf_version}")

    if(${pnetcdf_version} VERSION_GREATER "1.6.0")
      SET(STATUS_PNETCDF "ON")
      INCLUDE_DIRECTORIES(${PNETCDF_INCLUDE_DIR})
      SET(HAVE_LIBPNETCDF ON)
      # PnetCDF => parallel
      SET(STATUS_PARALLEL ON)
      SET(USE_PARALLEL ON)
      MESSAGE(STATUS "Using PnetCDF Library: ${PNETCDF}")
    ELSE()
      MESSAGE(WARNING "ENABLE_PNETCDF requires version 1.6.1 or later; found version ${pnetcdf_version}. PnetCDF is disabled")
    ENDIF()
  ENDIF(NOT PNETCDF)
ENDIF()

# Options to enable use of fill values for elements causing NC_ERANGE
SET(ENABLE_ERANGE_FILL AUTO CACHE STRING "AUTO")
OPTION(ENABLE_ERANGE_FILL "Enable use of fill value when out-of-range type conversion causes NC_ERANGE error." OFF)
IF(ENABLE_ERANGE_FILL) # enable or auto
   STRING(TOUPPER ${ENABLE_ERANGE_FILL} ENABLE_ERANGE_FILL)
   IF(ENABLE_ERANGE_FILL AND NOT ENABLE_ERANGE_FILL STREQUAL "AUTO")
      # explicitly enabled
      SET(ENABLE_ERANGE_FILL ON)
   ELSE()
      IF(NOT ENABLE_ERANGE_FILL STREQUAL "AUTO")
         SET(ENABLE_ERANGE_FILL OFF)
      ENDIF()
   ENDIF()
ENDIF(ENABLE_ERANGE_FILL)
# Now ENABLE_ERANGE_FILL is either AUTO, ON, or OFF

# Options to use a more relaxed coordinate argument boundary check
SET(ENABLE_ZERO_LENGTH_COORD_BOUND AUTO CACHE STRING "AUTO")
OPTION(ENABLE_ZERO_LENGTH_COORD_BOUND "Enable a more relaxed boundary error check NC_EINVALCOORDS to allow coordinate start argument equal to dimension size when argument count is zero." OFF)
IF(ENABLE_ZERO_LENGTH_COORD_BOUND) # enable or auto
   STRING(TOUPPER ${ENABLE_ZERO_LENGTH_COORD_BOUND} ENABLE_ZERO_LENGTH_COORD_BOUND)
   IF(ENABLE_ZERO_LENGTH_COORD_BOUND AND NOT ENABLE_ZERO_LENGTH_COORD_BOUND STREQUAL "AUTO")
      # explicitly enabled
      SET(ENABLE_ZERO_LENGTH_COORD_BOUND ON)
   ELSE()
      IF(NOT ENABLE_ZERO_LENGTH_COORD_BOUND STREQUAL "AUTO")
         SET(ENABLE_ZERO_LENGTH_COORD_BOUND OFF)
      ENDIF()
   ENDIF()
ENDIF(ENABLE_ZERO_LENGTH_COORD_BOUND)
# Now ENABLE_ZERO_LENGTH_COORD_BOUND is either AUTO, ON, or OFF

# check and conform with PnetCDF settings on ERANGE_FILL and RELAX_COORD_BOUND
IF(STATUS_PNETCDF)
  file(STRINGS "${pnetcdf_h}" enable_erange_fill_pnetcdf REGEX "^#define PNETCDF_ERANGE_FILL")
  string(REGEX REPLACE "[^0-9]" "" erange_fill_pnetcdf "${enable_erange_fill_pnetcdf}")
  IF("x${erange_fill_pnetcdf}" STREQUAL "x1")
     SET(erange_fill_pnetcdf "ON")
  ELSE()
     SET(erange_fill_pnetcdf "OFF")
  ENDIF()
  IF(ENABLE_ERANGE_FILL STREQUAL "AUTO") # not set on command line
     SET(ENABLE_ERANGE_FILL "${erange_fill_pnetcdf}")
  ELSE()
     # user explicitly set this option on command line
     IF(NOT ENABLE_ERANGE_FILL STREQUAL "${erange_fill_pnetcdf}")
        IF(ENABLE_ERANGE_FILL)
           MESSAGE(FATAL_ERROR "Enabling erange-fill conflicts with PnetCDF setting")
        ELSE()
           MESSAGE(FATAL_ERROR "Disabling erange-fill conflicts with PnetCDF setting")
        ENDIF()
     ENDIF()
  ENDIF()

  file(STRINGS "${pnetcdf_h}" relax_coord_bound_pnetcdf REGEX "^#define PNETCDF_RELAX_COORD_BOUND")
  string(REGEX REPLACE "[^0-9]" "" relax_coord_bound "${relax_coord_bound_pnetcdf}")
  IF("x${relax_coord_bound}" STREQUAL "x1")
     SET(relax_coord_bound_pnetcdf "ON")
  ELSE()
     SET(relax_coord_bound_pnetcdf "OFF")
  ENDIF()
  IF(ENABLE_ZERO_LENGTH_COORD_BOUND STREQUAL "AUTO") # not set on command line
     SET(ENABLE_ZERO_LENGTH_COORD_BOUND "${relax_coord_bound_pnetcdf}")
  ELSE()
     # user explicitly set this option on command line
     IF(NOT ENABLE_ZERO_LENGTH_COORD_BOUND STREQUAL "${relax_coord_bound_pnetcdf}")
        IF(ENABLE_ZERO_LENGTH_COORD_BOUND)
           MESSAGE(FATAL_ERROR "Enabling relax-coord-bound conflicts with PnetCDF setting")
        ELSE()
           MESSAGE(FATAL_ERROR "Disabling relax-coord-bound conflicts with PnetCDF setting")
        ENDIF()
     ENDIF()
  ENDIF()
ENDIF()

IF(ENABLE_ERANGE_FILL)
  MESSAGE(STATUS "Enabling use of fill value when NC_ERANGE")
  SET(M4FLAGS "-DERANGE_FILL" CACHE STRING "")
ENDIF()

IF(ENABLE_ZERO_LENGTH_COORD_BOUND)
  MESSAGE(STATUS "Enabling a more relaxed check for NC_EINVALCOORDS")
  ADD_DEFINITIONS(-DRELAX_COORD_BOUND)
ENDIF()

# Enable Parallel Tests.
OPTION(ENABLE_PARALLEL_TESTS "Enable Parallel IO Tests. Requires HDF5/NetCDF4 with parallel I/O Support." ${USE_PARALLEL})
IF(ENABLE_PARALLEL_TESTS AND USE_PARALLEL)
  SET(TEST_PARALLEL ON CACHE BOOL "")
  IF(USE_NETCDF4)
    SET(TEST_PARALLEL4 ON CACHE BOOL "")
  ENDIF()
ENDIF()

IF (ENABLE_PARALLEL_TESTS AND NOT USE_PARALLEL)
  MESSAGE(FATAL_ERROR "Parallel tests requested, but no parallel HDF5 installation detected.")
ENDIF()

# Enable special filter test; experimental when using cmake.
OPTION(ENABLE_FILTER_TESTING "Enable filter testing. Ignored if shared libraries or netCDF4 are not enabled" OFF)
IF(NOT ENABLE_NETCDF4)
  MESSAGE(WARNING "ENABLE_FILTER_TESTING requires netCDF-4. Disabling.")
  SET(ENABLE_FILTER_TESTING OFF CACHE BOOL "")
ENDIF()
IF(NOT BUILD_SHARED_LIBS)
  MESSAGE(WARNING "ENABLE_FILTER_TESTING requires shared libraries. Disabling.")
  SET(ENABLE_FILTER_TESTING OFF CACHE BOOL "")
ENDIF()

# Determine whether or not to generate documentation.
OPTION(ENABLE_DOXYGEN "Enable generation of doxygen-based documentation." OFF)
IF(ENABLE_DOXYGEN)
  FIND_PACKAGE(Doxygen REQUIRED)
  # Offer the option to build internal documentation.
  OPTION(ENABLE_INTERNAL_DOCS "Build internal documentation. This is of interest to developers only." OFF)
  IF(ENABLE_INTERNAL_DOCS)
    SET(BUILD_INTERNAL_DOCS YES CACHE STRING "")
  ELSE()
    SET(BUILD_INTERNAL_DOCS NO CACHE STRING "")
  ENDIF()

  ###
  #
  # If we are building release documentation, we need to set some
  # variables that will be used in the Doxygen.in template.
  ###
  OPTION(ENABLE_DOXYGEN_BUILD_RELEASE_DOCS "Build release documentation.  This is of interest only to the netCDF developers." OFF)
  IF(ENABLE_DOXYGEN_BUILD_RELEASE_DOCS)
    SET(DOXYGEN_CSS_FILE "${CMAKE_SOURCE_DIR}/docs/release.css" CACHE STRING "")
    SET(DOXYGEN_HEADER_FILE "${CMAKE_SOURCE_DIR}/docs/release_header.html" CACHE STRING "")
    SET(DOXYGEN_SEARCHENGINE "NO" CACHE STRING "")
    SET(ENABLE_DOXYGEN_SERVER_BASED_SEARCH NO CACHE STRING "")
  ELSE()
    SET(DOXYGEN_CSS_FILE "" CACHE STRING "")
    SET(DOXYGEN_HEADER_FILE "" CACHE STRING "")
    SET(DOXYGEN_SEARCHENGINE "YES" CACHE STRING "")

    # If not using release document configuration,
    # provide an option for server-based search.
    OPTION(ENABLE_DOXYGEN_SERVER_SIDE_SEARCH "Configure Doxygen with server-based search." OFF)
    IF(ENABLE_DOXYGEN_SERVER_SIDE_SEARCH)
      SET(DOXYGEN_SERVER_BASED_SEARCH "YES" CACHE STRING "")
    ELSE()
      SET(DOXYGEN_SERVER_BASED_SEARCH "NO" CACHE STRING "")
    ENDIF(ENABLE_DOXYGEN_SERVER_SIDE_SEARCH)


  ENDIF(ENABLE_DOXYGEN_BUILD_RELEASE_DOCS)
  # Option to turn on the TODO list in the doxygen-generated documentation.
  OPTION(DOXYGEN_ENABLE_TASKS "Turn on test, todo, bug lists in documentation. This is of interest to developers only." OFF)
  IF(DOXYGEN_ENABLE_TASKS)
    SET(SHOW_DOXYGEN_TAG_LIST YES CACHE STRING "")
  ELSE(DOXYGEN_ENABLE_TASKS)
    SET(SHOW_DOXYGEN_TODO_LIST NO CACHE STRING "")
  ENDIF(DOXYGEN_ENABLE_TASKS)

  OPTION(ENABLE_DOXYGEN_PDF_OUTPUT "[EXPERIMENTAL] Turn on PDF output for Doxygen-generated documentation." OFF)

  IF(ENABLE_DOXYGEN_PDF_OUTPUT)
    SET(NC_ENABLE_DOXYGEN_PDF_OUTPUT "YES" CACHE STRING "")
  ELSE()
    SET(NC_ENABLE_DOXYGEN_PDF_OUTPUT "NO" CACHE STRING "")
  ENDIF()

  FIND_PROGRAM(NC_DOT NAMES dot)
  # Specify whether or not 'dot' was found on the system path.
  IF(NC_DOT)
    SET(HAVE_DOT YES CACHE STRING "")
  ELSE(NC_DOT)
    SET(HAVE_DOT NO CACHE STRING "")
  ENDIF(NC_DOT)
ENDIF()

# Always enable DISKLESS
OPTION(ENABLE_DISKLESS "Enable in-memory files" ON)


# By default, MSVC has a stack size of 1000000.
# Allow a user to override this.
IF(MSVC)
  SET(NC_MSVC_STACK_SIZE 40000000 CACHE STRING "Default stack size for MSVC-based projects.")
  # By default, CMake sets the stack to 1000000.
  # Remove this limitation.
  # See here for more details:
  # http://www.cmake.org/pipermail/cmake/2009-April/028710.html
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:${NC_MSVC_STACK_SIZE}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /STACK:${NC_MSVC_STACK_SIZE}")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /STACK:${NC_MSVC_STACK_SIZE}")
ENDIF()

# Set some of the options as advanced.
MARK_AS_ADVANCED(ENABLE_INTERNAL_DOCS VALGRIND_TESTS ENABLE_COVERAGE_TESTS )
MARK_AS_ADVANCED(ENABLE_DAP_REMOTE_TESTS ENABLE_DAP_LONG_TESTS USE_REMOTE_CDASH)
MARK_AS_ADVANCED(ENABLE_DOXYGEN_BUILD_RELEASE_DOCS DOXYGEN_ENABLE_TASKS ENABLE_DOXYGEN_SERVER_SIDE_SEARCH)
MARK_AS_ADVANCED(ENABLE_SHARED_LIBRARY_VERSION)

################################
# Option checks
################################

# Library include checks
CHECK_INCLUDE_FILE("math.h"      HAVE_MATH_H)
CHECK_INCLUDE_FILE("unistd.h"  HAVE_UNISTD_H)
# Solve a compatibility issue in ncgen/, which checks
# for NO_UNISTD_H
IF(NOT HAVE_UNISTD_H)
  SET(YY_NO_UNISTD_H TRUE)
ENDIF()

CHECK_INCLUDE_FILE("alloca.h"  HAVE_ALLOCA_H)
CHECK_INCLUDE_FILE("malloc.h"    HAVE_MALLOC_H)
CHECK_INCLUDE_FILE("fcntl.h"   HAVE_FCNTL_H)
CHECK_INCLUDE_FILE("getopt.h"  HAVE_GETOPT_H)
CHECK_INCLUDE_FILE("locale.h"  HAVE_LOCALE_H)
CHECK_INCLUDE_FILE("stdint.h"  HAVE_STDINT_H)
CHECK_INCLUDE_FILE("stdio.h"   HAVE_STDIO_H)
IF(MSVC)
CHECK_INCLUDE_FILE("io.h"      HAVE_IO_H)
ENDIF(MSVC)
CHECK_INCLUDE_FILE("stdlib.h"  HAVE_STDLIB_H)
CHECK_INCLUDE_FILE("stdarg.h"    HAVE_STDARG_H)
CHECK_INCLUDE_FILE("strings.h"   HAVE_STRINGS_H)
CHECK_INCLUDE_FILE("signal.h"    HAVE_SIGNAL_H)
CHECK_INCLUDE_FILE("sys/param.h" HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILE("sys/stat.h"  HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE("sys/time.h"  HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE("sys/types.h" HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE("sys/mman.h"  HAVE_SYS_MMAN_H)
CHECK_INCLUDE_FILE("sys/resource.h" HAVE_SYS_RESOURCE_H)
CHECK_INCLUDE_FILE("fcntl.h"  HAVE_FCNTL_H)
CHECK_INCLUDE_FILE("inttypes.h"  HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE("pstdint.h"  HAVE_PSTDINT_H)
CHECK_INCLUDE_FILE("endian.h" HAVE_ENDIAN_H)
CHECK_INCLUDE_FILE("BaseTsd.h"  HAVE_BASETSD_H)
CHECK_INCLUDE_FILE("stddef.h"   HAVE_STDDEF_H)
CHECK_INCLUDE_FILE("string.h"  HAVE_STRING_H)
CHECK_INCLUDE_FILE("winsock2.h" HAVE_WINSOCK2_H)
CHECK_INCLUDE_FILE("ftw.h"  HAVE_FTW_H)
CHECK_INCLUDE_FILE("libgen.h" HAVE_LIBGEN_H)

# Symbol Exists
CHECK_SYMBOL_EXISTS(isfinite "math.h" HAVE_DECL_ISFINITE)
CHECK_SYMBOL_EXISTS(isnan "math.h" HAVE_DECL_ISNAN)
CHECK_SYMBOL_EXISTS(isinf "math.h" HAVE_DECL_ISINF)
CHECK_SYMBOL_EXISTS(st_blksize "sys/stat.h" HAVE_STRUCT_STAT_ST_BLKSIZE)
CHECK_SYMBOL_EXISTS(alloca "alloca.h" HAVE_ALLOCA)
CHECK_SYMBOL_EXISTS(snprintf "stdio.h" HAVE_SNPRINTF)

# Type checks
# Aliases for automake consistency
SET(SIZEOF_VOIDSTAR ${CMAKE_SIZEOF_VOID_P})
SET(SIZEOF_VOIDP ${SIZEOF_VOIDSTAR})
CHECK_TYPE_SIZE("char"      SIZEOF_CHAR)
CHECK_TYPE_SIZE("double"    SIZEOF_DOUBLE)
CHECK_TYPE_SIZE("float"     SIZEOF_FLOAT)
CHECK_TYPE_SIZE("int"       SIZEOF_INT)
CHECK_TYPE_SIZE("uint"      SIZEOF_UINT)
IF(SIZEOF_UINT)
  SET(HAVE_UINT TRUE)
ENDIF(SIZEOF_UINT)

CHECK_TYPE_SIZE("schar"      SIZEOF_SCHAR)
IF(SIZEOF_SCHAR)
  SET(HAVE_SCHAR TRUE)
ENDIF(SIZEOF_SCHAR)

CHECK_TYPE_SIZE("long"      SIZEOF_LONG)
CHECK_TYPE_SIZE("long long" SIZEOF_LONG_LONG)
IF(SIZEOF_LONG_LONG)
  SET(HAVE_LONG_LONG_INT TRUE)
ENDIF(SIZEOF_LONG_LONG)

CHECK_TYPE_SIZE("unsigned long long" SIZEOF_UNSIGNED_LONG_LONG)

CHECK_TYPE_SIZE("off_t"     SIZEOF_OFF_T)
CHECK_TYPE_SIZE("off64_t"   SIZEOF_OFF64_T)
CHECK_TYPE_SIZE("short"     SIZEOF_SHORT)
CHECK_TYPE_SIZE("ushort"    SIZEOF_USHORT)
IF(SIZEOF_USHORT)
  SET(HAVE_USHORT TRUE)
ENDIF(SIZEOF_USHORT)

CHECK_TYPE_SIZE("_Bool"     SIZEOF__BOOL)

CHECK_TYPE_SIZE("size_t"    SIZEOF_SIZE_T)

# Check whether to turn on or off CDF5 support.
SET(ENABLE_CDF5 AUTO CACHE STRING "AUTO")
OPTION(ENABLE_CDF5 "Enable CDF5 support" ON)
IF(SIZEOF_SIZE_T EQUAL 4)
   IF(ENABLE_CDF5) # enable or auto
      STRING(TOUPPER ${ENABLE_CDF5} ENABLE_CDF5)
      IF(ENABLE_CDF5 AND NOT ENABLE_CDF5 STREQUAL "AUTO") # explicitly enabled
         MESSAGE(FATAL_ERROR "Unable to support CDF5 feature because size_t is less than 8 bytes")
      ENDIF(ENABLE_CDF5 AND NOT ENABLE_CDF5 STREQUAL "AUTO")
      SET(ENABLE_CDF5 OFF) # cannot support CDF5
      SET(USE_CDF5 OFF CACHE BOOL "") # cannot support CDF5
   ENDIF(ENABLE_CDF5)
ELSE(SIZEOF_SIZE_T EQUAL 4)
   IF(ENABLE_CDF5) # explicitly set by user or not set
      SET(USE_CDF5 ON CACHE BOOL "")
   ELSE(ENABLE_CDF5) # explicitly disabled by user
      SET(USE_CDF5 OFF CACHE BOOL "")
   ENDIF(ENABLE_CDF5)
ENDIF(SIZEOF_SIZE_T EQUAL 4)

CHECK_TYPE_SIZE("ssize_t"   SIZEOF_SSIZE_T)
IF(SIZEOF_SSIZE_T)
  SET(HAVE_SSIZE_T TRUE)
ENDIF(SIZEOF_SSIZE_T)
CHECK_TYPE_SIZE("ptrdiff_t" SIZEOF_PTRDIFF_T)
IF(SIZEOF_PTRDIFF_T)
  SET(HAVE_PTRDIFF_T TRUE)
ENDIF(SIZEOF_PTRDIFF_T)
CHECK_TYPE_SIZE("uintptr_t" SIZEOF_UINTPTR_T)
IF(SIZEOF_UINTPTR_T)
  SET(HAVE_UINTPTR_T TRUE)
ENDIF(SIZEOF_UINTPTR_T)

# __int64 is used on Windows for large file support.
CHECK_TYPE_SIZE("__int64"   SIZEOF___INT_64)
CHECK_TYPE_SIZE("int64_t"   SIZEOF_INT64_T)
CHECK_TYPE_SIZE("uint64_t"  SIZEOF_UINT64_T)
CHECK_TYPE_SIZE("unsigned char"      SIZEOF_UCHAR)
CHECK_TYPE_SIZE("unsigned short int" SIZEOF_UNSIGNED_SHORT_INT)
CHECK_TYPE_SIZE("unsigned int"       SIZEOF_UNSIGNED_INT)
CHECK_TYPE_SIZE("long long"          SIZEOF_LONGLONG)
CHECK_TYPE_SIZE("unsigned long long" SIZEOF_ULONGLONG)

# On windows systems, we redefine off_t as __int64
# to enable LFS. This is true on 32 and 64 bit system.s
# We must redefine SIZEOF_OFF_T to match.
IF(MSVC AND SIZEOF___INT_64)
  SET(SIZEOF_OFF_T  ${SIZEOF___INT_64})
ENDIF()

# Check for various functions.
CHECK_FUNCTION_EXISTS(fsync HAVE_FSYNC)
CHECK_FUNCTION_EXISTS(strlcat   HAVE_STRLCAT)
CHECK_FUNCTION_EXISTS(strdup  HAVE_STRDUP)
CHECK_FUNCTION_EXISTS(strndup HAVE_STRNDUP)
CHECK_FUNCTION_EXISTS(strtoll HAVE_STRTOLL)
CHECK_FUNCTION_EXISTS(strtoull  HAVE_STRTOULL)
CHECK_FUNCTION_EXISTS(mkstemp HAVE_MKSTEMP)
CHECK_FUNCTION_EXISTS(mktemp HAVE_MKTEMP)
CHECK_FUNCTION_EXISTS(random HAVE_RANDOM)
CHECK_FUNCTION_EXISTS(gettimeofday  HAVE_GETTIMEOFDAY)
CHECK_FUNCTION_EXISTS(MPI_Comm_f2c  HAVE_MPI_COMM_F2C)
CHECK_FUNCTION_EXISTS(MPI_Info_f2c  HAVE_MPI_INFO_F2C)
CHECK_FUNCTION_EXISTS(memmove HAVE_MEMMOVE)
CHECK_FUNCTION_EXISTS(getpagesize HAVE_GETPAGESIZE)
CHECK_FUNCTION_EXISTS(sysconf HAVE_SYSCONF)
CHECK_FUNCTION_EXISTS(getrlimit HAVE_GETRLIMIT)
CHECK_FUNCTION_EXISTS(_filelengthi64 HAVE_FILE_LENGTH_I64)
CHECK_FUNCTION_EXISTS(mmap HAVE_MMAP)
CHECK_FUNCTION_EXISTS(mremap HAVE_MREMAP)

IF(NOT HAVE_MMMAP OR NOT HAVE_MREMAP)
  MESSAGE(WARNING "mmap or mremap not found: disabling MMAP support.")
  SET(ENABLE_MMAP OFF)
ENDIF()

IF(ENABLE_MMAP)
  SET(BUILD_MMAP ON)
  SET(USE_MMAP ON)
ENDIF(ENABLE_MMAP)

#CHECK_FUNCTION_EXISTS(alloca HAVE_ALLOCA)
#####
# End system inspection checks.
#####

################################
# Define Utility Macros
################################

# Macro to append files to the EXTRA_DIST files.
# Note: can only be used in subdirectories because of the use of PARENT_SCOPE
SET(EXTRA_DIST "")
MACRO(ADD_EXTRA_DIST files)
  FOREACH(F ${files})
    SET(EXTRA_DIST ${EXTRA_DIST} ${CMAKE_CURRENT_SOURCE_DIR}/${F})
    SET(EXTRA_DIST ${EXTRA_DIST} PARENT_SCOPE)
  ENDFOREACH()
ENDMACRO()

# A basic script used to convert m4 files
FIND_PROGRAM(NC_M4 NAMES m4 m4.exe)
IF(NC_M4)
  MESSAGE(STATUS "Found m4: ${NC_M4}")
  SET(HAVE_M4 TRUE)
ELSE()
  MESSAGE(STATUS "m4 not found.")
ENDIF()

MACRO(GEN_m4 filename)

  IF(HAVE_M4)

  # If m4 is available, remove generated file if it exists.
  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c)
    FILE(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c)
  ENDIF()

  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c
    COMMAND ${NC_M4}
    ARGS ${M4FLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/${filename}.c
    VERBATIM
    )
ENDIF(HAVE_M4)
ENDMACRO(GEN_m4)

# Binary tests, but ones which depend on value of 'TEMP_LARGE' being defined.
MACRO(add_bin_env_temp_large_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdf)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()

  ADD_TEST(${prefix}_${F} bash "-c" "TEMP_LARGE=${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
    SET_TARGET_PROPERTIES(${prefix}_${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${prefix}_${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${prefix}_${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
      ${CMAKE_CURRENT_BINARY_DIR})
   SET_TARGET_PROPERTIES(${prefix}_${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO
      ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF()
ENDMACRO()


# Tests which are binary, but depend on a particular environmental variable.
MACRO(add_bin_env_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F} netcdf)
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()

  ADD_TEST(${prefix}_${F} bash "-c" "TOPSRCDIR=${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${prefix}_${F}")
  IF(MSVC)
    SET_PROPERTY(TARGET ${prefix}_${F} PROPERTY FOLDER "tests")
  ENDIF()
ENDMACRO()

# Build a binary used by a script, but don't make a test out of it.
MACRO(build_bin_test F)
  ADD_EXECUTABLE(${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${F} netcdf ${ALL_TLL_LIBS})
  IF(MSVC)
    SET_TARGET_PROPERTIES(${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
    SET_TARGET_PROPERTIES(${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
      ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF()
ENDMACRO()

# Binary tests which are used by a script looking for a specific name.
MACRO(add_bin_test_no_prefix F)
  build_bin_test(${F})
  ADD_TEST(${F} ${EXECUTABLE_OUTPUT_PATH}/${F})
  IF(MSVC)
    SET_PROPERTY(TEST ${F} PROPERTY FOLDER "tests/")
    SET_TARGET_PROPERTIES(${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
      ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF()
ENDMACRO()

# Binary tests which are used by a script looking for a specific name.
MACRO(build_bin_test_no_prefix F)
  build_bin_test(${F})
  IF(MSVC)
    #SET_PROPERTY(TEST ${F} PROPERTY FOLDER "tests/")
    SET_TARGET_PROPERTIES(${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
      ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF()
ENDMACRO()

MACRO(add_bin_test prefix F)
  ADD_EXECUTABLE(${prefix}_${F} ${F}.c)
  TARGET_LINK_LIBRARIES(${prefix}_${F}
    ${ALL_TLL_LIBS}
    netcdf
    )
  IF(MSVC)
    SET_TARGET_PROPERTIES(${prefix}_${F}
      PROPERTIES LINK_FLAGS_DEBUG " /NODEFAULTLIB:MSVCRT"
      )
  ENDIF()
  ADD_TEST(${prefix}_${F}
           ${EXECUTABLE_OUTPUT_PATH}/${prefix}_${F}
           )
  IF(MSVC)
    SET_PROPERTY(TEST ${prefix}_${F} PROPERTY FOLDER "tests/")
    SET_TARGET_PROPERTIES(${prefix}_${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${prefix}_${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
      ${CMAKE_CURRENT_BINARY_DIR})
    SET_TARGET_PROPERTIES(${prefix}_${F} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
      ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF()
ENDMACRO()

# A cmake script to print out information at the end of the configuration step.
MACRO(print_conf_summary)
  MESSAGE("")
  MESSAGE("")
  MESSAGE("Configuration Summary:")
  MESSAGE("")
  MESSAGE(STATUS "Building Shared Libraries:     ${BUILD_SHARED_LIBS}")
  MESSAGE(STATUS "Building netCDF-4:             ${ENABLE_NETCDF_4}")
  MESSAGE(STATUS "Building DAP2 Support:         ${ENABLE_DAP2}")
  MESSAGE(STATUS "Building DAP4 Support:         ${ENABLE_DAP4}")
  MESSAGE(STATUS "Building Utilities:            ${BUILD_UTILITIES}")
  IF(CMAKE_PREFIX_PATH)
    MESSAGE(STATUS "CMake Prefix Path:             ${CMAKE_PREFIX_PATH}")
  ENDIF()
  MESSAGE("")

  IF(${STATUS_PNETCDF} OR ${STATUS_PARALLEL})
    MESSAGE("Building Parallel NetCDF")
    MESSAGE(STATUS "Using PnetCDF:       ${STATUS_PNETCDF}")
    MESSAGE(STATUS "Using Parallel IO:   ${STATUS_PARALLEL}")
    MESSAGE("")
  ENDIF()

  MESSAGE("Tests Enabled:              ${ENABLE_TESTS}")
  IF(ENABLE_TESTS)
    MESSAGE(STATUS "DAP Remote Tests:  ${ENABLE_DAP_REMOTE_TESTS}")
    MESSAGE(STATUS "Extra Tests:       ${ENABLE_EXTRA_TESTS}")
    MESSAGE(STATUS "Coverage Tests:    ${ENABLE_COVERAGE_TESTS}")
    MESSAGE(STATUS "Parallel Tests:    ${ENABLE_PARALLEL_TESTS}")
    MESSAGE(STATUS "Large File Tests:  ${ENABLE_LARGE_FILE_TESTS}")
    MESSAGE(STATUS "Extreme Numbers:   ${ENABLE_EXTREME_NUMBERS}")
  ENDIF()

  MESSAGE("")
  MESSAGE("Compiler:")
  MESSAGE("")
  MESSAGE(STATUS "Build Type:           ${CMAKE_BUILD_TYPE}")
  MESSAGE(STATUS "CMAKE_C_COMPILER:     ${CMAKE_C_COMPILER}")
  MESSAGE(STATUS "CMAKE_C_FLAGS:        ${CMAKE_C_FLAGS}")
  IF("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
    MESSAGE(STATUS "CMAKE_C_FLAGS_DEBUG:  ${CMAKE_C_FLAGS_DEBUG}")
  ENDIF()
  IF("${CMAKE_BUILD_TYPE}" STREQUAL "RELEASE")
    MESSAGE(STATUS "CMAKE_C_FLAGS_RELEASE:   ${CMAKE_C_FLAGS_RELEASE}")
  ENDIF()

  MESSAGE(STATUS "Linking against:      ${ALL_TLL_LIBS}")

  MESSAGE("")
ENDMACRO()
##
# Shell script Macro
##
# Determine if 'bash' is on the system.
##

OPTION(ENABLE_BASH_SCRIPT_TESTING "Detection is typically automatic, but this option can be used to force enable/disable bash-script based tests." ON)

IF(ENABLE_BASH_SCRIPT_TESTING)
  FIND_PROGRAM(HAVE_BASH bash)
  IF(HAVE_BASH)
    STRING(COMPARE EQUAL "${HAVE_BASH}" "C:/Windows/System32/bash.exe" IS_BASH_EXE)
    IF(NOT IS_BASH_EXE)
      MESSAGE(STATUS "Found bash: ${HAVE_BASH}")
    ELSE()
      MESSAGE(STATUS "Ignoring ${HAVE_BASH}")
      SET(HAVE_BASH "")
    ENDIF()
  ELSE()
    MESSAGE(STATUS "Bash shell not found; disabling shell script tests.")
  ENDIF()
ELSE(ENABLE_BASH_SCRIPT_TESTING)
  SET(HAVE_BASH "")
ENDIF(ENABLE_BASH_SCRIPT_TESTING)

MACRO(add_sh_test prefix F)
  IF(HAVE_BASH)
    ADD_TEST(${prefix}_${F} bash "-c" "export srcdir=${CMAKE_CURRENT_SOURCE_DIR};export TOPSRCDIR=${CMAKE_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}/${F}.sh")
  ENDIF()
ENDMACRO()



# A function used to create autotools-style 'yes/no' definitions.
# If a variable is set, it 'yes' is returned. Otherwise, 'no' is
# returned.
#
# Also creates a version of the ret_val prepended with 'NC',
# when feature is true, which is used to generate netcdf_meta.h.
FUNCTION(is_enabled feature ret_val)
  IF(${feature})
    SET(${ret_val} "yes" PARENT_SCOPE)
    SET("NC_${ret_val}" 1 PARENT_SCOPE)
  ELSE()
    SET(${ret_val} "no" PARENT_SCOPE)
    SET("NC_${ret_val}" 0 PARENT_SCOPE)
  ENDIF(${feature})
ENDFUNCTION()

# A function used to create autotools-style 'yes/no' definitions.
# If a variable is set, it 'yes' is returned. Otherwise, 'no' is
# returned.
#
# Also creates a version of the ret_val prepended with 'NC',
# when feature is true, which is used to generate netcdf_meta.h.
FUNCTION(is_disabled feature ret_val)
  IF(${feature})
    SET(${ret_val} "no" PARENT_SCOPE)
  ELSE()
    SET(${ret_val} "yes" PARENT_SCOPE)
    SET("NC_${ret_val}" 1 PARENT_SCOPE)
  ENDIF(${feature})
ENDFUNCTION()

################################
# End Macro Definitions
################################

# Create config.h file.
configure_file("${netCDF_SOURCE_DIR}/config.h.cmake.in"
  "${netCDF_BINARY_DIR}/config.h")
ADD_DEFINITIONS(-DHAVE_CONFIG_H)
INCLUDE_DIRECTORIES(${netCDF_BINARY_DIR})
# End autotools-style checs for config.h

#####
# Set core names of the libraries.
#####
SET(netCDF_LIB_CORENAME  "netcdf")

#####
# Set the true names of all the libraries, if customized by external project
#####
# Recurse into other subdirectories.
add_subdirectory("include")
add_subdirectory(libdispatch)
add_subdirectory(libsrc)

IF(ENABLE_FILTER_TESTING)
add_subdirectory(plugins)
ENDIF()

IF(USE_PNETCDF)
  add_subdirectory(libsrcp)
ENDIF(USE_PNETCDF)

IF(USE_HDF5)
  add_subdirectory(libsrc4)
  add_subdirectory(libhdf5)
ENDIF(USE_HDF5)

IF(USE_HDF4)
  add_subdirectory(libhdf4)
  add_subdirectory(hdf4_test)
ENDIF(USE_HDF4)

IF(ENABLE_DAP2)
  ADD_SUBDIRECTORY(oc2)
  ADD_SUBDIRECTORY(libdap2)
ENDIF()

IF(ENABLE_DAP4)
  ADD_SUBDIRECTORY(libdap4)
ENDIF()

add_subdirectory(liblib)

# For tests and utilities, we are no longer
# exporting symbols but rather importing them.
IF(BUILD_DLL)
  REMOVE_DEFINITIONS(-DDLL_EXPORT)
ENDIF()

# Enable Utilities.
IF(BUILD_UTILITIES)
  INCLUDE_DIRECTORIES(ncdump)
  ADD_SUBDIRECTORY(ncgen)
  ADD_SUBDIRECTORY(ncgen3)
  ADD_SUBDIRECTORY(ncdump)
ENDIF()

# Enable tests
IF(ENABLE_TESTS)
  IF(ENABLE_V2_API)
    ADD_SUBDIRECTORY(nctest)
  ENDIF()
  ADD_SUBDIRECTORY(nc_test)
  IF(USE_NETCDF4)
    INCLUDE_DIRECTORIES(h5_test)
    ADD_SUBDIRECTORY(nc_test4)
    ADD_SUBDIRECTORY(h5_test)
  ENDIF()
  IF(ENABLE_DAP2)
    ADD_SUBDIRECTORY(ncdap_test)
  ENDIF()
  IF(ENABLE_DAP4)
    ADD_SUBDIRECTORY(dap4_test)
  ENDIF()
  IF(ENABLE_EXAMPLES)
    ADD_SUBDIRECTORY(examples)
  ENDIF()
ENDIF()

# Code to generate an export header
#GENERATE_EXPORT_HEADER(netcdf
# BASE_NAME netcdf
# EXPORT_MACRO_NAME netcdf_EXPORT
# EXPORT_FILE_NAME netcdf_Export.h
# STATIC_DEFINE netcdf_BUILT_AS_STATIC
#)

#####
# Build doxygen documentation, if need be.
#####
ADD_SUBDIRECTORY(docs)

##
# Brute force, grab all of the dlls from the dependency directory,
# install them in the binary dir. Grab all of the .libs, put them
# in the libdir.
##
IF(MSVC)
  FILE(GLOB COPY_FILES ${CMAKE_PREFIX_PATH}/lib/*.lib)
  INSTALL(FILES ${COPY_FILES}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  COMPONENT dependencies)

  FILE(GLOB COPY_FILES ${CMAKE_PREFIX_PATH}/bin/*.dll)
  STRING(REGEX REPLACE "msv[.*].dll" "" COPY_FILES "${COPY_FILES}")
  INSTALL(FILES ${COPY_FILES}
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT dependencies)

ENDIF()

# Subdirectory CMakeLists.txt files should specify their own
# 'install' files.
# Including 'CPack' kicks everything off.
INCLUDE(InstallRequiredSystemLibraries)
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/FixBundle.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/FixBundle.cmake
  @ONLY
  )

###
# Create pkgconfig files.
###

IF(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  SET(CMAKE_INSTALL_LIBDIR lib)
ENDIF(NOT DEFINED CMAKE_INSTALL_LIBDIR)

# Set
SET(prefix ${CMAKE_INSTALL_PREFIX})
SET(exec_prefix ${CMAKE_INSTALL_PREFIX})
SET(libdir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
SET(includedir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
SET(CC ${CMAKE_C_COMPILER})

# Process all dependency libraries and create a string
# used when parsing netcdf.pc.in

FOREACH(_LIB ${ALL_TLL_LIBS})
  GET_FILENAME_COMPONENT(_LIB_NAME ${_LIB} NAME_WE)
  STRING(REGEX REPLACE "^lib" "" _NAME ${_LIB_NAME})
  LIST(APPEND NC_LIBS "-l${_NAME}")
  GET_FILENAME_COMPONENT(_LIB_DIR ${_LIB} PATH)
  LIST(APPEND LINKFLAGS "-L${_LIB_DIR}")
ENDFOREACH()

SET(NC_LIBS "-lnetcdf ${NC_LIBS}")

STRING(REPLACE ";" " " NC_LIBS "${NC_LIBS}")
STRING(REPLACE "-lhdf5::hdf5-shared" "-lhdf5" NC_LIBS ${NC_LIBS})
STRING(REPLACE "-lhdf5::hdf5_hl-shared" "-lhdf5_hl" NC_LIBS ${NC_LIBS})
STRING(REPLACE "-lhdf5::hdf5-static" "-lhdf5" NC_LIBS ${NC_LIBS})
STRING(REPLACE "-lhdf5::hdf5_hl-static" "-lhdf5_hl" NC_LIBS ${NC_LIBS})


STRING(REPLACE ";" " " LINKFLAGS "${LINKFLAGS}")

LIST(REMOVE_DUPLICATES NC_LIBS)
LIST(REMOVE_DUPLICATES LINKFLAGS)

configure_file(
  ${netCDF_SOURCE_DIR}/netcdf.pc.in
  ${netCDF_BINARY_DIR}/netcdf.pc @ONLY)


IF(NOT IS_DIRECTORY ${netCDF_BINARY_DIR}/tmp)
  FILE(MAKE_DIRECTORY ${netCDF_BINARY_DIR}/tmp)
ENDIF()

configure_file("${netCDF_SOURCE_DIR}/nc-config.cmake.in"
  "${netCDF_BINARY_DIR}/tmp/nc-config" @ONLY
  NEWLINE_STYLE LF)
FILE(COPY "${netCDF_BINARY_DIR}/tmp/nc-config"
  DESTINATION ${netCDF_BINARY_DIR}/
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

INSTALL(FILES ${netCDF_BINARY_DIR}/netcdf.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
  COMPONENT utilities)

INSTALL(PROGRAMS ${netCDF_BINARY_DIR}/nc-config
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT utilities)

###
# End pkgconfig, nc-config file creation.
###

##
# Print the configuration summary
##
print_conf_summary()

# Enable Makedist files.
ADD_MAKEDIST()
ENABLE_MAKEDIST(README.md COPYRIGHT RELEASE_NOTES.md INSTALL INSTALL.cmake test_prog.c lib_flags.am cmake CMakeLists.txt COMPILE.cmake.txt config.h.cmake.in cmake_uninstall.cmake.in netcdf-config-version.cmake.in netcdf-config.cmake.in FixBundle.cmake.in nc-config.cmake.in configure configure.ac install-sh config.h.in config.sub CTestConfig.cmake.in)

#####
# Configure and print the libnetcdf.settings file.
#####

# Set variables to mirror those used by autoconf.
# This way we don't need to maintain two separate template
# files.
SET(host_cpu "${cpu}")
SET(host_vendor "${osname}")
SET(host_os "${osrel}")
SET(abs_top_builddir "${CMAKE_BINARY_DIR}")
SET(abs_top_srcdir "${CMAKE_SOURCE_DIR}")

SET(CC_VERSION "${CMAKE_C_COMPILER}")

# Build *FLAGS for libnetcdf.settings.
SET(CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
SET(CPPFLAGS "${CMAKE_CPP_FLAGS} ${CMAKE_CPP_FLAGS_${CMAKE_BUILD_TYPE}}")
SET(LDFLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}")

is_disabled(BUILD_SHARED_LIBS enable_static)
is_enabled(BUILD_SHARED_LIBS enable_shared)

# Remove libnetcdf from NC_LIBS.
STRING(REPLACE "-lnetcdf " "" TMP_NC_LIBS "${NC_LIBS}")
SET(LIBS "${TMP_NC_LIBS}")

is_enabled(ENABLE_V2_API HAS_NC2)
is_enabled(ENABLE_NETCDF_4 HAS_NC4)
is_enabled(ENABLE_HDF4 HAS_HDF4)
is_enabled(ENABLE_NETCDF_4 HAS_HDF5)
is_enabled(USE_SZIP HAS_SZIP)
is_enabled(STATUS_PNETCDF HAS_PNETCDF)
is_enabled(STATUS_PARALLEL HAS_PARALLEL)
is_enabled(ENABLE_PARALLEL4 HAS_PARALLEL4)
is_enabled(ENABLE_DAP HAS_DAP)
is_enabled(ENABLE_DAP HAS_DAP2)
is_enabled(ENABLE_DAP4 HAS_DAP4)
is_enabled(ENABLE_DISKLESS HAS_DISKLESS)
is_enabled(USE_MMAP HAS_MMAP)
is_enabled(JNA HAS_JNA)
is_enabled(ENABLE_ZERO_LENGTH_COORD_BOUND RELAX_COORD_BOUND)
is_enabled(USE_CDF5 HAS_CDF5)
is_enabled(ENABLE_ERANGE_FILL HAS_ERANGE_FILL)

# Generate file from template.
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/libnetcdf.settings.in"
  "${CMAKE_CURRENT_BINARY_DIR}/libnetcdf.settings"
  @ONLY)

# Read in settings file, print out.
# Avoid using system-specific calls so that this
# might also work on Windows.
FILE(READ "${CMAKE_CURRENT_BINARY_DIR}/libnetcdf.settings"
  LIBNETCDF_SETTINGS)
MESSAGE(${LIBNETCDF_SETTINGS})

# Install libnetcdf.settings file into same location
# as the libraries.
INSTALL(FILES "${netCDF_BINARY_DIR}/libnetcdf.settings"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  COMPONENT libraries)

#####
# End libnetcdf.settings section.
#####

#####
# Create 'netcdf_meta.h' include file.
#####
configure_file(
  ${netCDF_SOURCE_DIR}/include/netcdf_meta.h.in
  ${netCDF_SOURCE_DIR}/include/netcdf_meta.h @ONLY)

FILE(COPY ${netCDF_SOURCE_DIR}/include/netcdf_meta.h
  DESTINATION ${netCDF_BINARY_DIR}/include/)

#####
# Build test_common.sh
#####
SET(EXTRA_DIST ${EXTRA_DIST} ${CMAKE_SOURCE_DIR}/test_common.in)
SET(TOPSRCDIR "${CMAKE_SOURCE_DIR}")
SET(TOPBUILDDIR "${CMAKE_BINARY_DIR}")
configure_file(${CMAKE_SOURCE_DIR}/test_common.in ${CMAKE_BINARY_DIR}/test_common.sh @ONLY NEWLINE_STYLE LF)

#####
# Build nc_test4/findplugin.sh
#####
SET(ISCMAKE "1")
configure_file(${CMAKE_SOURCE_DIR}/nc_test4/findplugin.in ${CMAKE_BINARY_DIR}/nc_test4/findplugin.sh @ONLY NEWLINE_STYLE LF)

IF(ENABLE_EXAMPLES)
configure_file(${CMAKE_SOURCE_DIR}/nc_test4/findplugin.in ${CMAKE_BINARY_DIR}/examples/C/findplugin.sh @ONLY NEWLINE_STYLE LF)
ENDIF()

#####
# Build ncdap_test|dap4_test/findtestserver[4].c
#####
configure_file(${CMAKE_SOURCE_DIR}/ncdap_test/findtestserver.c.in ${CMAKE_SOURCE_DIR}/ncdap_test/findtestserver.c @ONLY NEWLINE_STYLE LF)
configure_file(${CMAKE_SOURCE_DIR}/ncdap_test/findtestserver.c.in ${CMAKE_SOURCE_DIR}/dap4_test/findtestserver4.c @ONLY NEWLINE_STYLE LF)

####
# Export files
####

# Create CMake package configuration files. With these, other packages using
# cmake should be able to find netcdf using find_package and find_library.
# The EXPORT call is paired with one in liblib.
set(ConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/netCDF)

install(EXPORT netCDFTargets
  DESTINATION ${ConfigPackageLocation}
  COMPONENT headers
  )

include(CMakePackageConfigHelpers)
CONFIGURE_PACKAGE_CONFIG_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/netCDFConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/netCDFConfig.cmake"
  INSTALL_DESTINATION "${ConfigPackageLocation}"
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  PATH_VARS
  CMAKE_INSTALL_PREFIX
  CMAKE_INSTALL_INCLUDEDIR
  CMAKE_INSTALL_LIBDIR
  )

INSTALL(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/netCDFConfig.cmake"
  DESTINATION "${ConfigPackageLocation}"
  COMPONENT headers
  )

# Create export configuration
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/netCDF/netCDFConfigVersion.cmake"
  VERSION ${netCDF_VERSION}
  COMPATIBILITY SameMajorVersion
  )

install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/netCDF/netCDFConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}
  COMPONENT headers
  )

####
# End export files
####



# CPack inclusion must come last.
# INCLUDE(CPack)
INCLUDE(CMakeInstallation.cmake)
