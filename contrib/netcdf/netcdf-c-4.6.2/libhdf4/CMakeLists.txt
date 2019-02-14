## This is a CMake file, part of Unidata's netCDF package.
# Copyright 2018, see the COPYRIGHT file for more information.
#
# This builds the HDF4 dispatch layer.
#
# Ed Hartnett

# The source files for the HDF4 dispatch layer.
SET(libhdf4_SOURCES hdf4dispatch.c hdf4file.c hdf4func.c hdf4var.c)

# Build the HDF4 dispatch layer as a library that will be included in
# the netCDF library.
add_library(netcdfhdf4 OBJECT ${libhdf4_SOURCES})

# Remember to package this file for CMake builds.
ADD_EXTRA_DIST(${libhdf4_SOURCES} CMakeLists.txt)
