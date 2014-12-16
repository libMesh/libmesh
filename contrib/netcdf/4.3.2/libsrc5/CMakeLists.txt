SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

SET(libsrc5_SOURCES nc5dispatch.c)

add_library(netcdf5 OBJECT ${libsrc5_SOURCES})

ADD_EXTRA_DIST(${libsrc5_SOURCES} CMakeLists.txt)
