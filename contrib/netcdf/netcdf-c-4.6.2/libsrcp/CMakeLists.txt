SET(libsrcp_SOURCES ncpdispatch.c)

add_library(netcdfp OBJECT ${libsrcp_SOURCES})

ADD_EXTRA_DIST(${libsrcp_SOURCES} CMakeLists.txt)
