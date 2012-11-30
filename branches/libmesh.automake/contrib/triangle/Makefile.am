# generally we do not install contributed headers, but for triangle we
# must because it can optionally be included directly by our
# mesh_triangle_interface.h
includedir = $(prefix)/include/libmesh

include_HEADERS = triangle.h
pkg_sources     = triangle.c
pkg_cppflags    = -DTRILIBRARY -DANSI_DECLARATORS \
                  -I$(top_builddir)/include #for #include "libmesh/libmesh_config.h"

# import as much commonality as possible
include $(top_srcdir)/contrib/Make.common
