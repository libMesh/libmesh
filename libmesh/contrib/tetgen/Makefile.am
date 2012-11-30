# generally we do not install contributed headers, but for tetgen we
# must because it can optionally be included directly by our
# mesh_tetgen_wrapper.h
includedir = $(prefix)/include/libmesh

include_HEADERS = tetgen.h
pkg_sources     = tetgen.C predicates.C
pkg_cppflags    = -DTETLIBRARY -DLIBMESH_IS_COMPILING_TETGEN

# import as much commonality as possible
include $(top_srcdir)/contrib/Make.common
