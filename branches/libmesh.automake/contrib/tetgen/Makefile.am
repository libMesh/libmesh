AM_CPPFLAGS = $(CPPFLAGS_OPT) -DTETLIBRARY -DLIBMESH_IS_COMPILING_TETGEN
AM_CFLAGS   = $(CFLAGS_OPT)
AM_CXXFLAGS = $(CXXFLAGS_OPT)

# generally we do not install contributed headers, but for tetgen we
# must because it can optionally be included directly by our
# mesh_triangle_interface.h
includedir = $(prefix)/include/libmesh

noinst_LTLIBRARIES   = libtetgen.la
include_HEADERS      = tetgen.h
libtetgen_la_SOURCES = tetgen.C predicates.C
