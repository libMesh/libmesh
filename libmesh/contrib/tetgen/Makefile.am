include $(top_builddir)/Make.common.@METHOD@

AM_CPPFLAGS = $(libmesh_CPPFLAGS) -DTETLIBRARY -DLIBMESH_IS_COMPILING_TETGEN
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

# generally we do not install contributed headers, but for tetgen we
# must because it can optionally be included directly by our
# mesh_triangle_interface.h
includedir = $(prefix)/include/libmesh

noinst_LTLIBRARIES   = libtetgen.la
include_HEADERS      = tetgen.h
libtetgen_la_SOURCES = tetgen.C predicates.C
