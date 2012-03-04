include $(top_builddir)/Make.common

AM_CPPFLAGS = $(libmesh_CPPFLAGS) -DTRILIBRARY -DANSI_DECLARATORS
AM_CFLAGS   = $(libmesh_CFLAGS) 
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

# generally we do not install contributed headers, but for triangle we
# must because it can optionally be included directly by our
# mesh_triangle_interface.h
includedir = $(prefix)/include/contrib

noinst_LTLIBRARIES     = libtriangle.la
include_HEADERS        = triangle.h
libtriangle_la_SOURCES = triangle.c
