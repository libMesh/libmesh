include $(top_builddir)/Make.common.@METHOD@

AM_CPPFLAGS  = $(libmesh_CPPFLAGS) -DTRILIBRARY -DANSI_DECLARATORS
AM_CFLAGS    = $(libmesh_CFLAGS) 
AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CPPFLAGS += -I$(top_builddir)/include #for #include "libmesh/libmesh_config.h"

# generally we do not install contributed headers, but for triangle we
# must because it can optionally be included directly by our
# mesh_triangle_interface.h
includedir = $(prefix)/include/libmesh

noinst_LTLIBRARIES     = libtriangle.la
include_HEADERS        = triangle.h
libtriangle_la_SOURCES = triangle.c
