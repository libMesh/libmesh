#include $(top_builddir)/Make.common

AM_CPPFLAGS  = $(CPPFLAGS_OPT) -DTRILIBRARY -DANSI_DECLARATORS
AM_CFLAGS    = $(CFLAGS_OPT) 
AM_CXXFLAGS  = $(CXXFLAGS_OPT)
AM_CPPFLAGS += -I$(top_builddir)/include #for #include "libmesh/libmesh_config.h"

# generally we do not install contributed headers, but for triangle we
# must because it can optionally be included directly by our
# mesh_triangle_interface.h
includedir = $(prefix)/include/libmesh

noinst_LTLIBRARIES     = libtriangle.la
include_HEADERS        = triangle.h
libtriangle_la_SOURCES = triangle.c
