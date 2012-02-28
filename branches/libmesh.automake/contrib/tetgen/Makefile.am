include $(top_builddir)/Make.common

AM_CPPFLAGS = $(libmesh_CPPFLAGS) -DTETLIBRARY -DLIBMESH_IS_COMPILING_TETGEN
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

includedir = $(prefix)/include/contrib

noinst_LTLIBRARIES   = libtetgen.la
include_HEADERS      = tetgen.h
libtetgen_la_SOURCES = tetgen.C predicates.C
