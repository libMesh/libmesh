include $(top_builddir)/Make.common

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS) -DTRILIBRARY -DANSI_DECLARATORS
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

noinst_LTLIBRARIES   = libtriangle.la
libtriangle_la_SOURCES = triangle.c triangle.h
