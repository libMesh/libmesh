include $(top_builddir)/Make.common

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS) -DTETLIBRARY
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

noinst_LTLIBRARIES   = libtetgen.la
libtetgen_la_SOURCES = tetgen.C predicates.C tetgen.h
