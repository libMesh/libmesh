include $(top_builddir)/Make.common.@METHOD@

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

noinst_LTLIBRARIES     = libgzstream.la
libgzstream_la_SOURCES = gzstream.C gzstream.h
libgzstream_la_LIBADD  = -lz