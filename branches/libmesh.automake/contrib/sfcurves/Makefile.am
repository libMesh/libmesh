include $(top_builddir)/Make.common.@METHOD@

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

noinst_LTLIBRARIES     = libsfcurves.la
libsfcurves_la_SOURCES = sfcurves.h sfcurves_internal.h cmp.c hilbert.c morton.c
