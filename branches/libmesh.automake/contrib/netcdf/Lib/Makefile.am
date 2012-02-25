include $(top_builddir)/Make.common

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

noinst_LTLIBRARIES   = libnetcdf.la
libnetcdf_la_SOURCES = attr.c ncx.c putget.c dim.c error.c libvers.c nc.c string.c v1hpg.c var.c posixio.c v2i.c \
		       config.h fbits.h nc3convert.h nc.h ncio.h ncx.h netcdf3.h netcdf.h onstack.h rnd.h

