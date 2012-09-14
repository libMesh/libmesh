include $(top_builddir)/Make.common.@METHOD@

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

noinst_LTLIBRARIES   = libnetcdf.la

# generally we do not install contributed headers, but Moose and
# some client apps expect exodusII.h to be there when libMesh
# is present. since we've implicitly guaranteed this for ~10 years
# no need to break that trust now.  Since exodusII.h includes netcdf.h,
# here we are...
includedir      = $(prefix)/include/libmesh
include_HEADERS = netcdf.h

# the files required to actually build libnetcdf.la
libnetcdf_la_SOURCES = attr.c ncx.c putget.c dim.c error.c libvers.c nc.c string.c v1hpg.c var.c posixio.c v2i.c \
		       config.h fbits.h nc3convert.h nc.h ncio.h ncx.h netcdf3.h netcdf.h onstack.h rnd.h

