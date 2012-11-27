AM_CPPFLAGS = $(CPPFLAGS_OPT)
AM_CFLAGS   = $(CFLAGS_OPT)
AM_CXXFLAGS = $(CXXFLAGS_OPT)

noinst_LTLIBRARIES = libcontrib.la


###########################################################
# add contributed libraries that are conditionally compiled
libcontrib_la_LIBADD  =             # start empty, append below
libcontrib_la_LIBS    =
libcontrib_la_SOURCES = cppsource.C # Dummy C++ source to cause C++ linking.

SUBDIRS =
LIBS =

if LIBMESH_INSTALL_INTERNAL_BOOST
   SUBDIRS += boost/include
endif

if LIBMESH_ENABLE_GMV
  SUBDIRS += gmv
  libcontrib_la_LIBADD += gmv/libgmv.la
endif

if LIBMESH_ENABLE_GZSTREAMS
  SUBDIRS += gzstream
  libcontrib_la_LIBADD += gzstream/libgzstream.la
endif

if LIBMESH_ENABLE_NETCDF
  SUBDIRS += netcdf/Lib
  libcontrib_la_LIBADD += netcdf/Lib/libnetcdf.la
endif

if LIBMESH_ENABLE_EXODUS
  SUBDIRS += exodusii/Lib
  libcontrib_la_LIBADD += exodusii/Lib/libexodusII.la
endif

if LIBMESH_ENABLE_NEMESIS
  SUBDIRS += nemesis/Lib
  libcontrib_la_LIBADD += nemesis/Lib/libnemesis.la
endif

if LIBMESH_ENABLE_TECIO
  SUBDIRS += tecplot/tecio
  libcontrib_la_LIBADD += tecplot/tecio/libtecio.la
endif

# I don't think this actually copies any sombols into
# libcontrib.la because we aren't using them here...
if LIBMESH_ENABLE_TECPLOT
  libcontrib_la_LIBADD += @TECPLOT_LIBRARY@
endif

if LIBMESH_ENABLE_TRIANGLE
  SUBDIRS += triangle
  libcontrib_la_LIBADD += triangle/libtriangle.la
endif

if LIBMESH_ENABLE_TETGEN
  SUBDIRS += tetgen
  libcontrib_la_LIBADD += tetgen/libtetgen.la
endif

if LIBMESH_ENABLE_LIBHILBERT
  SUBDIRS += libHilbert
  libcontrib_la_LIBADD += libHilbert/libHilbert.la
endif

if LIBMESH_ENABLE_LASPACK
  SUBDIRS += laspack
  libcontrib_la_LIBADD += laspack/liblaspack.la
endif

if LIBMESH_ENABLE_METIS
  SUBDIRS += metis
  libcontrib_la_LIBADD += metis/libmetis.la
endif

if LIBMESH_ENABLE_PARMETIS
  SUBDIRS += parmetis
  libcontrib_la_LIBADD += parmetis/libparmetis.la
endif

if LIBMESH_ENABLE_SFC
  SUBDIRS += sfcurves
  libcontrib_la_LIBADD += sfcurves/libsfcurves.la
endif

if LIBMESH_ENABLE_FPARSER
  SUBDIRS += fparser
  libcontrib_la_LIBADD += fparser/libfparser.la
endif
