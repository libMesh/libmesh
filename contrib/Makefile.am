AUTOMAKE_OPTIONS = subdir-objects

###########################################################
# add contributed libraries that are conditionally compiled
SUBDIRS    =
opt_libs   =
dbg_libs   =
devel_libs =
prof_libs  =
oprof_libs =

# Doesn't seem to be used directly in this Makefile; add
# it in case someone includes this Makefile?
AM_LDFLAGS   = $(libmesh_LDFLAGS)

if LIBMESH_ENABLE_CAPNPROTO
  SUBDIRS    += capnproto
  opt_libs   += capnproto/libopt.la
  dbg_libs   += capnproto/libdbg.la
  devel_libs += capnproto/libdevel.la
  prof_libs  += capnproto/libprof.la
  oprof_libs += capnproto/liboprof.la
endif

if LIBMESH_INSTALL_INTERNAL_BOOST
   SUBDIRS += boost/include
endif

if LIBMESH_INSTALL_INTERNAL_EIGEN
   SUBDIRS += eigen/eigen
endif

if LIBMESH_ENABLE_GMV
  SUBDIRS    += gmv
  opt_libs   += gmv/libopt.la
  dbg_libs   += gmv/libdbg.la
  devel_libs += gmv/libdevel.la
  prof_libs  += gmv/libprof.la
  oprof_libs += gmv/liboprof.la
endif

if LIBMESH_ENABLE_GZSTREAMS
  SUBDIRS    += gzstream
  opt_libs   += gzstream/libopt.la
  dbg_libs   += gzstream/libdbg.la
  devel_libs += gzstream/libdevel.la
  prof_libs  += gzstream/libprof.la
  oprof_libs += gzstream/liboprof.la
endif

if LIBMESH_ENABLE_NETCDF
# NetCDF version 3, build minimal subset
if LIBMESH_ENABLE_NETCDF_V3
  SUBDIRS    += netcdf/v3
  opt_libs   += netcdf/v3/libopt.la
  dbg_libs   += netcdf/v3/libdbg.la
  devel_libs += netcdf/v3/libdevel.la
  prof_libs  += netcdf/v3/libprof.la
  oprof_libs += netcdf/v3/liboprof.la
endif

# NetCDF version 4, use their autotools package
if LIBMESH_ENABLE_NETCDF_V4
  SUBDIRS    += netcdf/v4
  opt_libs   += netcdf/v4/liblib/libnetcdf.la
  dbg_libs   += netcdf/v4/liblib/libnetcdf.la
  devel_libs += netcdf/v4/liblib/libnetcdf.la
  prof_libs  += netcdf/v4/liblib/libnetcdf.la
  oprof_libs += netcdf/v4/liblib/libnetcdf.la
endif
endif


if LIBMESH_ENABLE_EXODUS
if LIBMESH_ENABLE_EXODUS_V509
  SUBDIRS    += exodusii/v5.09
  opt_libs   += exodusii/v5.09/libopt.la
  dbg_libs   += exodusii/v5.09/libdbg.la
  devel_libs += exodusii/v5.09/libdevel.la
  prof_libs  += exodusii/v5.09/libprof.la
  oprof_libs += exodusii/v5.09/liboprof.la
endif

if LIBMESH_ENABLE_EXODUS_V522
  SUBDIRS    += exodusii/v5.22/exodus
  opt_libs   += exodusii/v5.22/exodus/libopt.la
  dbg_libs   += exodusii/v5.22/exodus/libdbg.la
  devel_libs += exodusii/v5.22/exodus/libdevel.la
  prof_libs  += exodusii/v5.22/exodus/libprof.la
  oprof_libs += exodusii/v5.22/exodus/liboprof.la
endif
endif


if LIBMESH_ENABLE_NEMESIS
if LIBMESH_ENABLE_NEMESIS_V309
  SUBDIRS    += nemesis/v3.09
  opt_libs   += nemesis/v3.09/libopt.la
  dbg_libs   += nemesis/v3.09/libdbg.la
  devel_libs += nemesis/v3.09/libdevel.la
  prof_libs  += nemesis/v3.09/libprof.la
  oprof_libs += nemesis/v3.09/liboprof.la
endif

if LIBMESH_ENABLE_NEMESIS_V522
  SUBDIRS    += nemesis/v5.22/nemesis
  opt_libs   += nemesis/v5.22/nemesis/libopt.la
  dbg_libs   += nemesis/v5.22/nemesis/libdbg.la
  devel_libs += nemesis/v5.22/nemesis/libdevel.la
  prof_libs  += nemesis/v5.22/nemesis/libprof.la
  oprof_libs += nemesis/v5.22/nemesis/liboprof.la
endif
endif


if LIBMESH_ENABLE_TECIO
  SUBDIRS    += tecplot/tecio
  opt_libs   += tecplot/tecio/libopt.la
  dbg_libs   += tecplot/tecio/libdbg.la
  devel_libs += tecplot/tecio/libdevel.la
  prof_libs  += tecplot/tecio/libprof.la
  oprof_libs += tecplot/tecio/liboprof.la
endif

# I don't think this actually copies any symbols into
# libcontrib.la because we aren't using them here...
if LIBMESH_ENABLE_TECPLOT
  SUBDIRS    += tecplot/binary
  opt_libs   += tecplot/binary/libtecio.la
  dbg_libs   += tecplot/binary/libtecio.la
  devel_libs += tecplot/binary/libtecio.la
  prof_libs  += tecplot/binary/libtecio.la
  oprof_libs += tecplot/binary/libtecio.la
endif

if LIBMESH_ENABLE_TRIANGLE
  SUBDIRS    += triangle
  opt_libs   += triangle/libopt.la
  dbg_libs   += triangle/libdbg.la
  devel_libs += triangle/libdevel.la
  prof_libs  += triangle/libprof.la
  oprof_libs += triangle/liboprof.la
endif

if LIBMESH_ENABLE_TETGEN
  SUBDIRS    += tetgen
  opt_libs   += tetgen/libopt.la
  dbg_libs   += tetgen/libdbg.la
  devel_libs += tetgen/libdevel.la
  prof_libs  += tetgen/libprof.la
  oprof_libs += tetgen/liboprof.la
endif

if LIBMESH_ENABLE_QHULL
  SUBDIRS    += qhull/qhull
  opt_libs   += qhull/qhull/libopt.la
  dbg_libs   += qhull/qhull/libdbg.la
  devel_libs += qhull/qhull/libdevel.la
  prof_libs  += qhull/qhull/libprof.la
  oprof_libs += qhull/qhull/liboprof.la
endif

if LIBMESH_ENABLE_LIBHILBERT
  SUBDIRS    += libHilbert
  opt_libs   += libHilbert/libopt.la
  dbg_libs   += libHilbert/libdbg.la
  devel_libs += libHilbert/libdevel.la
  prof_libs  += libHilbert/libprof.la
  oprof_libs += libHilbert/liboprof.la
endif

if LIBMESH_ENABLE_LASPACK
  SUBDIRS    += laspack
  opt_libs   += laspack/libopt.la
  dbg_libs   += laspack/libdbg.la
  devel_libs += laspack/libdevel.la
  prof_libs  += laspack/libprof.la
  oprof_libs += laspack/liboprof.la
endif

if LIBMESH_ENABLE_METIS
if BUILD_METIS
  SUBDIRS    += metis
  opt_libs   += metis/libopt.la
  dbg_libs   += metis/libdbg.la
  devel_libs += metis/libdevel.la
  prof_libs  += metis/libprof.la
  oprof_libs += metis/liboprof.la
endif
endif

if LIBMESH_ENABLE_PARMETIS
if BUILD_METIS
  SUBDIRS    += parmetis
  opt_libs   += parmetis/libopt.la
  dbg_libs   += parmetis/libdbg.la
  devel_libs += parmetis/libdevel.la
  prof_libs  += parmetis/libprof.la
  oprof_libs += parmetis/liboprof.la
endif
endif

if LIBMESH_ENABLE_SFC
  SUBDIRS    += sfcurves
  opt_libs   += sfcurves/libopt.la
  dbg_libs   += sfcurves/libdbg.la
  devel_libs += sfcurves/libdevel.la
  prof_libs  += sfcurves/libprof.la
  oprof_libs += sfcurves/liboprof.la
endif

if LIBMESH_ENABLE_FPARSER
  # ATM always compile fparser with optimized flags,
  # and reuse in all methods
  SUBDIRS    += fparser
  opt_libs   += fparser/libopt.la
  dbg_libs   += fparser/libdbg.la
  devel_libs += fparser/libdevel.la
  prof_libs  += fparser/libprof.la
  oprof_libs += fparser/liboprof.la
endif

if LIBMESH_ENABLE_NANOFLANN
  SUBDIRS += nanoflann
  # header-only library - no library dependencies
endif

if LIBMESH_ENABLE_METAPHYSICL
  SUBDIRS += metaphysicl/0.2.0
  # header-only library for the parts we use - no library dependencies
endif



######################################################################
#
# Building the flavors
#
noinst_LTLIBRARIES =

if LIBMESH_OPT_MODE
  noinst_LTLIBRARIES        += libcontrib_opt.la
  libcontrib_opt_la_SOURCES  = cppsource.C
  libcontrib_opt_la_CPPFLAGS = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
  libcontrib_opt_la_CXXFLAGS = $(CXXFLAGS_OPT)
  libcontrib_opt_la_CFLAGS   = $(CFLAGS_OPT)
  libcontrib_opt_la_LIBADD   = $(opt_libs)
endif

if LIBMESH_DBG_MODE
  noinst_LTLIBRARIES        += libcontrib_dbg.la
  libcontrib_dbg_la_SOURCES  = cppsource.C
  libcontrib_dbg_la_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
  libcontrib_dbg_la_CXXFLAGS = $(CXXFLAGS_DBG)
  libcontrib_dbg_la_CFLAGS   = $(CFLAGS_DBG)
  libcontrib_dbg_la_LIBADD   = $(dbg_libs)
endif

if LIBMESH_DEVEL_MODE
  noinst_LTLIBRARIES          += libcontrib_devel.la
  libcontrib_devel_la_SOURCES  = cppsource.C
  libcontrib_devel_la_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
  libcontrib_devel_la_CXXFLAGS = $(CXXFLAGS_DEVEL)
  libcontrib_devel_la_CFLAGS   = $(CFLAGS_DEVEL)
  libcontrib_devel_la_LIBADD   = $(devel_libs)
endif

if LIBMESH_PROF_MODE
  noinst_LTLIBRARIES         += libcontrib_prof.la
  libcontrib_prof_la_SOURCES  = cppsource.C
  libcontrib_prof_la_CPPFLAGS = $(CPPFLAGS_PROF) $(AM_CPPFLAGS)
  libcontrib_prof_la_CXXFLAGS = $(CXXFLAGS_PROF)
  libcontrib_prof_la_CFLAGS   = $(CFLAGS_PROF)
  libcontrib_prof_la_LIBADD   = $(prof_libs)
endif

if LIBMESH_OPROF_MODE
  noinst_LTLIBRARIES          += libcontrib_oprof.la
  libcontrib_oprof_la_SOURCES  = cppsource.C
  libcontrib_oprof_la_CPPFLAGS = $(CPPFLAGS_OPROF) $(AM_CPPFLAGS)
  libcontrib_oprof_la_CXXFLAGS = $(CXXFLAGS_OPROF)
  libcontrib_oprof_la_CFLAGS   = $(CFLAGS_OPROF)
  libcontrib_oprof_la_LIBADD   = $(oprof_libs)
endif
