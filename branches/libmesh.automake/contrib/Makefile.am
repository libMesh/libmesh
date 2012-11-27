AUTOMAKE_OPTIONS = subdir-objects 

###########################################################
# add contributed libraries that are conditionally compiled
SUBDIRS    =
opt_libs   =
dbg_libs   =
devel_libs =

if LIBMESH_INSTALL_INTERNAL_BOOST
   SUBDIRS += boost/include
endif


if LIBMESH_ENABLE_GMV
  SUBDIRS    += gmv
  opt_libs   += gmv/libopt.la 
  dbg_libs   += gmv/libdbg.la 
  devel_libs += gmv/libdevel.la 
endif

if LIBMESH_ENABLE_GZSTREAMS
  SUBDIRS    += gzstream
  opt_libs   += gzstream/libopt.la
  dbg_libs   += gzstream/libdbg.la
  devel_libs += gzstream/libdevel.la
endif

if LIBMESH_ENABLE_NETCDF
  SUBDIRS    += netcdf/Lib
  opt_libs   += netcdf/Lib/libopt.la
  dbg_libs   += netcdf/Lib/libdbg.la
  devel_libs += netcdf/Lib/libdevel.la
endif

if LIBMESH_ENABLE_EXODUS
  SUBDIRS    += exodusii/Lib
  opt_libs   += exodusii/Lib/libopt.la
  dbg_libs   += exodusii/Lib/libdbg.la
  devel_libs += exodusii/Lib/libdevel.la
endif

if LIBMESH_ENABLE_NEMESIS
  SUBDIRS    += nemesis/Lib
  opt_libs   += nemesis/Lib/libopt.la 
  dbg_libs   += nemesis/Lib/libdbg.la 
  devel_libs += nemesis/Lib/libdevel.la 
endif

if LIBMESH_ENABLE_TECIO
  SUBDIRS    += tecplot/tecio
  opt_libs   += tecplot/tecio/libopt.la
  dbg_libs   += tecplot/tecio/libdbg.la
  devel_libs += tecplot/tecio/libdevel.la
endif

# I don't think this actually copies any symbols into
# libcontrib.la because we aren't using them here...
if LIBMESH_ENABLE_TECPLOT
  opt_libs   += @TECPLOT_LIBRARY@
  dbg_libs   += @TECPLOT_LIBRARY@
  devel_libs += @TECPLOT_LIBRARY@
endif

if LIBMESH_ENABLE_TRIANGLE
  SUBDIRS    += triangle
  opt_libs   += triangle/libopt.la
  dbg_libs   += triangle/libdbg.la
  devel_libs += triangle/libdevel.la
endif

if LIBMESH_ENABLE_TETGEN
  SUBDIRS    += tetgen
  opt_libs   += tetgen/libopt.la
  dbg_libs   += tetgen/libdbg.la
  devel_libs += tetgen/libdevel.la
endif

if LIBMESH_ENABLE_LIBHILBERT
  SUBDIRS    += libHilbert
  opt_libs   += libHilbert/libopt.la
  dbg_libs   += libHilbert/libdbg.la
  devel_libs += libHilbert/libdevel.la
endif

if LIBMESH_ENABLE_LASPACK
  SUBDIRS    += laspack
  opt_libs   += laspack/libopt.la
  dbg_libs   += laspack/libdbg.la
  devel_libs += laspack/libdevel.la
endif

if LIBMESH_ENABLE_METIS
  SUBDIRS    += metis
  opt_libs   += metis/libopt.la
  dbg_libs   += metis/libdbg.la
  devel_libs += metis/libdevel.la
endif

if LIBMESH_ENABLE_PARMETIS
  SUBDIRS    += parmetis
  opt_libs   += parmetis/libopt.la
  dbg_libs   += parmetis/libdbg.la
  devel_libs += parmetis/libdevel.la
endif

if LIBMESH_ENABLE_SFC
  SUBDIRS    += sfcurves
  opt_libs   += sfcurves/libopt.la
  dbg_libs   += sfcurves/libdbg.la
  devel_libs += sfcurves/libdevel.la
endif

if LIBMESH_ENABLE_FPARSER
  # ATM always compile fparser with optimized flags,
  # and reuse in all methods
  SUBDIRS    += fparser
  opt_libs   += fparser/libfparser.la
  dbg_libs   += fparser/libfparser.la
  devel_libs += fparser/libfparser.la
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
