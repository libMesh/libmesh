#include $(top_builddir)/Make.common

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS = $(CPPFLAGS_OPT) -DLIBMESH_IS_COMPILING_PARMETIS
AM_CFLAGS   = $(CFLAGS_OPT)
AM_CXXFLAGS = $(CXXFLAGS_OPT)

# include paths - note that parmetis needs metis!
AM_CPPFLAGS += -I$(top_srcdir)/contrib/metis/GKlib \
               -I$(top_srcdir)/contrib/metis/include \
               -I$(srcdir)/include \
               -I$(srcdir)/Lib 

AM_CPPFLAGS += -I$(top_builddir)/include #for #include "libmesh/libmesh_config.h"
AM_CPPFLAGS += $(libmesh_optional_INCLUDES) # for a PETSc-provided MPI, if required

noinst_LTLIBRARIES     = libparmetis.la

libparmetis_la_SOURCES = \
	include/parmetis.h \
	Lib/defs.h \
	Lib/macros.h \
	Lib/parmetislib.h \
	Lib/proto.h \
	Lib/rename.h \
	Lib/struct.h \
	Lib/akwayfm.c \
	Lib/ametis.c \
	Lib/balancemylink.c \
	Lib/comm.c \
	Lib/csrmatch.c \
	Lib/ctrl.c \
	Lib/debug.c \
	Lib/diffutil.c \
	Lib/frename.c \
	Lib/gkmetis.c \
	Lib/gkmpi.c \
	Lib/graph.c \
	Lib/initbalance.c \
	Lib/initmsection.c \
	Lib/initpart.c \
	Lib/kmetis.c \
	Lib/kwayrefine.c \
	Lib/match.c \
	Lib/mdiffusion.c \
	Lib/mesh.c \
	Lib/mmetis.c \
	Lib/move.c \
	Lib/msetup.c \
	Lib/node_refine.c \
	Lib/ometis.c \
	Lib/pspases.c \
	Lib/redomylink.c \
	Lib/remap.c \
	Lib/renumber.c \
	Lib/rmetis.c \
	Lib/selectq.c \
	Lib/serial.c \
	Lib/stat.c \
	Lib/timer.c \
	Lib/util.c \
	Lib/wave.c \
	Lib/weird.c \
	Lib/wspace.c \
	Lib/xyzpart.c
