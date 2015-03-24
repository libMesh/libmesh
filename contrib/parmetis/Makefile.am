# include paths - note that parmetis needs metis!
pkg_cppflags  = -DLIBMESH_IS_COMPILING_PARMETIS
pkg_cppflags += -I$(top_srcdir)/contrib/metis/GKlib \
                -I$(top_srcdir)/contrib/metis/include \
                -I$(srcdir)/include \
                -I$(srcdir)/libparmetis \
                -I$(top_srcdir)/contrib/metis/libmetis
pkg_cppflags += -I$(top_builddir)/include #for #include "libmesh/libmesh_config.h"
pkg_cppflags += $(libmesh_optional_INCLUDES) # for a PETSc-provided MPI, if required
pkg_sources   = \
	include/parmetis.h \
	libparmetis/akwayfm.c \
	libparmetis/ametis.c \
	libparmetis/balancemylink.c \
	libparmetis/comm.c \
	libparmetis/csrmatch.c \
	libparmetis/ctrl.c \
	libparmetis/debug.c \
	libparmetis/defs.h \
	libparmetis/diffutil.c \
	libparmetis/frename.c \
	libparmetis/gkmetis.c \
	libparmetis/gkmpi.c \
	libparmetis/graph.c \
	libparmetis/initbalance.c \
	libparmetis/initmsection.c \
	libparmetis/initpart.c \
	libparmetis/kmetis.c \
	libparmetis/kwayrefine.c \
	libparmetis/macros.h \
	libparmetis/match.c \
	libparmetis/mdiffusion.c \
	libparmetis/mesh.c \
	libparmetis/mmetis.c \
	libparmetis/move.c \
	libparmetis/msetup.c \
	libparmetis/node_refine.c \
	libparmetis/ometis.c \
	libparmetis/parmetislib.h \
	libparmetis/proto.h \
	libparmetis/pspases.c \
	libparmetis/redomylink.c \
	libparmetis/remap.c \
	libparmetis/rename.h \
	libparmetis/renumber.c \
	libparmetis/rmetis.c \
	libparmetis/selectq.c \
	libparmetis/serial.c \
	libparmetis/stat.c \
	libparmetis/struct.h \
	libparmetis/timer.c \
	libparmetis/util.c \
	libparmetis/wave.c \
	libparmetis/weird.c \
	libparmetis/wspace.c \
	libparmetis/xyzpart.c

# import as much commonality as possible
include $(top_srcdir)/contrib/Make.common
