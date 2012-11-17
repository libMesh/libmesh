include $(top_builddir)/Make.common.@METHOD@

AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = $(libmesh_CPPFLAGS) -DLIBMESH_IS_COMPILING_METIS
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

AM_CPPFLAGS += -I$(srcdir)/include  -I$(srcdir)/GKlib -I$(srcdir)/Lib
AM_CPPFLAGS += -I$(top_builddir)/include #for #include "libmesh/libmesh_config.h"

noinst_LTLIBRARIES    = libmetis.la
libmetis_la_SOURCES = \
	include/metis.h \
	GKlib/GKlib.h \
	GKlib/gk_arch.h \
	GKlib/gk_defs.h \
	GKlib/gk_externs.h \
	GKlib/gk_getopt.h \
	GKlib/gk_macros.h \
	GKlib/gk_mkblas.h \
	GKlib/gk_mkmemory.h \
	GKlib/gk_mkpqueue.h \
	GKlib/gk_mkrandom.h \
	GKlib/gk_mksort.h \
	GKlib/gk_mkutils.h \
	GKlib/gk_proto.h \
	GKlib/gk_struct.h \
	GKlib/gk_types.h \
	GKlib/gkregex.h \
	GKlib/ms_inttypes.h \
	GKlib/ms_stat.h \
	GKlib/ms_stdint.h \
	GKlib/b64.c \
	GKlib/blas.c \
	GKlib/csr.c \
	GKlib/error.c \
	GKlib/evaluate.c \
	GKlib/fkvkselect.c \
	GKlib/fs.c \
	GKlib/getopt.c \
	GKlib/gkregex.c \
	GKlib/htable.c \
	GKlib/io.c \
	GKlib/itemsets.c \
	GKlib/mcore.c \
	GKlib/memory.c \
	GKlib/omp.c \
	GKlib/pdb.c \
	GKlib/pqueue.c \
	GKlib/random.c \
	GKlib/seq.c \
	GKlib/sort.c \
	GKlib/string.c \
	GKlib/timers.c \
	GKlib/tokenizer.c \
	GKlib/util.c \
	Lib/defs.h \
	Lib/gklib_defs.h \
	Lib/gklib_rename.h \
	Lib/macros.h \
	Lib/metislib.h \
	Lib/proto.h \
	Lib/rename.h \
	Lib/stdheaders.h \
	Lib/struct.h \
	Lib/auxapi.c \
	Lib/balance.c \
	Lib/bucketsort.c \
	Lib/checkgraph.c \
	Lib/coarsen.c \
	Lib/compress.c \
	Lib/contig.c \
	Lib/debug.c \
	Lib/fm.c \
	Lib/fortran.c \
	Lib/frename.c \
	Lib/gklib.c \
	Lib/graph.c \
	Lib/initpart.c \
	Lib/kmetis.c \
	Lib/kwayfm.c \
	Lib/kwayrefine.c \
	Lib/mcutil.c \
	Lib/mesh.c \
	Lib/meshpart.c \
	Lib/minconn.c \
	Lib/mincover.c \
	Lib/mmd.c \
	Lib/ometis.c \
	Lib/options.c \
	Lib/parmetis.c \
	Lib/pmetis.c \
	Lib/refine.c \
	Lib/separator.c \
	Lib/sfm.c \
	Lib/srefine.c \
	Lib/stat.c \
	Lib/timing.c \
	Lib/util.c \
	Lib/wspace.c
