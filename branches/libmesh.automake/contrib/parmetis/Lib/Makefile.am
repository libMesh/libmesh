include $(top_builddir)/Make.common.@METHOD@

AM_CPPFLAGS = $(libmesh_CPPFLAGS) -DLIBMESH_IS_COMPILING_PARMETIS
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

noinst_LTLIBRARIES     = libparmetis.la
libparmetis_la_SOURCES = \
	defs.h \
	macros.h \
	parmetis.h \
	parmetislib.h \
	proto.h \
	rename.h \
	stdheaders.h \
	struct.h \
	akwayfm.c \
	ametis.c \
	backcompat.c \
	balancemylink.c \
	coarsen.c \
	comm.c \
	csrmatch.c \
	debug.c \
	diffutil.c \
	fpqueue.c \
	frename.c \
	gkmetis.c \
	grsetup.c \
	iidxsort.c \
	iintsort.c \
	ikeysort.c \
	ikeyvalsort.c \
	initbalance.c \
	initmsection.c \
	initpart.c \
	kmetis.c \
	kwaybalance.c \
	kwayfm.c \
	kwayrefine.c \
	lmatch.c \
	match.c \
	mdiffusion.c \
	memory.c \
	mesh.c \
	mmetis.c \
	move.c \
	node_refine.c \
	ometis.c \
	order.c \
	pspases.c \
	redomylink.c \
	remap.c \
	rmetis.c \
	selectq.c \
	serial.c \
	setup.c \
	stat.c \
	timer.c \
	util.c \
	wave.c \
	weird.c \
	xyzpart.c
