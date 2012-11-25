include $(top_builddir)/Make.common.@METHOD@

AM_CPPFLAGS  = $(libmesh_CPPFLAGS)
AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CPPFLAGS += -I$(top_builddir)/include #for #include "libmesh/libmesh_config.h"

noinst_LTLIBRARIES = liblaspack.la

liblaspack_la_SOURCES  = \
	copyrght.h \
	eigenval.h \
	elcmp.h \
	errhandl.h \
	factor.h \
	itersolv.h \
	laspack_config.h \
	lastypes.h \
	matrix.h \
	mlsolv.h \
	operats.h \
	precond.h \
	qmatrix.h \
	qvector.h \
	rtc.h \
	version.h \
	eigenval.C \
	errhandl.C \
	factor.C \
	itersolv.C \
	matrix.C \
	mlsolv.C \
	operats.C \
	precond.C \
	qmatrix.C \
	qvector.C \
	rtc.C
