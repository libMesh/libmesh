include $(top_builddir)/Make.common.@METHOD@

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS += -I$(srcdir)/include

noinst_LTLIBRARIES    = libHilbert.la

libHilbert_la_SOURCES = \
	src/BigBitVec.C \
	src/FixBitVec.C \
	src/Hilbert.C \
	src/hilbert_mpi_ops.C \
	include/Hilbert/GetLocation.hpp \
	include/Hilbert/SetBits.hpp \
	include/Hilbert/Common.hpp \
	include/Hilbert/BigBitVec.hpp \
	include/Hilbert/FixBitVec.hpp \
	include/Hilbert/Algorithm.hpp \
	include/Hilbert/GetBits.hpp \
	include/Hilbert/GrayCodeRank.hpp \
	include/Hilbert/SetLocation.hpp \
	include/Hilbert.hpp \
	include/hilbert.h 
