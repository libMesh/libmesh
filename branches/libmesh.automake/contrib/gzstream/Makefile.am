# include $(top_builddir)/Make.common

AM_CPPFLAGS = $(CPPFLAGS_OPT)
AM_CFLAGS   = $(CFLAGS_OPT)
AM_CXXFLAGS = $(CXXFLAGS_OPT)

noinst_LTLIBRARIES     = libgzstream.la
libgzstream_la_SOURCES = gzstream.C gzstream.h
libgzstream_la_LIBADD  = -lz
