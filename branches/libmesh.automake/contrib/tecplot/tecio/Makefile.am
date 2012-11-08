include $(top_builddir)/Make.common.@METHOD@

AM_CPPFLAGS = $(libmesh_CPPFLAGS) @TECIO_CPPFLAGS@ -DUSEENUM -DTHREED -DMAKEARCHIVE -DNO_ASSERTS
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

AUTOMAKE_OPTIONS = subdir-objects 


noinst_LTLIBRARIES = libtecio.la

libtecio_la_SOURCES = \
	include/TECIO.h \
	include/TECXXX.h \
	tecsrc/ALLOC.h \
	tecsrc/ARRLIST.h \
	tecsrc/AUXDATA.h \
	tecsrc/CHARTYPE.h \
	tecsrc/DATAIO.h \
	tecsrc/DATAIO4.h \
	tecsrc/DATASET.h \
	tecsrc/DATASET0.h \
	tecsrc/DATASHR.h \
	tecsrc/DATAUTIL.h \
	tecsrc/FACE.h \
	tecsrc/FILESTREAM.h \
	tecsrc/GEOM.h \
	tecsrc/GEOM2.h \
	tecsrc/GLOBAL.h \
	tecsrc/INPUT.h \
	tecsrc/MASTER.h \
	tecsrc/Q_MSG.h \
	tecsrc/Q_UNICODE.h \
	tecsrc/SET.h \
	tecsrc/STRLIST.h \
	tecsrc/STRUTIL.h \
	tecsrc/SYSTEM.h \
	tecsrc/TASSERT.h \
	tecsrc/TECIO.h \
	tecsrc/TECXXX.h \
	tecsrc/TEXT.h \
	tecsrc/TranslatedString.h \
	tecsrc/VERSION.h \
	tecsrc/W__BASE.h \
	tecsrc/stdafx.h \
	tecsrc/TranslatedString.cpp \
	tecsrc/alloc.cpp \
	tecsrc/arrlist.cpp \
	tecsrc/auxdata.cpp \
	tecsrc/dataio.cpp \
	tecsrc/dataio4.cpp \
	tecsrc/dataset.cpp \
	tecsrc/dataset0.cpp \
	tecsrc/datautil.cpp \
	tecsrc/filestream.cpp \
	tecsrc/geom2.cpp \
	tecsrc/q_msg.cpp \
	tecsrc/q_unicode.cpp \
	tecsrc/set.cpp \
	tecsrc/strlist.cpp \
	tecsrc/strutil.cpp \
	tecsrc/tassert.cpp \
	tecsrc/tecxxx.cpp 
