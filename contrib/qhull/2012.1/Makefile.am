# the -Dqh_QHpointer is required for linking the C++ API
pkg_cppflags = -I$(srcdir)/src -Dqh_QHpointer
pkg_libs     = $(QHULL_LIBS)

# libqhull c sources
libqhull_sources = \
     src/libqhull/geom.h \
     src/libqhull/io.h \
     src/libqhull/libqhull.h \
     src/libqhull/mem.h \
     src/libqhull/merge.h \
     src/libqhull/poly.h \
     src/libqhull/qhull_a.h \
     src/libqhull/qset.h \
     src/libqhull/random.h \
     src/libqhull/stat.h \
     src/libqhull/user.h \
     src/libqhull/geom.c \
     src/libqhull/geom2.c \
     src/libqhull/global.c \
     src/libqhull/io.c \
     src/libqhull/libqhull.c \
     src/libqhull/mem.c \
     src/libqhull/merge.c \
     src/libqhull/poly.c \
     src/libqhull/poly2.c \
     src/libqhull/qset.c \
     src/libqhull/random.c \
     src/libqhull/rboxlib.c \
     src/libqhull/stat.c \
     src/libqhull/user.c \
     src/libqhull/usermem.c \
     src/libqhull/userprintf.c \
     src/libqhull/userprintf_rbox.c

# libqhull c++ sources
libqhullcpp_sources = \
     src/libqhull/geom.h \
     src/libqhull/io.h \
     src/libqhull/libqhull.h \
     src/libqhull/mem.h \
     src/libqhull/merge.h \
     src/libqhull/poly.h \
     src/libqhull/qhull_a.h \
     src/libqhull/qset.h \
     src/libqhull/random.h \
     src/libqhull/stat.h \
     src/libqhull/user.h \
     src/libqhull/geom.c \
     src/libqhull/geom2.c \
     src/libqhull/global.c \
     src/libqhull/io.c \
     src/libqhull/libqhull.c \
     src/libqhull/mem.c \
     src/libqhull/merge.c \
     src/libqhull/poly.c \
     src/libqhull/poly2.c \
     src/libqhull/qset.c \
     src/libqhull/random.c \
     src/libqhull/rboxlib.c \
     src/libqhull/stat.c \
     src/libqhull/user.c \
     src/libqhull/usermem.c \
     src/libqhullcpp/Coordinates.h \
     src/libqhullcpp/functionObjects.h \
     src/libqhullcpp/PointCoordinates.h \
     src/libqhullcpp/Qhull.h \
     src/libqhullcpp/QhullError.h \
     src/libqhullcpp/QhullFacet.h \
     src/libqhullcpp/QhullFacetList.h \
     src/libqhullcpp/QhullFacetSet.h \
     src/libqhullcpp/QhullHyperplane.h \
     src/libqhullcpp/QhullIterator.h \
     src/libqhullcpp/QhullLinkedList.h \
     src/libqhullcpp/QhullPoint.h \
     src/libqhullcpp/QhullPoints.h \
     src/libqhullcpp/QhullPointSet.h \
     src/libqhullcpp/QhullQh.h \
     src/libqhullcpp/QhullRidge.h \
     src/libqhullcpp/QhullSet.h \
     src/libqhullcpp/QhullSets.h \
     src/libqhullcpp/QhullStat.h \
     src/libqhullcpp/QhullVertex.h \
     src/libqhullcpp/QhullVertexSet.h \
     src/libqhullcpp/RboxPoints.h \
     src/libqhullcpp/RoadError.h \
     src/libqhullcpp/RoadLogEvent.h \
     src/libqhullcpp/UsingLibQhull.h \
     src/libqhullcpp/Coordinates.cpp \
     src/libqhullcpp/PointCoordinates.cpp \
     src/libqhullcpp/Qhull.cpp \
     src/libqhullcpp/QhullFacet.cpp \
     src/libqhullcpp/QhullFacetList.cpp \
     src/libqhullcpp/QhullFacetSet.cpp \
     src/libqhullcpp/QhullHyperplane.cpp \
     src/libqhullcpp/QhullPoint.cpp \
     src/libqhullcpp/QhullPoints.cpp \
     src/libqhullcpp/QhullPointSet.cpp \
     src/libqhullcpp/QhullQh.cpp \
     src/libqhullcpp/QhullRidge.cpp \
     src/libqhullcpp/QhullSet.cpp \
     src/libqhullcpp/QhullStat.cpp \
     src/libqhullcpp/QhullVertex.cpp \
     src/libqhullcpp/QhullVertexSet.cpp \
     src/libqhullcpp/RboxPoints.cpp \
     src/libqhullcpp/RoadError.cpp \
     src/libqhullcpp/RoadLogEvent.cpp \
     src/libqhullcpp/UsingLibQhull.cpp

pkg_sources = $(libqhullcpp_sources)

# import as much commonality as possible
include $(top_srcdir)/contrib/Make.common






##############################################################
# What follows are the qhull executables.
# We build them as part of the 'check' targets but do not
# actually install them

check_LIBRARIES = libqhull.a libqhull_p.a libqhullcpp.a

# static C library, *without* -Dqh_QHpointer
libqhull_a_CPPFLAGS = -I$(srcdir)/src
libqhull_a_SOURCES  = $(libqhull_sources)

# static C library, *with* -Dqh_QHpointer
libqhull_p_a_CPPFLAGS = -I$(srcdir)/src -Dqh_QHpointer
libqhull_p_a_SOURCES  = $(libqhull_sources)

# static C++ library, requires -Dqh_QHpointer
libqhullcpp_a_CPPFLAGS = -I$(srcdir)/src -Dqh_QHpointer
libqhullcpp_a_SOURCES  = $(libqhullcpp_sources)


check_PROGRAMS = # empty, append below

check_PROGRAMS += rbox
rbox_SOURCES    = src/rbox/rbox.c
rbox_CPPFLAGS   = -I$(srcdir)/src/libqhull
rbox_LDADD      = libqhull.a -lm

check_PROGRAMS  += qconvex
qconvex_SOURCES  = src/qconvex/qconvex.c
qconvex_CPPFLAGS = -I$(srcdir)/src/libqhull
qconvex_LDADD    = libqhull.a -lm

check_PROGRAMS    += qdelaunay
qdelaunay_SOURCES  = src/qdelaunay/qdelaun.c
qdelaunay_CPPFLAGS = -I$(srcdir)/src/libqhull
qdelaunay_LDADD    = libqhull.a -lm

check_PROGRAMS += qhalf
qhalf_SOURCES   = src/qhalf/qhalf.c
qhalf_CPPFLAGS  = -I$(srcdir)/src/libqhull
qhalf_LDADD     = libqhull.a -lm

check_PROGRAMS   += qvoronoi
qvoronoi_SOURCES  = src/qvoronoi/qvoronoi.c
qvoronoi_CPPFLAGS = -I$(srcdir)/src/libqhull
qvoronoi_LDADD    = libqhull.a -lm

check_PROGRAMS += qhull
qhull_SOURCES   = src/qhull/unix.c
qhull_CPPFLAGS  = -I$(srcdir)/src/libqhull
qhull_LDADD     = libqhull.a -lm

check_PROGRAMS   += testqset
testqset_SOURCES  = src/testqset/testqset.c src/libqhull/qset.c src/libqhull/mem.c
testqset_CPPFLAGS = -I$(srcdir)/src/libqhull -lm

check_PROGRAMS  += user_eg
user_eg_SOURCES  = src/user_eg/user_eg.c
user_eg_CPPFLAGS = -I$(srcdir)/src/libqhull -Dqh_QHpointer
user_eg_LDADD    = libqhull_p.a -lm

check_PROGRAMS   += user_eg2
user_eg2_SOURCES  = src/user_eg2/user_eg2.c
user_eg2_CPPFLAGS = -I$(srcdir)/src/libqhull
user_eg2_LDADD    = libqhull.a -lm

check_PROGRAMS   += user_eg3
user_eg3_SOURCES  = src/user_eg3/user_eg3.cpp
user_eg3_CPPFLAGS = -I$(srcdir)/src \
                    -I$(srcdir)/src/libqhullcpp \
                    -I$(srcdir)/src/libqhull \
                    -Dqh_QHpointer
user_eg3_LDADD    = libqhullcpp.a libqhull_p.a -lm

# All these targets were pulled straight
# from the upstream Makefile
check-local:
	@echo ==============================
	@echo == Test qset.c with mem.c ====
	@echo ==============================
	-./testqset 10000
	@echo ==============================
	@echo == Run the qhull smoketest ===
	@echo ==============================
	-./rbox D4 | ./qhull Tv
	@echo ==============================
	@echo ========= qconvex ============
	@echo ==============================
	-./rbox 10 | ./qconvex Tv
	@echo ==============================
	@echo ========= qdelaunay ==========
	@echo ==============================
	-./rbox 10 | ./qdelaunay Tv
	@echo ==============================
	@echo ========= qhalf ==============
	@echo ==============================
	-./rbox 10 | ./qconvex FQ FV n Tv | ./qhalf Tv
	@echo ==============================
	@echo ========= qvoronoi ===========
	@echo ==============================
	-./rbox 10 | ./qvoronoi Tv
	@echo ==============================
	@echo ========= user_eg ============
	@echo ==============================
	-./user_eg
	@echo ==============================
	@echo ========= user_eg2 ===========
	@echo ==============================
	-./user_eg2
	@echo ==============================
	@echo ========= user_eg3 ===========
	@echo ==============================
	-./user_eg3 rbox "10 D2"  "2 D2" qhull "s p" facets
