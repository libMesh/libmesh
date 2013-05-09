# generally we do not install contributed headers, but Moose and
# some client apps expect exodusII.h to be there when libMesh
# is present. since we've implicitly guaranteed this for ~10 years
# no need to break that trust now.
includedir      = $(prefix)/include/libmesh
include_HEADERS =

pkg_cppflags = -I$(srcdir)/src -Dqh_QHpointer

pkg_sources = # empty, append below

# libqhull c sources
pkg_sources  += \
     src/libqhull/geom.c \
     src/libqhull/geom.h \
     src/libqhull/geom2.c \
     src/libqhull/global.c \
     src/libqhull/index.htm \
     src/libqhull/io.c \
     src/libqhull/io.h \
     src/libqhull/libqhull.c \
     src/libqhull/libqhull.h \
     src/libqhull/mem.c \
     src/libqhull/mem.h \
     src/libqhull/merge.c \
     src/libqhull/merge.h \
     src/libqhull/poly.c \
     src/libqhull/poly.h \
     src/libqhull/poly2.c \
     src/libqhull/qhull_a.h \
     src/libqhull/qset.c \
     src/libqhull/qset.h \
     src/libqhull/random.c \
     src/libqhull/random.h \
     src/libqhull/rboxlib.c \
     src/libqhull/stat.c \
     src/libqhull/stat.h \
     src/libqhull/user.c \
     src/libqhull/user.h \
     src/libqhull/usermem.c \
     src/libqhull/userprintf.c \
     src/libqhull/userprintf_rbox.c

# libqhull c++ sources
pkg_sources  += \
     src/libqhullcpp/Coordinates.cpp \
     src/libqhullcpp/Coordinates.h \
     src/libqhullcpp/functionObjects.h \
     src/libqhullcpp/PointCoordinates.cpp \
     src/libqhullcpp/PointCoordinates.h \
     src/libqhullcpp/Qhull.cpp \
     src/libqhullcpp/Qhull.h \
     src/libqhullcpp/QhullError.h \
     src/libqhullcpp/QhullFacet.cpp \
     src/libqhullcpp/QhullFacet.h \
     src/libqhullcpp/QhullFacetList.cpp \
     src/libqhullcpp/QhullFacetList.h \
     src/libqhullcpp/QhullFacetSet.cpp \
     src/libqhullcpp/QhullFacetSet.h \
     src/libqhullcpp/QhullHyperplane.cpp \
     src/libqhullcpp/QhullHyperplane.h \
     src/libqhullcpp/QhullIterator.h \
     src/libqhullcpp/QhullLinkedList.h \
     src/libqhullcpp/QhullPoint.cpp \
     src/libqhullcpp/QhullPoint.h \
     src/libqhullcpp/QhullPoints.cpp \
     src/libqhullcpp/QhullPoints.h \
     src/libqhullcpp/QhullPointSet.cpp \
     src/libqhullcpp/QhullPointSet.h \
     src/libqhullcpp/QhullQh.cpp \
     src/libqhullcpp/QhullQh.h \
     src/libqhullcpp/QhullRidge.cpp \
     src/libqhullcpp/QhullRidge.h \
     src/libqhullcpp/QhullSet.cpp \
     src/libqhullcpp/QhullSet.h \
     src/libqhullcpp/QhullSets.h \
     src/libqhullcpp/QhullStat.cpp \
     src/libqhullcpp/QhullStat.h \
     src/libqhullcpp/QhullVertex.cpp \
     src/libqhullcpp/QhullVertex.h \
     src/libqhullcpp/QhullVertexSet.cpp \
     src/libqhullcpp/QhullVertexSet.h \
     src/libqhullcpp/RboxPoints.cpp \
     src/libqhullcpp/RboxPoints.h \
     src/libqhullcpp/RoadError.cpp \
     src/libqhullcpp/RoadError.h \
     src/libqhullcpp/RoadLogEvent.cpp \
     src/libqhullcpp/RoadLogEvent.h \
     src/libqhullcpp/UsingLibQhull.cpp \
     src/libqhullcpp/UsingLibQhull.h

# import as much commonality as possible
include $(top_srcdir)/contrib/Make.common
