include $(top_builddir)/Make.common

AUTOMAKE_OPTIONS = foreign subdir-objects 
ACLOCAL_AMFLAGS  = -I m4
#SUBDIRS          = include src

# # additional files which must be included in 'make dist'
# EXTRA_DIST = \
# 	   LICENSE \
# 	   data \
# 	   dist_version

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh    \
           build-aux/missing

# Additional files to be deleted by 'make distclean'
DISTCLEANFILES =  \
	   _configs.sed

# # Additional files to be deleted by 'make maintainer-clean'
# MAINTAINERCLEANFILES =  \
# 	   $(AUX_DIST)                              \
# 	   Makefile.in                              \
# 	   aclocal.m4                               \
# 	   aminclude.am                             \
# 	   autom4te.cache/*                         \
# 	   config.h	                            \
# 	   config.h.in	                            \
# 	   config.log                               \
# 	   config.status                            \
# 	   config.sub                               \
# 	   configure                                \
# 	   data/Makefile.in                         \
# 	   src/Makefile.in                          \
# 	   src/physics/properties/Makefile.in       \
# 	   test/Makefile.in                         \
# 	   test/common/Makefile.in 		    \
# 	   test/comp_ns/Makefile.in		    \
# 	   test/unit/Makefile.in                    






SUBDIRS = src include

AM_CPPFLAGS = $(libmesh_CPPFLAGS)
AM_CFLAGS   = $(libmesh_CFLAGS)
AM_CXXFLAGS = $(libmesh_CXXFLAGS)

INCLUDES = $(libmesh_INCLUDE)

###########################################################
# building up libmesh.la - actually delegate the compilation
# of the libMesh sources to our ./src subdirectoty.
#
lib_LTLIBRARIES = libmesh.la                 # !!! the final library that we will install
libmesh_la_SOURCES =                         # !!! no local sources for this library             
libmesh_la_LIBADD  = src/libmesh_local.la    # !!! but import 



###########################################################
# add contributed libraries that are conditionally compiled

if ENABLE_GMV
  SUBDIRS += contrib/gmv
  libmesh_la_LIBADD += contrib/gmv/libgmv.la
endif

if ENABLE_GZSTREAMS
  SUBDIRS += contrib/gzstream
  libmesh_la_LIBADD += contrib/gzstream/libgzstream.la
endif

if ENABLE_NETCDF
  SUBDIRS += contrib/netcdf/Lib
  libmesh_la_LIBADD += contrib/netcdf/Lib/libnetcdf.la
endif

if ENABLE_EXODUS
  SUBDIRS += contrib/exodusii/Lib
  libmesh_la_LIBADD += contrib/exodusii/Lib/libexodusII.la
endif

if ENABLE_NEMESIS
  SUBDIRS += contrib/nemesis/Lib
  libmesh_la_LIBADD += contrib/nemesis/Lib/libnemesis.la
endif

if ENABLE_TRIANGLE
  SUBDIRS += contrib/triangle
  libmesh_la_LIBADD += contrib/triangle/libtriangle.la
endif

if ENABLE_TETGEN
  SUBDIRS += contrib/tetgen
  libmesh_la_LIBADD += contrib/tetgen/libtetgen.la
endif

if ENABLE_LIBHILBERT
  SUBDIRS += contrib/libHilbert
  libmesh_la_LIBADD += contrib/libHilbert/libHilbert.la
endif

if ENABLE_LASPACK
  SUBDIRS += contrib/laspack
  libmesh_la_LIBADD += contrib/laspack/liblaspack.la
endif

if ENABLE_METIS
  SUBDIRS += contrib/metis/Lib
  libmesh_la_LIBADD += contrib/metis/Lib/libmetis.la
endif

if ENABLE_PARMETIS
  SUBDIRS += contrib/parmetis/Lib
  libmesh_la_LIBADD += contrib/parmetis/Lib/libparmetis.la
endif

if ENABLE_SFC
  SUBDIRS += contrib/sfcurves
  libmesh_la_LIBADD += contrib/sfcurves/libsfcurves.la
endif
###########################################################





###########################################################
# Utility programs
bin_PROGRAMS = # empty, append below
bin_PROGRAMS     += meshtool
meshtool_SOURCES  = src/apps/meshtool.C
meshtool_LDADD    = libmesh.la $(libmesh_LIBS)

bin_PROGRAMS    += compare
compare_SOURCES  = src/apps/compare.C
compare_LDADD    = libmesh.la $(libmesh_LIBS)

bin_PROGRAMS     += meshbcid
meshbcid_SOURCES  = src/apps/meshbcid.C
meshbcid_LDADD    = libmesh.la $(libmesh_LIBS)

bin_PROGRAMS   += meshid
meshid_SOURCES  = src/apps/meshid.C
meshid_LDADD    = libmesh.la $(libmesh_LIBS)

bin_PROGRAMS     += meshnorm
meshnorm_SOURCES  = src/apps/meshnorm.C
meshnorm_LDADD    = libmesh.la $(libmesh_LIBS)

bin_PROGRAMS       += projection
projection_SOURCES  = src/apps/projection.C
projection_LDADD    = libmesh.la $(libmesh_LIBS)

bin_PROGRAMS += amr
amr_SOURCES   = src/apps/amr.C
amr_LDADD     = libmesh.la $(libmesh_LIBS)




#
.PHONY: doc svnexpand

#
# doxygen documentation
doc:
	mkdir -p $(top_srcdir)/doc/html
	$(DOXYGEN) $(top_builddir)/doc/Doxyfile

#
# Mark source files for SVN keyword expansion
svnexpand:
	svn -R propset svn:keywords "Date Author Revision HeadURL Id" $(top_srcdir)/include $(top_srcdir)/src


#
# special things to do when running 'make dist'
dist-hook:
	rm -rf `find $(distdir) -name .svn`


#
# Run unit tests only
check-unit:
	@$(MAKE) && cd $(top_builddir)/test/unit && $(MAKE) check

#
# Run regression tests only.  Don't forget to make compare first though!
check-regression:
	@$(MAKE) 
	@cd $(top_builddir)/test/common  && $(MAKE) check 
	@cd $(top_builddir)/test/comp_ns && $(MAKE) check




#
# 
#
.PHONY: echo echo_cxxflags echo_include echo_ldflags

#
# Useful for checking make rules
#
echo:
	@echo "C++ compiler:"
	@echo "$(libmesh_CXX)"
	@echo
	@echo "Source Files:"
	@echo "$(srcfiles)"
	@echo
	@echo "Object Files:"
	@echo "$(objects)"
	@echo
	@echo "Application object Files:"
	@echo "$(appobjects)"
	@echo
	@echo "Target:"
	@echo "$(target)"
	@echo
	@echo "Examples Source Files:"
	@echo "$(examplesrcfiles)"
	@echo
	@echo "libmesh_CPPFLAGS:"
	@echo "$(libmesh_CPPFLAGS)"
	@echo
	@echo "libmesh_CFLAGS:"
	@echo "$(libmesh_CFLAGS)"
	@echo
	@echo "libmesh_CXXFLAGS:"
	@echo "$(libmesh_CXXFLAGS)"
	@echo
	@echo "libmesh_CXXSHAREDFLAG:"
	@echo "$(libmesh_CXXSHAREDFLAG)"
	@echo
	@echo "libmesh_INCLUDE:"
	@echo "$(libmesh_INCLUDE)"
	@echo
	@echo "libmesh_LIBS:"
	@echo "$(libmesh_LIBS)"
	@echo
	@echo "libmesh_LDFLAGS:"
	@echo "$(libmesh_LDFLAGS)"
	@echo
	@echo "libmesh_DLFLAGS:"
	@echo "$(libmesh_DLFLAGS)"
	@echo
	@echo "EXAMPLES:"
	@echo "$(examplesrcfiles)"
	@echo

#
# Print the name of the C++ compiler, padded with whitespace
#
echo_cxx:
	@echo  " " $(libmesh_CXX) " "

#
# Print the name of the C compiler, padded with whitespace
#
echo_cc:
	@echo  " " $(libmesh_CC) " "

#
# Print the flags used for the C preprocessor, padded with whitespace
#
echo_cppflags:
	@echo " " $(libmesh_CPPFLAGS) " "

#
# Print the flags used for C++ compilation, padded with whitespace
#
echo_cxxflags:
	@echo " " $(libmesh_CXXFLAGS) " "

#
# Print the flags used for C compilation, padded with whitespace
#
echo_cflags:
	@echo " " $(libmesh_CFLAGS) " "

#
# Print C++ compiler include path, padded with whitespace
#
echo_include:
	@echo " " $(libmesh_INCLUDE) " "

#
# Print the flags used to link, padded with whitespace
#
echo_ldflags:
	@echo " " $(libmesh_LIBS) $(libmesh_LDFLAGS) $(libmesh_DLFLAGS) " "

