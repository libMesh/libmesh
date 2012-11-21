include $(top_builddir)/Make.common.@METHOD@

bindir=@bindir@/@host@_@METHOD@
libdir=@libdir@/@host@_@METHOD@

SUBDIRS = include #src

AUTOMAKE_OPTIONS = foreign subdir-objects 
ACLOCAL_AMFLAGS  = -I m4
AM_CPPFLAGS      = $(libmesh_CPPFLAGS) -DLIBMESH_IS_COMPILING_ITSELF
AM_CFLAGS        = $(libmesh_CFLAGS)
AM_CXXFLAGS      = $(libmesh_CXXFLAGS)
AM_CPPFLAGS     += $(libmesh_contrib_INCLUDES) \
                   $(libmesh_INCLUDE) \
                   -I$(top_builddir)/include # required for libmesh_version.h

LIBS             = $(libmesh_LIBS)

# additional files which must be included in 'make dist'
EXTRA_DIST = license \
             reference_elements \
             include/include_HEADERS \
	     include/rebuild_include_HEADERS.sh \
	     src/libmesh_la_SOURCES \
	     src/rebuild_libmesh_la_SOURCES.sh \
	     contrib/tecplot/include \
             contrib/tecplot/lib \
             contrib/tecplot/README \
	     contrib/bin/strip_dup_incl_paths.pl \
	     contrib/bin/strip_dup_libs.pl

nobase_data_DATA = \
	reference_elements/2D/one_quad.xda \
	reference_elements/2D/one_quad8.xda \
	reference_elements/2D/one_quad9.xda \
	reference_elements/2D/one_tri.xda \
	reference_elements/2D/one_tri6.xda \
	reference_elements/3D/one_hex.xda \
	reference_elements/3D/one_hex20.xda \
	reference_elements/3D/one_hex27.xda \
	reference_elements/3D/one_prism.xda \
	reference_elements/3D/one_prism15.xda \
	reference_elements/3D/one_prism18.xda \
	reference_elements/3D/one_pyramid.xda \
	reference_elements/3D/one_tet.xda \
	reference_elements/3D/one_tet10.xda

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh \
           build-aux/missing

# Support for pkgconfig
pkgconfigdir     = $(libdir)/pkgconfig
pkgconfig_DATA   = Make.common.@METHOD@ contrib/utils/libmesh.pc
configdir        = $(sysconfdir)/libmesh
config_DATA      = contrib/utils/Make.common contrib/utils/libmesh.pc

# install the old make_dependencies.pl script for backwards compatibility
contribbindir      = $(prefix)/contrib/bin
contribbin_SCRIPTS = contrib/bin/make_dependencies.pl
EXTRA_DIST        += contrib/bin/make_dependencies.pl

# when doing 'make install, symlik etc/libmesh/Make.common
install-data-hook:
	cd $(DESTDIR)$(prefix) && rm -f Make.common && $(LN_S) etc/libmesh/Make.common Make.common

# Additional files to be deleted by 'make distclean'
# _configs.sed seems to be leftover by prefix_config
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







# A convenience library to hold proper libMesh 
# objects.  This will get appended with any contributed 
# sources to create the final library.
lib_LTLIBRARIES = libmesh.la

libmesh_la_LIBADD = # start empty, append below

# definition of libmesh_la_SOURCES - get from auto-maintained list
include src/libmesh_la_SOURCES


# ###########################################################
# # building up libmesh.la - actually delegate the compilation
# # of the libMesh sources to our ./src subdirectoty.
# #
# # see file:///Users/benkirk/Desktop/automake.html#Libtool-Convenience-Libraries 
# #     file:///Users/benkirk/Desktop/automake.html#How-the-Linker-is-Chosen
# # without sources we may be using the c linker?
# lib_LTLIBRARIES    = libmesh.la              # !!! the final library that we will install
# libmesh_la_SOURCES =                         # !!! no local sources for this library             
# libmesh_la_LIBADD  = src/libmesh_local.la    # !!! but import from our sources
# libmesh_la_LINK    = $(CXXLINK)              # manually specify the linker because there are no 
#                                              # listed sources to infer it from

# # This rule will be run after install
# install-exec-hook:
# 	cd $(DESTDIR)$(bindir) && \
#  	  $(MKDIR_P) $(hosttype)_$(METHOD) && \
# 	  mv -f meshtool $(hosttype)_$(METHOD) && \
# 	  $(LN_S) $(hosttype)_$(METHOD)/meshtool meshtool

###########################################################
# add contributed libraries that are conditionally compiled

if LIBMESH_INSTALL_INTERNAL_BOOST
   SUBDIRS += contrib/boost/include
endif

if LIBMESH_ENABLE_GMV
  SUBDIRS += contrib/gmv
  libmesh_la_LIBADD += contrib/gmv/libgmv.la
endif

if LIBMESH_ENABLE_GZSTREAMS
  SUBDIRS += contrib/gzstream
  libmesh_la_LIBADD += contrib/gzstream/libgzstream.la
endif

if LIBMESH_ENABLE_NETCDF
  SUBDIRS += contrib/netcdf/Lib
  libmesh_la_LIBADD += contrib/netcdf/Lib/libnetcdf.la
endif

if LIBMESH_ENABLE_EXODUS
  SUBDIRS += contrib/exodusii/Lib
  libmesh_la_LIBADD += contrib/exodusii/Lib/libexodusII.la
endif

if LIBMESH_ENABLE_NEMESIS
  SUBDIRS += contrib/nemesis/Lib
  libmesh_la_LIBADD += contrib/nemesis/Lib/libnemesis.la
endif

if LIBMESH_ENABLE_TECIO
  SUBDIRS += contrib/tecplot/tecio
  libmesh_la_LIBADD += contrib/tecplot/tecio/libtecio.la
endif

if LIBMESH_ENABLE_TECPLOT
  libmesh_la_LIBADD += @TECPLOT_LIBRARY@
endif

if LIBMESH_ENABLE_TRIANGLE
  SUBDIRS += contrib/triangle
  libmesh_la_LIBADD += contrib/triangle/libtriangle.la
endif

if LIBMESH_ENABLE_TETGEN
  SUBDIRS += contrib/tetgen
  libmesh_la_LIBADD += contrib/tetgen/libtetgen.la
endif

if LIBMESH_ENABLE_LIBHILBERT
  SUBDIRS += contrib/libHilbert
  libmesh_la_LIBADD += contrib/libHilbert/libHilbert.la
endif

if LIBMESH_ENABLE_LASPACK
  SUBDIRS += contrib/laspack
  libmesh_la_LIBADD += contrib/laspack/liblaspack.la
endif

if LIBMESH_ENABLE_METIS
  SUBDIRS += contrib/metis
  libmesh_la_LIBADD += contrib/metis/libmetis.la
endif

if LIBMESH_ENABLE_PARMETIS
  SUBDIRS += contrib/parmetis
  libmesh_la_LIBADD += contrib/parmetis/libparmetis.la
endif

if LIBMESH_ENABLE_SFC
  SUBDIRS += contrib/sfcurves
  libmesh_la_LIBADD += contrib/sfcurves/libsfcurves.la
endif

if LIBMESH_ENABLE_FPARSER
  SUBDIRS += contrib/fparser
  libmesh_la_LIBADD += contrib/fparser/libfparser.la
endif
###########################################################



###########################################################
# Utility programs
bin_PROGRAMS = # empty, append below
bin_SCRIPTS  = # empty, append below
CLEANFILES   = # empty, append below

bin_PROGRAMS     += meshtool
meshtool_SOURCES  = src/apps/meshtool.C
meshtool_LDADD    = libmesh.la

bin_PROGRAMS    += compare
compare_SOURCES  = src/apps/compare.C
compare_LDADD    = libmesh.la

bin_PROGRAMS     += meshbcid
meshbcid_SOURCES  = src/apps/meshbcid.C
meshbcid_LDADD    = libmesh.la

bin_PROGRAMS   += meshid
meshid_SOURCES  = src/apps/meshid.C
meshid_LDADD    = libmesh.la

bin_PROGRAMS     += meshnorm
meshnorm_SOURCES  = src/apps/meshnorm.C
meshnorm_LDADD    = libmesh.la

bin_PROGRAMS       += projection
projection_SOURCES  = src/apps/projection.C
projection_LDADD    = libmesh.la

bin_PROGRAMS += amr
amr_SOURCES   = src/apps/amr.C
amr_LDADD     = libmesh.la

bin_PROGRAMS                  += output_libmesh_version
output_libmesh_version_SOURCES = src/apps/output_libmesh_version.C
output_libmesh_version_LDADD   = libmesh.la

bin_PROGRAMS     += meshplot
meshplot_SOURCES  = src/apps/meshplot.C
meshplot_LDADD    = libmesh.la

bin_PROGRAMS               += solution_components
solution_components_SOURCES = src/apps/solution_components.C
solution_components_LDADD   = libmesh.la

# create a METHOD-specific configure script - will need to clean it too!
bin_SCRIPTS += contrib/bin/@host@_@METHOD@/libmesh-config
CLEANFILES  += contrib/bin/@host@_@METHOD@/libmesh-config

contrib/bin/@host@_@METHOD@/libmesh-config : contrib/bin/libmesh-config-METHOD
	@echo "copying $< to $@"
	@$(MKDIR_P) contrib/bin/@host@_@METHOD@
	@cp $< $@

# Support for libmesh-config-METHOD in $(exec_prefix)/bin
libmesh_configdir      = $(exec_prefix)/bin
libmesh_config_SCRIPTS = contrib/bin/libmesh-config

###########################################################
# Examples
SUBDIRS += examples


#
.PHONY: examples_doc doc svnexpand

#
# documentation for the examples
examples_doc:
	@cd $(top_builddir)/examples && $(MAKE) doc

#
# doxygen documentation - be sure to generate symlinks in include/libmesh!
doc html:
	$(MAKE) examples_doc
	@cd $(top_builddir)/include/libmesh && $(MAKE)
	$(MKDIR_P) $(top_builddir)/doc/html
	$(DOXYGEN) $(top_builddir)/doc/Doxyfile

# Upload the web page to sourceforge.  We need a way to specify usernames
# other than $USER when connecting to sourceforge servers... Please set the
# environment variable: $LIBMESH_SVN_USER.
upload:
	chmod -R g+w $(top_builddir)/doc/html/* $(top_srcdir)/doc/html/* $(top_srcdir)/doc/latex/*/*
	rsync -rltzve ssh --exclude '.svn' $(top_srcdir)/doc/html/ $(LIBMESH_SVN_USER),libmesh@web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -rltzve ssh --exclude '.svn' $(top_builddir)/doc/html/ $(LIBMESH_SVN_USER),libmesh@web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -rltzve ssh --exclude '.svn' $(top_srcdir)/doc/latex/howto $(LIBMESH_SVN_USER),libmesh@web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -rltzve ssh --exclude '.svn' $(top_srcdir)/doc/latex/xda_format $(LIBMESH_SVN_USER),libmesh@web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	chmod -R g-w $(top_builddir)/doc/html/* $(top_srcdir)/doc/html/* $(top_srcdir)/doc/latex/*/*

upload_test:
	chmod -R g+w $(top_builddir)/doc/html/* $(top_srcdir)/doc/html/* $(top_srcdir)/doc/latex/*/*
	rsync -nrltzve ssh --exclude '.svn' $(top_srcdir)/doc/html/ $(LIBMESH_SVN_USER),libmesh@web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -nrltzve ssh --exclude '.svn' $(top_builddir)/doc/html/ $(LIBMESH_SVN_USER),libmesh@web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -nrltzve ssh --exclude '.svn' $(top_srcdir)/doc/latex/howto $(LIBMESH_SVN_USER),libmesh@web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	rsync -nrltzve ssh --exclude '.svn' $(top_srcdir)/doc/latex/xda_format $(LIBMESH_SVN_USER),libmesh@web.sourceforge.net:/home/groups/l/li/libmesh/htdocs/
	chmod -R g-w $(top_builddir)/doc/html/* $(top_srcdir)/doc/html/* $(top_srcdir)/doc/latex/*/*

# #
# # Mark source files for SVN keyword expansion
# svnexpand:
# 	svn -R propset svn:keywords "Date Author Revision HeadURL Id" $(top_srcdir)/include $(top_srcdir)/src


#
# special things to do when running 'make dist'
dist-hook:
	rm -rf `find $(distdir) -name .svn`


#
# support 'make run_examples'
run_examples:
	@$(MAKE) && cd $(top_builddir)/examples && $(MAKE) check
