include $(top_builddir)/Make.common

#bindir=@bindir@/@host@_@METHOD@
#libdir=@libdir@/@host@_@METHOD@

SUBDIRS = include #src

AUTOMAKE_OPTIONS = foreign subdir-objects 
ACLOCAL_AMFLAGS  = -I m4
AM_CFLAGS        = $(libmesh_CFLAGS)

# AM_CPPFLAGS are method-independent cppflags that 
# we use when compiling libmesh proper, or its utility
# codes.  These will be augmented with method-specific
# cppflags later on
AM_CPPFLAGS      = -DLIBMESH_IS_COMPILING_ITSELF \
		   $(libmesh_contrib_INCLUDES) \
                   $(libmesh_INCLUDE) \
                   -I$(top_builddir)/include # required for libmesh_version.h

LIBS             = $(libmesh_LIBS)

if LIBMESH_ENABLE_TECPLOT
  LIBS += @TECPLOT_LIBRARY@
endif

# additional files which must be included in 'make dist'
EXTRA_DIST = license \
             reference_elements \
             include/include_HEADERS \
	     include/rebuild_include_HEADERS.sh \
	     src/libmesh_SOURCES \
	     src/rebuild_libmesh_SOURCES.sh \
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
pkgconfig_DATA   = Make.common \
                   contrib/utils/libmesh.pc \
                   contrib/utils/libmesh-opt.pc \
                   contrib/utils/libmesh-dbg.pc \
                   contrib/utils/libmesh-devel.pc

configdir        = $(sysconfdir)/libmesh
config_DATA      = contrib/utils/Make.common \
                   contrib/utils/libmesh.pc \
                   contrib/utils/libmesh-opt.pc \
                   contrib/utils/libmesh-dbg.pc \
                   contrib/utils/libmesh-devel.pc

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



SUBDIRS += contrib

include src/libmesh_SOURCES

# A convenience library to hold proper libMesh 
# objects.  This will get appended with any contributed 
# sources to create the final library.
lib_LTLIBRARIES = 

if LIBMESH_DBG_MODE
  lib_LTLIBRARIES            += dbg/libmesh_dbg.la
  dbg_libmesh_dbg_la_SOURCES  = $(libmesh_SOURCES)
  dbg_libmesh_dbg_la_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS) 
  dbg_libmesh_dbg_la_CXXFLAGS = $(CXXFLAGS_DBG)
  dbg_libmesh_dbg_la_CFLAGS   = $(CFLAGS_DBG)
  dbg_libmesh_dbg_la_LIBADD   = contrib/libcontrib.la $(LIBS)
endif

if LIBMESH_DEVEL_MODE
  lib_LTLIBRARIES                += devel/libmesh_devel.la
  devel_libmesh_devel_la_SOURCES  = $(libmesh_SOURCES)
  devel_libmesh_devel_la_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS) 
  devel_libmesh_devel_la_CXXFLAGS = $(CXXFLAGS_DEVEL)
  devel_libmesh_devel_la_CFLAGS   = $(CFLAGS_DEVEL)
  devel_libmesh_devel_la_LIBADD   = contrib/libcontrib.la $(LIBS)
endif

if LIBMESH_OPT_MODE
  lib_LTLIBRARIES        += opt/libmesh.la
  opt_libmesh_la_SOURCES  = $(libmesh_SOURCES)
  opt_libmesh_la_CPPFLAGS = $(CPPFLAGS_OPT) $(AM_CPPFLAGS) 
  opt_libmesh_la_CXXFLAGS = $(CXXFLAGS_OPT)
  opt_libmesh_la_CFLAGS   = $(CFLAGS_OPT)
  opt_libmesh_la_LIBADD   = contrib/libcontrib.la $(LIBS)
endif


if LIBMESH_ENABLE_CPPUNIT
 SUBDIRS += tests
endif

###########################################################


bin_PROGRAMS = # empty, append below
bin_SCRIPTS  = # empty, append below
CLEANFILES   = # empty, append below

###########################################################
# Utility programs
opt_programs   = # empty, append below
devel_programs = # empty, append below
dbg_programs   = # empty, append below

# meshtool
opt_programs           += opt/meshtool
opt_meshtool_SOURCES    = src/apps/meshtool.C
opt_meshtool_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
opt_meshtool_CXXFLAGS   = $(CXXFLAGS_OPT)
opt_meshtool_LDADD      = opt/libmesh.la

devel_programs         += devel/meshtool-devel
devel_meshtool_devel_SOURCES  = src/apps/meshtool.C
devel_meshtool_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
devel_meshtool_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
devel_meshtool_devel_LDADD    = devel/libmesh_devel.la

dbg_programs           += dbg/meshtool-dbg
dbg_meshtool_dbg_SOURCES    = src/apps/meshtool.C
dbg_meshtool_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
dbg_meshtool_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
dbg_meshtool_dbg_LDADD      = dbg/libmesh_dbg.la

# compare
opt_programs          += opt/compare
opt_compare_SOURCES    = src/apps/compare.C
opt_compare_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
opt_compare_CXXFLAGS   = $(CXXFLAGS_OPT)
opt_compare_LDADD      = opt/libmesh.la

devel_programs        += devel/compare-devel
devel_compare_devel_SOURCES  = src/apps/compare.C
devel_compare_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
devel_compare_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
devel_compare_devel_LDADD    = devel/libmesh_devel.la

dbg_programs          += dbg/compare-dbg
dbg_compare_dbg_SOURCES    = src/apps/compare.C
dbg_compare_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
dbg_compare_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
dbg_compare_dbg_LDADD      = dbg/libmesh_dbg.la

# meshbcid
opt_programs           += opt/meshbcid
opt_meshbcid_SOURCES    = src/apps/meshbcid.C
opt_meshbcid_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
opt_meshbcid_CXXFLAGS   = $(CXXFLAGS_OPT)
opt_meshbcid_LDADD      = opt/libmesh.la

devel_programs         += devel/meshbcid-devel
devel_meshbcid_devel_SOURCES  = src/apps/meshbcid.C
devel_meshbcid_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
devel_meshbcid_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
devel_meshbcid_devel_LDADD    = devel/libmesh_devel.la

dbg_programs           += dbg/meshbcid-dbg
dbg_meshbcid_dbg_SOURCES    = src/apps/meshbcid.C
dbg_meshbcid_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
dbg_meshbcid_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
dbg_meshbcid_dbg_LDADD      = dbg/libmesh_dbg.la

# meshid
opt_programs         += opt/meshid
opt_meshid_SOURCES    = src/apps/meshid.C
opt_meshid_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
opt_meshid_CXXFLAGS   = $(CXXFLAGS_OPT)
opt_meshid_LDADD      = opt/libmesh.la

devel_programs        += devel/meshid-devel
devel_meshid_devel_SOURCES  = src/apps/meshid.C
devel_meshid_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
devel_meshid_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
devel_meshid_devel_LDADD    = devel/libmesh_devel.la

dbg_programs         += dbg/meshid-dbg
dbg_meshid_dbg_SOURCES    = src/apps/meshid.C
dbg_meshid_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
dbg_meshid_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
dbg_meshid_dbg_LDADD      = dbg/libmesh_dbg.la

# meshnorm
opt_programs           += opt/meshnorm
opt_meshnorm_SOURCES    = src/apps/meshnorm.C
opt_meshnorm_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
opt_meshnorm_CXXFLAGS   = $(CXXFLAGS_OPT)
opt_meshnorm_LDADD      = opt/libmesh.la

devel_programs         += devel/meshnorm-devel
devel_meshnorm_devel_SOURCES  = src/apps/meshnorm.C
devel_meshnorm_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
devel_meshnorm_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
devel_meshnorm_devel_LDADD    = devel/libmesh_devel.la

dbg_programs           += dbg/meshnorm-dbg
dbg_meshnorm_dbg_SOURCES    = src/apps/meshnorm.C
dbg_meshnorm_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
dbg_meshnorm_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
dbg_meshnorm_dbg_LDADD      = dbg/libmesh_dbg.la

# projection
opt_programs             += opt/projection
opt_projection_SOURCES    = src/apps/projection.C
opt_projection_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
opt_projection_CXXFLAGS   = $(CXXFLAGS_OPT)
opt_projection_LDADD      = opt/libmesh.la

devel_programs           += devel/projection-devel
devel_projection_devel_SOURCES  = src/apps/projection.C
devel_projection_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
devel_projection_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
devel_projection_devel_LDADD    = devel/libmesh_devel.la

dbg_programs             += dbg/projection-dbg
dbg_projection_dbg_SOURCES    = src/apps/projection.C
dbg_projection_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
dbg_projection_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
dbg_projection_dbg_LDADD      = dbg/libmesh_dbg.la

# output_libmesh_version
opt_programs                         += opt/output_libmesh_version
opt_output_libmesh_version_SOURCES    = src/apps/output_libmesh_version.C
opt_output_libmesh_version_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
opt_output_libmesh_version_CXXFLAGS   = $(CXXFLAGS_OPT)
opt_output_libmesh_version_LDADD      = opt/libmesh.la

devel_programs                       += devel/output_libmesh_version-devel
devel_output_libmesh_version_devel_SOURCES  = src/apps/output_libmesh_version.C
devel_output_libmesh_version_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
devel_output_libmesh_version_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
devel_output_libmesh_version_devel_LDADD    = devel/libmesh_devel.la

dbg_programs                         += dbg/output_libmesh_version-dbg
dbg_output_libmesh_version_dbg_SOURCES    = src/apps/output_libmesh_version.C
dbg_output_libmesh_version_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
dbg_output_libmesh_version_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
dbg_output_libmesh_version_dbg_LDADD      = dbg/libmesh_dbg.la

# meshplot
opt_programs           += opt/meshplot
opt_meshplot_SOURCES    = src/apps/meshplot.C
opt_meshplot_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
opt_meshplot_CXXFLAGS   = $(CXXFLAGS_OPT)
opt_meshplot_LDADD      = opt/libmesh.la

devel_programs         += devel/meshplot-devel
devel_meshplot_devel_SOURCES  = src/apps/meshplot.C
devel_meshplot_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
devel_meshplot_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
devel_meshplot_devel_LDADD    = devel/libmesh_devel.la

dbg_programs           += dbg/meshplot-dbg
dbg_meshplot_dbg_SOURCES    = src/apps/meshplot.C
dbg_meshplot_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
dbg_meshplot_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
dbg_meshplot_dbg_LDADD      = dbg/libmesh_dbg.la

# solution_components
opt_programs                      += opt/solution_components
opt_solution_components_SOURCES    = src/apps/solution_components.C
opt_solution_components_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
opt_solution_components_CXXFLAGS   = $(CXXFLAGS_OPT)
opt_solution_components_LDADD      = opt/libmesh.la

devel_programs                    += devel/solution_components-devel
devel_solution_components_devel_SOURCES  = src/apps/solution_components.C
devel_solution_components_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
devel_solution_components_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
devel_solution_components_devel_LDADD    = devel/libmesh_devel.la

dbg_programs                      += dbg/solution_components-dbg
dbg_solution_components_dbg_SOURCES    = src/apps/solution_components.C
dbg_solution_components_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
dbg_solution_components_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
dbg_solution_components_dbg_LDADD      = dbg/libmesh_dbg.la


if LIBMESH_OPT_MODE
  bin_PROGRAMS += $(opt_programs)
endif

if LIBMESH_DEVEL_MODE
  bin_PROGRAMS += $(devel_programs)
endif

if LIBMESH_DBG_MODE
  bin_PROGRAMS += $(dbg_programs)
endif


# Support for libmesh-config in $(exec_prefix)/bin
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
