#bindir=@bindir@/@host@_@METHOD@
#libdir=@libdir@/@host@_@METHOD@

SUBDIRS = include #src

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS  = -I m4 -I m4/autoconf-submodule
AM_CFLAGS        = $(libmesh_CFLAGS)
AM_LDFLAGS       = $(libmesh_LDFLAGS)

# We might have turned on -Werror and/or paranoid warnings
CXXFLAGS_DBG += $(ANY_WERROR_FLAG) $(ANY_PARANOID_FLAGS)
CXXFLAGS_DEVEL += $(ANY_WERROR_FLAG) $(ANY_PARANOID_FLAGS)
CXXFLAGS_OPROF += $(ANY_WERROR_FLAG) $(ANY_PARANOID_FLAGS)
CXXFLAGS_OPT += $(ANY_WERROR_FLAG) $(ANY_PARANOID_FLAGS)
CXXFLAGS_PROF += $(ANY_WERROR_FLAG) $(ANY_PARANOID_FLAGS)

# We need to set the no-undefined flag for linking
# with libtool, because Windows requires all shared
# libraries to include all used symbols. This flag
# should do no harm on other platforms.
AM_LDFLAGS += -no-undefined

# As discussed on this thread [0], automake passes -h to tar in the
# "make dist" target, and this has the effect of creating hard links
# out of soft links in the resulting archives.  Since not all
# filesystems support hard links well (for example AFS, see libmesh
# #2046) this option causes the directory links to be created as
# normal files instead. The --hard-dereference option requires at least tar 1.29 to
# work, but we generally completely control the systems where we run
# "make dist" and therefore we don't need the most portable "make
# dist" possible, but rather a "make dist" target that gives us the
# most portable tarballs possible.
# [0]: http://gnu-automake.7480.n7.nabble.com/bug-19616-dist-tarball-contains-hardlinks-td21667.html
TAR_OPTIONS = --hard-dereference
export TAR_OPTIONS

# AM_CPPFLAGS are method-independent cppflags that
# we use when compiling libmesh proper, or its utility
# codes.  These will be augmented with method-specific
# cppflags later on
AM_CPPFLAGS      = -DLIBMESH_IS_COMPILING_ITSELF \
                   $(libmesh_contrib_INCLUDES) \
                   $(libmesh_optional_INCLUDES) \
                   -I$(top_builddir)/include # required for libmesh_version.h

LIBS             = $(libmesh_optional_LIBS) $(libmesh_precision_LIBS)

# additional files which must be included in 'make dist'
EXTRA_DIST = reference_elements \
             include/include_HEADERS \
             include/rebuild_include_HEADERS.sh \
             src/libmesh_SOURCES \
             src/rebuild_libmesh_SOURCES.sh \
             contrib/tecplot/README \
             contrib/bin/strip_dup_incl_paths.pl \
             contrib/bin/strip_dup_libs.pl \
             CITATION \
             contrib/metaphysicl/README \
             m4/autoconf-submodule/ax_split_version.m4 \
             contrib/timpi/m4/autoconf-submodule/ax_split_version.m4 \
             contrib/metaphysicl/m4/autoconf-submodule/ax_split_version.m4

# If our submodules have not all been initialized, throw an
# error during "make dist". We need the submodules (and their
# submodules) to be properly initalized so that the sources can be
# copied into the distribution, and trying to do it "on the fly" is
# too complicated.
contrib/metaphysicl/README:
	@echo ""
	@echo "Error: You must run git submodule update --init --recursive first"
	@echo ""
	@exit 1

m4/autoconf-submodule/ax_split_version.m4:
	@echo ""
	@echo "Error: You must run git submodule update --init --recursive first"
	@echo ""
	@exit 1

contrib/timpi/m4/autoconf-submodule/ax_split_version.m4:
	@echo ""
	@echo "Error: You must run git submodule update --init --recursive first"
	@echo ""
	@exit 1

contrib/metaphysicl/m4/autoconf-submodule/ax_split_version.m4:
	@echo ""
	@echo "Error: You must run git submodule update --init --recursive first"
	@echo ""
	@exit 1

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
        reference_elements/3D/one_pyramid13.xda \
        reference_elements/3D/one_pyramid14.xda \
        reference_elements/3D/one_tet.xda \
        reference_elements/3D/one_tet10.xda

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh \
           build-aux/missing

# Support for pkgconfig
pkgconfigdir     = $(libdir)/pkgconfig
pkgconfig_DATA   = contrib/utils/Make.common \
                   contrib/utils/libmesh.pc

configdir        = $(sysconfdir)/libmesh
config_DATA      = contrib/utils/Make.common \
                   contrib/utils/libmesh.pc

# Install pkgconfig files for enabled build modes
if LIBMESH_OPT_MODE
  pkgconfig_DATA += contrib/utils/libmesh-opt.pc
  config_DATA += contrib/utils/libmesh-opt.pc
endif
if LIBMESH_DBG_MODE
  pkgconfig_DATA += contrib/utils/libmesh-dbg.pc
  config_DATA += contrib/utils/libmesh-dbg.pc
endif
if LIBMESH_DEVEL_MODE
  pkgconfig_DATA += contrib/utils/libmesh-devel.pc
  config_DATA += contrib/utils/libmesh-devel.pc
endif
if LIBMESH_PROF_MODE
  pkgconfig_DATA += contrib/utils/libmesh-prof.pc
  config_DATA += contrib/utils/libmesh-prof.pc
endif
if LIBMESH_OPROF_MODE
  pkgconfig_DATA += contrib/utils/libmesh-oprof.pc
  config_DATA += contrib/utils/libmesh-oprof.pc
endif

# install the libtool script so that dependent apps can also use it
contribbindir      = $(prefix)/contrib/bin
contribbin_SCRIPTS = libtool

# Support for libmesh-config in $(exec_prefix)/bin
libmesh_configdir      = $(exec_prefix)/bin
libmesh_config_SCRIPTS = contrib/bin/libmesh-config

# This rule extends the built-in make install rule, and can
# be used to:
#  when doing 'make install prefix=/other/path', update generated files
#  when doing 'make install', symlik etc/libmesh/Make.common
#  when doing 'make install', update libmesh_config
#  fix up buggy libtool postdeps under specific compiler/mpi configurations
install-data-hook:
	@if (test "x$(prefix)" != "x@prefix@") ; then \
          oldprefix="@prefix@" ; \
          newprefix="$(prefix)" ; \
          for genfile in `echo $(DESTDIR)$(prefix)/lib/pkgconfig/*.pc $(DESTDIR)$(prefix)/etc/libmesh/*.pc $(DESTDIR)$(prefix)/etc/libmesh/Make.common $(DESTDIR)$(prefix)/lib/pkgconfig/Make.common $(DESTDIR)$(prefix)/bin/libmesh-config` ; do \
            echo " " ; \
            echo " *** replacing $$oldprefix" ; \
            echo " ***      with $$newprefix" ; \
            echo " *** in generated file $$genfile" ; \
            echo " " ; \
            cd $(DESTDIR)$(prefix) && cat $$genfile | $(SED) "s,$$oldprefix,$$newprefix,g" > $$genfile.replaced ; \
            cd $(DESTDIR)$(prefix) && mv $$genfile.replaced $$genfile ; \
          done ; \
        fi
	cd $(DESTDIR)$(prefix) && rm -f Make.common && $(LN_S) -f etc/libmesh/Make.common Make.common
	cat $(DESTDIR)$(libmesh_configdir)/libmesh-config | $(SED) "s/has_been_installed=no/has_been_installed=yes/g" > $(DESTDIR)$(libmesh_configdir)/libmesh-config.installed
	mv $(DESTDIR)$(libmesh_configdir)/libmesh-config.installed $(DESTDIR)$(libmesh_configdir)/libmesh-config && chmod +x $(DESTDIR)$(libmesh_configdir)/libmesh-config
	cat $(DESTDIR)$(prefix)/contrib/bin/libtool | $(SED) "s/-l //g" > $(DESTDIR)$(prefix)/contrib/bin/libtool.fixed
	mv $(DESTDIR)$(prefix)/contrib/bin/libtool.fixed $(DESTDIR)$(prefix)/contrib/bin/libtool && chmod +x $(DESTDIR)$(prefix)/contrib/bin/libtool


# this rule invokes emacs on each source file to remove trailing whitespace.
beautify:
	for file in `find $(srcdir)/{src,include,examples,tests} -name "*.C" -o -name "*.h" -o -name Makefile.am -type f` ; do \
          echo $$file ; \
          emacs -batch $$file --eval '(delete-trailing-whitespace)' -f save-buffer 2>/dev/null ; \
        done


# Additional files to be deleted by 'make distclean'
# _configs.sed seems to be leftover by prefix_config
DISTCLEANFILES =  \
           _configs.sed \
	   a.out  # intel-14.0 seems to leave an a.out from configure with mpich, clean that up for distclean

# # Additional files to be deleted by 'make maintainer-clean'
# MAINTAINERCLEANFILES =  \
#          $(AUX_DIST)                              \
#          Makefile.in                              \
#          aclocal.m4                               \
#          aminclude.am                             \
#          autom4te.cache/*                         \
#          config.h                                 \
#          config.h.in                              \
#          config.log                               \
#          config.status                            \
#          config.sub                               \
#          configure                                \
#          data/Makefile.in                         \
#          src/Makefile.in                          \
#          src/physics/properties/Makefile.in       \
#          test/Makefile.in                         \
#          test/common/Makefile.in                  \
#          test/comp_ns/Makefile.in                 \
#          test/unit/Makefile.in



SUBDIRS += contrib

include src/libmesh_SOURCES

# A convenience library to hold proper libMesh
# objects.  This will get appended with any contributed
# sources to create the final library.
lib_LTLIBRARIES =

if LIBMESH_DBG_MODE
  lib_LTLIBRARIES        += libmesh_dbg.la
  libmesh_dbg_la_SOURCES  = $(libmesh_SOURCES)
  libmesh_dbg_la_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
  libmesh_dbg_la_CXXFLAGS = $(CXXFLAGS_DBG)
  libmesh_dbg_la_CFLAGS   = $(CFLAGS_DBG)
  libmesh_dbg_la_LIBADD   = contrib/libcontrib_dbg.la $(LIBS)
endif

if LIBMESH_DEVEL_MODE
  lib_LTLIBRARIES          += libmesh_devel.la
  libmesh_devel_la_SOURCES  = $(libmesh_SOURCES)
  libmesh_devel_la_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
  libmesh_devel_la_CXXFLAGS = $(CXXFLAGS_DEVEL)
  libmesh_devel_la_CFLAGS   = $(CFLAGS_DEVEL)
  libmesh_devel_la_LIBADD   = contrib/libcontrib_devel.la $(LIBS)
endif

if LIBMESH_OPT_MODE
  lib_LTLIBRARIES        += libmesh_opt.la
  libmesh_opt_la_SOURCES  = $(libmesh_SOURCES)
  libmesh_opt_la_CPPFLAGS = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
  libmesh_opt_la_CXXFLAGS = $(CXXFLAGS_OPT)
  libmesh_opt_la_CFLAGS   = $(CFLAGS_OPT)
  libmesh_opt_la_LIBADD   = contrib/libcontrib_opt.la $(LIBS)
endif

if LIBMESH_PROF_MODE
  lib_LTLIBRARIES          += libmesh_prof.la
  libmesh_prof_la_SOURCES  = $(libmesh_SOURCES)
  libmesh_prof_la_CPPFLAGS = $(CPPFLAGS_PROF) $(AM_CPPFLAGS)
  libmesh_prof_la_CXXFLAGS = $(CXXFLAGS_PROF)
  libmesh_prof_la_CFLAGS   = $(CFLAGS_PROF)
  libmesh_prof_la_LIBADD   = contrib/libcontrib_prof.la $(LIBS)
endif

if LIBMESH_OPROF_MODE
  lib_LTLIBRARIES          += libmesh_oprof.la
  libmesh_oprof_la_SOURCES  = $(libmesh_SOURCES)
  libmesh_oprof_la_CPPFLAGS = $(CPPFLAGS_OPROF) $(AM_CPPFLAGS)
  libmesh_oprof_la_CXXFLAGS = $(CXXFLAGS_OPROF)
  libmesh_oprof_la_CFLAGS   = $(CFLAGS_OPROF)
  libmesh_oprof_la_LIBADD   = contrib/libcontrib_oprof.la $(LIBS)
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
prof_programs  = # empty, append below
oprof_programs = # empty, append below

# fparser_parse
opt_programs               += fparser_parse-opt
fparser_parse_opt_SOURCES    = src/apps/fparser_parse.C
fparser_parse_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
fparser_parse_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
fparser_parse_opt_LDADD      = libmesh_opt.la

devel_programs             += fparser_parse-devel
fparser_parse_devel_SOURCES  = src/apps/fparser_parse.C
fparser_parse_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
fparser_parse_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
fparser_parse_devel_LDADD    = libmesh_devel.la

dbg_programs               += fparser_parse-dbg
fparser_parse_dbg_SOURCES    = src/apps/fparser_parse.C
fparser_parse_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
fparser_parse_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
fparser_parse_dbg_LDADD      = libmesh_dbg.la

# getpot_parse
opt_programs               += getpot_parse-opt
getpot_parse_opt_SOURCES    = src/apps/getpot_parse.C
getpot_parse_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
getpot_parse_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
getpot_parse_opt_LDADD      = libmesh_opt.la

devel_programs             += getpot_parse-devel
getpot_parse_devel_SOURCES  = src/apps/getpot_parse.C
getpot_parse_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
getpot_parse_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
getpot_parse_devel_LDADD    = libmesh_devel.la

dbg_programs               += getpot_parse-dbg
getpot_parse_dbg_SOURCES    = src/apps/getpot_parse.C
getpot_parse_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
getpot_parse_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
getpot_parse_dbg_LDADD      = libmesh_dbg.la

# amr
opt_programs    += amr-opt
amr_opt_SOURCES  = src/apps/amr.C
amr_opt_CPPFLAGS = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
amr_opt_CXXFLAGS = $(CXXFLAGS_OPT)
amr_opt_LDADD    = libmesh_opt.la

devel_programs    += amr-devel
amr_devel_SOURCES  = src/apps/amr.C
amr_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
amr_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
amr_devel_LDADD    = libmesh_devel.la

dbg_programs    += amr-dbg
amr_dbg_SOURCES  = src/apps/amr.C
amr_dbg_CPPFLAGS = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
amr_dbg_CXXFLAGS = $(CXXFLAGS_DBG)
amr_dbg_LDADD    = libmesh_dbg.la

# meshtool
opt_programs           += meshtool-opt
meshtool_opt_SOURCES    = src/apps/meshtool.C
meshtool_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
meshtool_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
meshtool_opt_LDADD      = libmesh_opt.la

devel_programs         += meshtool-devel
meshtool_devel_SOURCES  = src/apps/meshtool.C
meshtool_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
meshtool_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
meshtool_devel_LDADD    = libmesh_devel.la

dbg_programs           += meshtool-dbg
meshtool_dbg_SOURCES    = src/apps/meshtool.C
meshtool_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
meshtool_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
meshtool_dbg_LDADD      = libmesh_dbg.la

# calculator
opt_programs          += calculator-opt
calculator_opt_SOURCES    = src/apps/calculator.C
calculator_opt_SOURCES   += src/apps/L2system.C src/apps/L2system.h
calculator_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
calculator_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
calculator_opt_LDADD      = libmesh_opt.la

devel_programs        += calculator-devel
calculator_devel_SOURCES  = src/apps/calculator.C
calculator_devel_SOURCES += src/apps/L2system.C src/apps/L2system.h
calculator_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
calculator_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
calculator_devel_LDADD    = libmesh_devel.la

dbg_programs          += calculator-dbg
calculator_dbg_SOURCES    = src/apps/calculator.C
calculator_dbg_SOURCES   += src/apps/L2system.C src/apps/L2system.h
calculator_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
calculator_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
calculator_dbg_LDADD      = libmesh_dbg.la

# compare
opt_programs          += compare-opt
compare_opt_SOURCES    = src/apps/compare.C
compare_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
compare_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
compare_opt_LDADD      = libmesh_opt.la

devel_programs        += compare-devel
compare_devel_SOURCES  = src/apps/compare.C
compare_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
compare_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
compare_devel_LDADD    = libmesh_devel.la

dbg_programs          += compare-dbg
compare_dbg_SOURCES    = src/apps/compare.C
compare_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
compare_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
compare_dbg_LDADD      = libmesh_dbg.la

# meshbcid
opt_programs           += meshbcid-opt
meshbcid_opt_SOURCES    = src/apps/meshbcid.C
meshbcid_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
meshbcid_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
meshbcid_opt_LDADD      = libmesh_opt.la

devel_programs         += meshbcid-devel
meshbcid_devel_SOURCES  = src/apps/meshbcid.C
meshbcid_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
meshbcid_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
meshbcid_devel_LDADD    = libmesh_devel.la

dbg_programs           += meshbcid-dbg
meshbcid_dbg_SOURCES    = src/apps/meshbcid.C
meshbcid_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
meshbcid_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
meshbcid_dbg_LDADD      = libmesh_dbg.la

# meshid
opt_programs         += meshid-opt
meshid_opt_SOURCES    = src/apps/meshid.C
meshid_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
meshid_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
meshid_opt_LDADD      = libmesh_opt.la

devel_programs        += meshid-devel
meshid_devel_SOURCES  = src/apps/meshid.C
meshid_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
meshid_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
meshid_devel_LDADD    = libmesh_devel.la

dbg_programs         += meshid-dbg
meshid_dbg_SOURCES    = src/apps/meshid.C
meshid_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
meshid_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
meshid_dbg_LDADD      = libmesh_dbg.la

# meshavg
opt_programs           += meshavg-opt
meshavg_opt_SOURCES    = src/apps/meshavg.C
meshavg_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
meshavg_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
meshavg_opt_LDADD      = libmesh_opt.la

devel_programs         += meshavg-devel
meshavg_devel_SOURCES  = src/apps/meshavg.C
meshavg_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
meshavg_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
meshavg_devel_LDADD    = libmesh_devel.la

dbg_programs           += meshavg-dbg
meshavg_dbg_SOURCES    = src/apps/meshavg.C
meshavg_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
meshavg_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
meshavg_dbg_LDADD      = libmesh_dbg.la

# meshdiff
opt_programs           += meshdiff-opt
meshdiff_opt_SOURCES    = src/apps/meshdiff.C
meshdiff_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
meshdiff_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
meshdiff_opt_LDADD      = libmesh_opt.la

devel_programs         += meshdiff-devel
meshdiff_devel_SOURCES  = src/apps/meshdiff.C
meshdiff_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
meshdiff_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
meshdiff_devel_LDADD    = libmesh_devel.la

dbg_programs           += meshdiff-dbg
meshdiff_dbg_SOURCES    = src/apps/meshdiff.C
meshdiff_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
meshdiff_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
meshdiff_dbg_LDADD      = libmesh_dbg.la

# meshnorm
opt_programs           += meshnorm-opt
meshnorm_opt_SOURCES    = src/apps/meshnorm.C
meshnorm_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
meshnorm_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
meshnorm_opt_LDADD      = libmesh_opt.la

devel_programs         += meshnorm-devel
meshnorm_devel_SOURCES  = src/apps/meshnorm.C
meshnorm_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
meshnorm_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
meshnorm_devel_LDADD    = libmesh_devel.la

dbg_programs           += meshnorm-dbg
meshnorm_dbg_SOURCES    = src/apps/meshnorm.C
meshnorm_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
meshnorm_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
meshnorm_dbg_LDADD      = libmesh_dbg.la

# projection
opt_programs             += projection-opt
projection_opt_SOURCES    = src/apps/projection.C
projection_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
projection_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
projection_opt_LDADD      = libmesh_opt.la

devel_programs           += projection-devel
projection_devel_SOURCES  = src/apps/projection.C
projection_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
projection_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
projection_devel_LDADD    = libmesh_devel.la

dbg_programs             += projection-dbg
projection_dbg_SOURCES    = src/apps/projection.C
projection_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
projection_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
projection_dbg_LDADD      = libmesh_dbg.la

# output_libmesh_version
opt_programs                         += output_libmesh_version-opt
output_libmesh_version_opt_SOURCES    = src/apps/output_libmesh_version.C
output_libmesh_version_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
output_libmesh_version_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
output_libmesh_version_opt_LDADD      = libmesh_opt.la

devel_programs                       += output_libmesh_version-devel
output_libmesh_version_devel_SOURCES  = src/apps/output_libmesh_version.C
output_libmesh_version_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
output_libmesh_version_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
output_libmesh_version_devel_LDADD    = libmesh_devel.la

dbg_programs                         += output_libmesh_version-dbg
output_libmesh_version_dbg_SOURCES    = src/apps/output_libmesh_version.C
output_libmesh_version_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
output_libmesh_version_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
output_libmesh_version_dbg_LDADD      = libmesh_dbg.la

# meshplot
opt_programs           += meshplot-opt
meshplot_opt_SOURCES    = src/apps/meshplot.C
meshplot_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
meshplot_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
meshplot_opt_LDADD      = libmesh_opt.la

devel_programs         += meshplot-devel
meshplot_devel_SOURCES  = src/apps/meshplot.C
meshplot_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
meshplot_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
meshplot_devel_LDADD    = libmesh_devel.la

dbg_programs           += meshplot-dbg
meshplot_dbg_SOURCES    = src/apps/meshplot.C
meshplot_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
meshplot_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
meshplot_dbg_LDADD      = libmesh_dbg.la

# solution_components
opt_programs                      += solution_components-opt
solution_components_opt_SOURCES    = src/apps/solution_components.C
solution_components_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
solution_components_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
solution_components_opt_LDADD      = libmesh_opt.la

devel_programs                    += solution_components-devel
solution_components_devel_SOURCES  = src/apps/solution_components.C
solution_components_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
solution_components_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
solution_components_devel_LDADD    = libmesh_devel.la

dbg_programs                      += solution_components-dbg
solution_components_dbg_SOURCES    = src/apps/solution_components.C
solution_components_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
solution_components_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
solution_components_dbg_LDADD      = libmesh_dbg.la

# splitter
opt_programs           += splitter-opt
splitter_opt_SOURCES    = src/apps/splitter.C
splitter_opt_CPPFLAGS   = $(CPPFLAGS_OPT) $(AM_CPPFLAGS)
splitter_opt_CXXFLAGS   = $(CXXFLAGS_OPT)
splitter_opt_LDADD      = libmesh_opt.la

devel_programs         += splitter-devel
splitter_devel_SOURCES  = src/apps/splitter.C
splitter_devel_CPPFLAGS = $(CPPFLAGS_DEVEL) $(AM_CPPFLAGS)
splitter_devel_CXXFLAGS = $(CXXFLAGS_DEVEL)
splitter_devel_LDADD    = libmesh_devel.la

dbg_programs           += splitter-dbg
splitter_dbg_SOURCES    = src/apps/splitter.C
splitter_dbg_CPPFLAGS   = $(CPPFLAGS_DBG) $(AM_CPPFLAGS)
splitter_dbg_CXXFLAGS   = $(CXXFLAGS_DBG)
splitter_dbg_LDADD      = libmesh_dbg.la

if LIBMESH_OPT_MODE
  bin_PROGRAMS += $(opt_programs)
endif

if LIBMESH_DEVEL_MODE
  bin_PROGRAMS += $(devel_programs)
endif

if LIBMESH_DBG_MODE
  bin_PROGRAMS += $(dbg_programs)
endif

###########################################################
# Examples
if LIBMESH_ENABLE_EXAMPLES
  SUBDIRS += examples
endif

###########################################################
# Documentation
.PHONY: examples_doc doc

SUBDIRS += doc

#
# documentation for the examples
examples_doc:
	@cd $(top_builddir)/examples && $(MAKE) doc

#
# doxygen documentation - be sure to generate symlinks in include/libmesh!
doc html:
#	$(MAKE) examples_doc
	$(MAKE) -C doc/html doc

# Upload the web page to github.io.
upload:
	-rm -rf libMesh.github.io && git clone https://github.com/libMesh/libMesh.github.io.git
	cd libMesh.github.io && git rm -rf doxygen
	rsync -aLv --exclude=src --exclude=Makefile --exclude=.linkstamp $(top_builddir)/doc/html/ libMesh.github.io/
	gitcommit=`git log | head -1` && cd libMesh.github.io && git add * && git ci -m"updating web page from \"make upload\", libmesh.git at $$gitcommit"
# uncomment below when ready
#	cd libMesh.github.io && git push
#	-rm -rf libMesh.github.io


#
# special things to do when running 'make dist'
dist-hook:
	rm -rf `find $(distdir) -name libmesh_config.h`


#
# support 'make run_examples'
run_examples:
	@$(MAKE) && cd $(top_builddir)/examples && $(MAKE) check

#
# support top-level 'make test_headers'
test_headers:
	@$(MAKE) -C $(top_builddir)/include/libmesh
	@cd $(top_builddir)/include && $(MAKE) test_headers

# -------------------------------------------
# Optional support for code coverage analysis
# -------------------------------------------

if CODE_COVERAGE_ENABLED

lcov_dir=$(top_builddir)/doc/lcov/html
lcov_bin=$(top_srcdir)/contrib/lcov/lcov-1.10/bin

# General philosophy is to maintain code coverage for the
# base library as generated by "make check" tests.

lcov-report:
	@mkdir -p $(lcov_dir)
	$(lcov_bin)/lcov --compat-libtool --directory . --capture --output-file $(lcov_dir)/lcov.info
	$(lcov_bin)/lcov --list-full-path -l $(lcov_dir)/lcov.info | grep -v "`cd -P $(top_srcdir)/src && pwd`" | cut -d\| -f1 > $(lcov_dir)/remove
	$(lcov_bin)/lcov -q -r $(lcov_dir)/lcov.info `cat $(lcov_dir)/remove` > $(lcov_dir)/lcov.cleaned.info
	@rm $(lcov_dir)/remove
	@mv $(lcov_dir)/lcov.cleaned.info $(lcov_dir)/lcov.info
	$(lcov_bin)/genhtml -t "libMesh" -o $(lcov_dir) $(lcov_dir)/lcov.info

lcov-reset:
	@rm -rf $(lcov_dir)
	@find . -name "*.gcda" -exec rm {} \;
	$(lcov_bin)/lcov --directory . --zerocounters

coverage: lcov-reset check lcov-report

CLEANFILES += src/apps/*.gcda src/apps/*.gcno

endif
