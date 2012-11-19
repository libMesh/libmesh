include $(top_builddir)/Make.common.@METHOD@

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS  = $(libmesh_CPPFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CPPFLAGS += $(libmesh_INCLUDE) -I$(top_builddir)/include
LIBS         = $(libmesh_LIBS)

if LIBMESH_ENABLE_FPARSER
  AM_CPPFLAGS += -I$(top_srcdir)/contrib/fparser
endif

BUILT_SOURCES = mesh.xda

mesh.xda: $(top_srcdir)/examples/adaptivity/adaptivity_ex2/mesh.xda
	$(AM_V_GEN)$(LN_S) $< $@

DISTCLEANFILES = $(BUILT_SOURCES)

check_PROGRAMS        = transient_ex1
check_SCRIPTS         = run.sh
transient_ex1dir      = $(examples_install_path)/transient/ex1
transient_ex1_DATA    = transient_ex1.C exact_solution.C ../../adaptivity/adaptivity_ex2/mesh.xda
transient_ex1_SOURCES = $(transient_ex1_DATA) run.sh
transient_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.gmv

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


# how to make documentation
CLEANFILES += $(check_PROGRAMS).php stdout.log
doc: $(check_PROGRAMS).php
html: $(check_PROGRAMS).php
$(check_PROGRAMS).php: $(check_PROGRAMS) \
                        $(top_srcdir)/contrib/bin/ex2html.sh \
                       Makefile
	$(MAKE) -s run > stdout.log
	$(AM_V_GEN) PATH=$(abs_top_srcdir)/contrib/bin:$(PATH) \
	  ex2html.sh $(check_PROGRAMS) $(abs_srcdir)

