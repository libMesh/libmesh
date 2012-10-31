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

BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/introduction/introduction_ex3/exact_solution.C 
	$(AM_V_GEN)$(LN_S) $< $@

check_PROGRAMS         = subdomains_ex2
check_SCRIPTS          = run.sh
subdomains_ex2dir      = $(examples_install_path)/subdomains/ex2
subdomains_ex2_DATA    = subdomains_ex2.C exact_solution.C
subdomains_ex2_SOURCES = $(subdomains_ex2_DATA) run.sh
subdomains_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = out_*.e gnuplot_*

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


