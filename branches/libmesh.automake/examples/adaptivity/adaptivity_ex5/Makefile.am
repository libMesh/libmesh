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
exact_solution.C: $(top_srcdir)/examples/transient/transient_ex1/exact_solution.C 
	$(AM_V_GEN)$(LN_S) $< $@

check_PROGRAMS         = adaptivity_ex5
check_SCRIPTS          = run.sh
adaptivity_ex5dir      = $(examples_install_path)/adaptivity/ex5
adaptivity_ex5_DATA    = adaptivity_ex5.C exact_solution.C
adaptivity_ex5_SOURCES = $(adaptivity_ex5_DATA) run.sh
adaptivity_ex5_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = solution_read_in.gmv solution_read_in.e saved_solution.xdr saved_solution.gmv saved_solution.e saved_mesh.xdr out_*.e out.gmv.*

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


