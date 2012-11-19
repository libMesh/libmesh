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
exact_solution.C: $(top_srcdir)/examples/transient/transient_ex1/exact_solution.C 
	$(AM_V_GEN)$(LN_S) $< $@

DISTCLEANFILES = $(BUILT_SOURCES)

check_PROGRAMS         = adaptivity_ex2
check_SCRIPTS          = run.sh
adaptivity_ex2dir      = $(examples_install_path)/adaptivity/ex2
adaptivity_ex2_DATA    = adaptivity_ex2.C exact_solution.C mesh.xda
adaptivity_ex2_SOURCES = $(adaptivity_ex2_DATA) run.sh
adaptivity_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = solution_read_in.gmv saved_solution.xda saved_solution.gmv saved_mesh.xda out.gmv.*

# also need to link mesh.xda for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f mesh.xda && $(LN_S) $(srcdir)/mesh.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += mesh.xda .linkstamp
endif

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

