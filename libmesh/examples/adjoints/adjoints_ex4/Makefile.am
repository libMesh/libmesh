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

check_PROGRAMS       = adjoints_ex4
check_SCRIPTS        = run.sh
adjoints_ex4dir      = $(examples_install_path)/adjoints/ex4
adjoints_ex4_DATA    = L-shaped.C element_qoi_derivative.C general.in side_postprocess.C L-shaped.h adjoints_ex4.C femparameters.C l-shaped.in side_qoi_derivative.C element_postprocess.C femparameters.h lshaped.xda
adjoints_ex4_SOURCES = $(adjoints_ex4_DATA) run.sh
adjoints_ex4_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = primal.out.* adjoint_*.out.*

# also need links for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f general.in  && $(LN_S) $(srcdir)/general.in .
	-rm -f lshaped.xda && $(LN_S) $(srcdir)/lshaped.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in lshaped.xda .linkstamp
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

