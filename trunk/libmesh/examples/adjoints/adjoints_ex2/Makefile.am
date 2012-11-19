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

# need to link lshaped.xda
BUILT_SOURCES = lshaped.xda
lshaped.xda:
	$(AM_V_GEN)$(LN_S) $(top_srcdir)/examples/adjoints/adjoints_ex1/lshaped.xda .

CLEANFILES = lshaped.xda

check_PROGRAMS       = adjoints_ex2
check_SCRIPTS        = run.sh
adjoints_ex2dir      = $(examples_install_path)/adjoints/ex2
adjoints_ex2_DATA    = adjoints_ex2.C ../adjoints_ex1/lshaped.xda L-shaped.C L-shaped.h femparameters.C femparameters.h general.in l-shaped.in L-qoi.C L-qoi.h
adjoints_ex2_SOURCES = $(adjoints_ex2_DATA) run.sh
adjoints_ex2_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES += primal.out.* adjoint_0.out.*

# also need to link input for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f general.in  && $(LN_S) $(srcdir)/general.in .
	-rm -f l-shaped.in && $(LN_S) $(srcdir)/l-shaped.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in l-shaped.in .linkstamp
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

