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

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
