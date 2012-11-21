check_PROGRAMS       = adjoints_ex4
check_SCRIPTS        = run.sh
adjoints_ex4dir      = $(examples_install_path)/adjoints/ex4
adjoints_ex4_DATA    = L-shaped.C element_qoi_derivative.C general.in side_postprocess.C L-shaped.h adjoints_ex4.C femparameters.C l-shaped.in side_qoi_derivative.C element_postprocess.C femparameters.h lshaped.xda run.sh $(top_srcdir)/contrib/utils/Makefile
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

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
