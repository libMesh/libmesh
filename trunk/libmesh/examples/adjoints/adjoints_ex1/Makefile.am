check_PROGRAMS       = adjoints_ex1
check_SCRIPTS        = run.sh
adjoints_ex1dir      = $(examples_install_path)/adjoints/ex1
adjoints_ex1_DATA    = adjoints_ex1.C L-shaped.C L-shaped.h element_postprocess.C element_qoi_derivative.C femparameters.C femparameters.h general.in l-shaped.in lshaped.xda side_postprocess.C side_qoi_derivative.C
adjoints_ex1_SOURCES = $(adjoints_ex1_DATA) run.sh
adjoints_ex1_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = primal_error.gmv primal.out.* dual_error.gmv adjoint_?.out.*

# also need to link mesh.xda for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f general.in  && $(LN_S) $(srcdir)/general.in .
	-rm -f l-shaped.in && $(LN_S) $(srcdir)/l-shaped.in .
	-rm -f lshaped.xda && $(LN_S) $(srcdir)/lshaped.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in l-shaped.in lshaped.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
