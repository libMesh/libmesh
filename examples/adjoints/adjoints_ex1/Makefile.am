example_name  = adjoints_ex1
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adjoints/ex1
data          = adjoints_ex1.C L-shaped.C L-shaped.h element_postprocess.C element_qoi_derivative.C femparameters.C femparameters.h general.in l-shaped.in lshaped.xda side_postprocess.C side_qoi_derivative.C run.sh
sources       = $(data) run.sh

CLEANFILES = primal_error.gmv primal.out.* dual_error.gmv adjoint_?.out.* primal.e* adjoint_0.e* adjoint_1.e*

# also need to link mesh.xda for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f general.in  && $(LN_S) -f $(srcdir)/general.in .
	-rm -f l-shaped.in && $(LN_S) -f $(srcdir)/l-shaped.in .
	-rm -f lshaped.xda && $(LN_S) -f $(srcdir)/lshaped.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in l-shaped.in lshaped.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
