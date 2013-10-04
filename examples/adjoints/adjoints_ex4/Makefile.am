example_name  = adjoints_ex4
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adjoints/ex4
data          = L-shaped.C element_qoi_derivative.C general.in side_postprocess.C L-shaped.h adjoints_ex4.C femparameters.C l-shaped.in side_qoi_derivative.C element_postprocess.C femparameters.h lshaped.xda run.sh
sources       = $(data) run.sh

CLEANFILES = primal.out.* adjoint_*.out.*

# also need links for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f general.in  && $(LN_S) $(srcdir)/general.in .
	-rm -f lshaped.xda && $(LN_S) $(srcdir)/lshaped.xda .
	-rm -f l-shaped.in && $(LN_S) $(srcdir)/l-shaped.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in lshaped.xda l-shaped.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
