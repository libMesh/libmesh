example_name = eigenproblems_ex3
install_dir  = $(examples_install_path)/eigenproblems/ex3
data         = eigenproblems_ex3.C drum1_mesh.e drum2_mesh.e run.sh
sources      = $(data) run.sh

CLEANFILES                = drum1_evals.txt drum2_evals.txt drum1_evec.e drum2_evec.e

# conditionally enable this test
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f drum1_mesh.e && $(LN_S) -f $(srcdir)/drum1_mesh.e .
	-rm -f drum2_mesh.e && $(LN_S) -f $(srcdir)/drum2_mesh.e .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += drum*_mesh.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
