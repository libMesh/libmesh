check_PROGRAMS            = eigenproblems_ex3
eigenproblems_ex3dir      = $(examples_install_path)/eigenproblems/ex3
eigenproblems_ex3_DATA    = eigenproblems_ex3.C drum1_mesh.e drum2_mesh.e run.sh $(top_srcdir)/contrib/utils/Makefile
eigenproblems_ex3_SOURCES = $(eigenproblems_ex3_DATA) run.sh
eigenproblems_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES                = drum1_evals.txt drum2_evals.txt drum1_evec.e drum2_evec.e

# conditionally enable this test
if LIBMESH_ENABLE_SLEPC
  check_SCRIPTS = run.sh
endif

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f drum1_mesh.e && $(LN_S) $(srcdir)/drum1_mesh.e .
	-rm -f drum2_mesh.e && $(LN_S) $(srcdir)/drum2_mesh.e .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += drum*_mesh.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
