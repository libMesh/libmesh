example_name  = reduced_basis_ex7
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/reduced_basis/ex7
data          = assembly.h horn.geo horn.msh rb_classes.h reduced_basis_ex7.C reduced_basis_ex7.in run.sh
sources       = $(data) run.sh

CLEANFILES =

clean-local:
	-rm -rf rb_data eim_data

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f reduced_basis_ex7.in && $(LN_S) -f $(srcdir)/reduced_basis_ex7.in .
	-rm -f horn.msh             && $(LN_S) -f $(srcdir)/horn.msh .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += reduced_basis_ex7.in horn.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
