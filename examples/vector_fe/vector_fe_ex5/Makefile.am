example_name  = vector_fe_ex5
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/vector_fe/ex5
data          = vector_fe_ex5.C exact_solution.C assembly.C vector_fe_ex5.in run.sh
sources       = $(data) run.sh

CLEANFILES = out.e

if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -rf vector_fe_ex5.in && $(LN_S) -f $(srcdir)/vector_fe_ex5.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += vector_fe_ex5.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
