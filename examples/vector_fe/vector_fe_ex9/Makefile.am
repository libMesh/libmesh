example_name  = vector_fe_ex9
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/vector_fe/ex9
data          = vector_fe_ex9.C vector_fe_ex9.in hdg_problem.h hdg_problem.C exact_soln.h run.sh
sources       = $(data)

CLEANFILES = out.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f vector_fe_ex9.in && $(LN_S) -f $(srcdir)/vector_fe_ex9.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += vector_fe_ex9.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
