example_name  = vector_fe_ex10
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/vector_fe/ex10
data          = grad_div_exact_solution.h solution_function.h vector_fe_ex10.C vector_fe_ex10.in run.sh
sources       = $(data)

CLEANFILES = out.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f vector_fe_ex10.in && $(LN_S) -f $(srcdir)/vector_fe_ex10.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += vector_fe_ex10.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
