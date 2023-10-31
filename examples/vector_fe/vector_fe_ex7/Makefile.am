example_name  = vector_fe_ex7
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/vector_fe/ex7
data          = div_grad_exact_solution.h solution_function.h vector_fe_ex7.C vector_fe_ex7.in run.sh
sources       = $(data)

CLEANFILES = out.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f vector_fe_ex7.in && $(LN_S) -f $(srcdir)/vector_fe_ex7.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += vector_fe_ex7.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
