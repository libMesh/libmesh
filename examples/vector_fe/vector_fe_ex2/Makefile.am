example_name  = vector_fe_ex2
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/vector_fe/ex2
data          = laplace.in laplace_system.C solution_function.h vector_fe_ex2.in laplace_exact_solution.h laplace_system.h vector_fe_ex2.C run.sh
sources       = $(data) run.sh

CLEANFILES = out.e out.inp

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f vector_fe_ex2.in && $(LN_S) $(srcdir)/vector_fe_ex2.in .
	-rm -f laplace.in && $(LN_S) $(srcdir)/laplace.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += vector_fe_ex2.in laplace.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
