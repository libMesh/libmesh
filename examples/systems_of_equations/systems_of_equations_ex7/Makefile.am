example_name  = systems_of_equations_ex7
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/systems_of_equations/ex7
data          = systems_of_equations_ex7.C run.sh systems_of_equations_ex7.in
sources       = $(data) run.sh

CLEANFILES = solution_0.exo

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f systems_of_equations_ex7.in && $(LN_S) -f $(srcdir)/systems_of_equations_ex7.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += systems_of_equations_ex7.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
