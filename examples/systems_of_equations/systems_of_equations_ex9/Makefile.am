example_name  = systems_of_equations_ex9
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/systems_of_equations/ex9
data          = systems_of_equations_ex9.C run.sh systems_of_equations_ex9.exo
sources       = $(data) run.sh

CLEANFILES = solution.exo

if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f systems_of_equations_ex9.exo && $(LN_S) -f $(srcdir)/systems_of_equations_ex9.exo .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += systems_of_equations_ex9.exo .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
