BUILT_SOURCES = femparameters.C # femparameters.h
DISTCLEANFILES = $(BUILT_SOURCES)

femparameters.C: $(top_srcdir)/examples/adjoints/adjoints_ex1/femparameters.C
	$(AM_V_GEN)$(LN_S) -f $< $@

#femparameters.h: $(top_srcdir)/examples/adjoints/adjoints_ex1/femparameters.h
#	$(AM_V_GEN)$(LN_S) -f $< .

example_name  = adjoints_ex5
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adjoints/ex5
data          = heatsystem.C element_qoi_derivative.C factoryfunction.C general.in initial.C heatsystem.h adjoints_ex5.C femparameters.C heat.in femparameters.h adjoint_initial.C adjoint_initial.h initial.h run.sh
sources       = $(data) run.sh

CLEANFILES = *.gmv.* primal.e* dual.e*

# also need links for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f general.in  && $(LN_S) -f $(srcdir)/general.in .
	-rm -f heat.in && $(LN_S) -f $(srcdir)/heat.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in heat.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
