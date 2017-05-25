BUILT_SOURCES = femparameters.C
DISTCLEANFILES = $(BUILT_SOURCES)

femparameters.C: $(top_srcdir)/examples/adjoints/adjoints_ex1/femparameters.C
	$(AM_V_GEN)$(LN_S) -f $< $@

example_name  = adjoints_ex6
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adjoints/ex6
data          = adjoints_ex6.C element_postprocess.C femparameters.C poisson.C general.in Makefile.in poisson.in femparameters.h poisson.h
sources       = $(data) run.sh

CLEANFILES = *.gmv.* primal.e* dual.e*

# also need links for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f general.in  && $(LN_S) -f $(srcdir)/general.in .
	-rm -f poisson.in && $(LN_S) -f $(srcdir)/poisson.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in poisson.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
