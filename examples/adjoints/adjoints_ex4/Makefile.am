BUILT_SOURCES = femparameters.C # femparameters.h
DISTCLEANFILES = $(BUILT_SOURCES)

femparameters.C: $(top_srcdir)/examples/adjoints/adjoints_ex1/femparameters.C
	$(AM_V_GEN)$(LN_S) -f $< $@

#femparameters.h: $(top_srcdir)/examples/adjoints/adjoints_ex1/femparameters.h
#	$(AM_V_GEN)$(LN_S) -f $< .

example_name  = adjoints_ex4
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adjoints/ex4
data          = L-shaped.C element_qoi_derivative.C general.in side_postprocess.C L-shaped.h adjoints_ex4.C femparameters.C l-shaped.in side_qoi_derivative.C element_postprocess.C femparameters.h lshaped.xda run.sh
sources       = $(data) run.sh

CLEANFILES = primal.out.* adjoint_*.out.* primal.e* adjoint_0.e* adjoint_1.e*

# also need links for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES += .linkstamp
.linkstamp:
	-rm -f general.in  && $(LN_S) -f $(srcdir)/general.in .
	-rm -f lshaped.xda && $(LN_S) -f $(srcdir)/lshaped.xda .
	-rm -f l-shaped.in && $(LN_S) -f $(srcdir)/l-shaped.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in lshaped.xda l-shaped.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
