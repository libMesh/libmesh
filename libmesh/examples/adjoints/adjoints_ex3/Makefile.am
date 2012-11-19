include $(top_builddir)/Make.common.@METHOD@

AUTOMAKE_OPTIONS = subdir-objects 

AM_CPPFLAGS  = $(libmesh_CPPFLAGS)
AM_CFLAGS    = $(libmesh_CFLAGS)
AM_CXXFLAGS  = $(libmesh_CXXFLAGS)
AM_CPPFLAGS += $(libmesh_INCLUDE) -I$(top_builddir)/include
LIBS         = $(libmesh_LIBS)

if LIBMESH_ENABLE_FPARSER
  AM_CPPFLAGS += -I$(top_srcdir)/contrib/fparser
endif

check_PROGRAMS       = adjoints_ex3
check_SCRIPTS        = run.sh
adjoints_ex3dir      = $(examples_install_path)/adjoints/ex3
adjoints_ex3_DATA    = H-qoi.C H-qoi.h H_channel_quads.xda adjoints_ex3.C coupled_system.C coupled_system.h coupled_system.in domain.C domain.h femparameters.C femparameters.h general.in initial.C initial.h mysystems.h output.C output.h
adjoints_ex3_SOURCES = $(adjoints_ex3_DATA) run.sh
adjoints_ex3_LDADD   = $(top_builddir)/libmesh.la

CLEANFILES = primal.out.gmv.000.05 primal.out.gmv.000.04 primal.out.gmv.000.03 primal.out.gmv.000.02 primal.out.gmv.000.01 primal.out.gmv.000.00 out_timesteps.m out_time.m out_solvesteps.m out_clocktime.m out_activemesh.m error.gmv.004.00_with_pressure error.gmv.004.00_total error.gmv.004.00_non_pressure error.gmv.004.00_convection_diffusion_y error.gmv.004.00_convection_diffusion_x error.gmv.003.00_with_pressure error.gmv.003.00_total error.gmv.003.00_non_pressure error.gmv.003.00_convection_diffusion_y error.gmv.003.00_convection_diffusion_x error.gmv.002.00_with_pressure error.gmv.002.00_total error.gmv.002.00_non_pressure error.gmv.002.00_convection_diffusion_y error.gmv.002.00_convection_diffusion_x error.gmv.001.00_with_pressure error.gmv.001.00_total error.gmv.001.00_non_pressure error.gmv.001.00_convection_diffusion_y error.gmv.001.00_convection_diffusion_x error.gmv.000.00_with_pressure error.gmv.000.00_total error.gmv.000.00_non_pressure error.gmv.000.00_convection_diffusion_y error.gmv.000.00_convection_diffusion_x complete adjoint.out.gmv.000.04 adjoint.out.gmv.000.03 adjoint.out.gmv.000.02 adjoint.out.gmv.000.01 adjoint.out.gmv.000.00 

# also need links for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f general.in  && $(LN_S) $(srcdir)/general.in .
	-rm -f H_channel_quads.xda && $(LN_S) $(srcdir)/H_channel_quads.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in H_channel_quads.xda .linkstamp
endif

######################################################################
#
# Running the tests
#
TESTS_ENVIRONMENT = LIBMESH_DIR=$(abs_top_srcdir)
TESTS = $(check_SCRIPTS)


# Convenience method: 'make run'
run: $(check_PROGRAMS)
	LIBMESH_DIR=$(abs_top_srcdir) $(srcdir)/$(check_SCRIPTS)


# how to make documentation
CLEANFILES += $(check_PROGRAMS).php stdout.log
doc: $(check_PROGRAMS).php
html: $(check_PROGRAMS).php
$(check_PROGRAMS).php: $(check_PROGRAMS) \
                        $(top_srcdir)/contrib/bin/ex2html.sh \
                       Makefile
	$(MAKE) -s run > stdout.log
	$(AM_V_GEN) PATH=$(abs_top_srcdir)/contrib/bin:$(PATH) \
	  ex2html.sh $(check_PROGRAMS) $(abs_srcdir)

