example_name  = adjoints_ex3
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/adjoints/ex3
data          = H-qoi.C H-qoi.h H_channel_quads.xda adjoints_ex3.C coupled_system.C coupled_system.h coupled_system.in domain.C domain.h femparameters.C femparameters.h general.in initial.C initial.h mysystems.h output.C output.h run.sh
sources       = $(data) run.sh

CLEANFILES = primal.out.gmv.000.05 primal.out.gmv.000.04 primal.out.gmv.000.03 primal.out.gmv.000.02 primal.out.gmv.000.01 primal.out.gmv.000.00 out_timesteps.m out_time.m out_solvesteps.m out_clocktime.m out_activemesh.m error.gmv.004.00_with_pressure error.gmv.004.00_total error.gmv.004.00_non_pressure error.gmv.004.00_convection_diffusion_y error.gmv.004.00_convection_diffusion_x error.gmv.003.00_with_pressure error.gmv.003.00_total error.gmv.003.00_non_pressure error.gmv.003.00_convection_diffusion_y error.gmv.003.00_convection_diffusion_x error.gmv.002.00_with_pressure error.gmv.002.00_total error.gmv.002.00_non_pressure error.gmv.002.00_convection_diffusion_y error.gmv.002.00_convection_diffusion_x error.gmv.001.00_with_pressure error.gmv.001.00_total error.gmv.001.00_non_pressure error.gmv.001.00_convection_diffusion_y error.gmv.001.00_convection_diffusion_x error.gmv.000.00_with_pressure error.gmv.000.00_total error.gmv.000.00_non_pressure error.gmv.000.00_convection_diffusion_y error.gmv.000.00_convection_diffusion_x complete adjoint.out.gmv.000.04 adjoint.out.gmv.000.03 adjoint.out.gmv.000.02 adjoint.out.gmv.000.01 adjoint.out.gmv.000.00

# also need links for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f general.in          && $(LN_S) $(srcdir)/general.in .
	-rm -f coupled_system.in   && $(LN_S) $(srcdir)/coupled_system.in .
	-rm -f H_channel_quads.xda && $(LN_S) $(srcdir)/H_channel_quads.xda .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += general.in coupled_system.in H_channel_quads.xda .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
