// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1>Miscellaneous Example 7 - The PetscDMNonlinearSolver and
// Variational Inequality (VI) problems</h1>
// \author Dmitry Karpeyev
// \date 2012
//
// In this example, LibMesh interfaces directly with PETSc's
// variational inequality solver through PetscDMNonlinearSolver
// (available in PETSc-3.3.0 or above).

// Example include files
#include "biharmonic.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Print usage information if requested on command line
void print_help(int argc, char ** argv);

int main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);
  if (on_command_line("--help"))
    print_help(argc, argv);
  else
    {
#if !defined(LIBMESH_ENABLE_SECOND_DERIVATIVES)
      libmesh_example_requires(false, "--enable-second");
#elif !defined(LIBMESH_ENABLE_PERIODIC)
      libmesh_example_requires(false, "--enable-periodic");
#endif

      // This is a PETSc-specific solver
      libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

      const int dim = command_line_value("dim", 1);

      // Skip higher-dimensional examples on a lower-dimensional libMesh build
      libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

      // This example only works with ReplicatedMesh.
      ReplicatedMesh mesh(init.comm());
      Biharmonic biharmonic(mesh);
      biharmonic.viewParameters();
      biharmonic.init();
      biharmonic.run();
    }
  return 0;
}





void print_help(int, char ** argv)
{
  libMesh::out << "This example solves the Cahn-Hillard equation with chemical potential f:\n"
               << "    u_t = \\div(M(u)\\grad f(u))\n"
               << "Here we have\n"
               << "    u, -1 <= u <= 1        -- relative concentration (difference of two concentrations in a binary mixture) \n"
               << "    M, M >= 0              -- mobility of the mixture\n"
               << "    f = \\delta E/\\delta u  -- variational derivative of the free energy functional E\n"
               << "    E = \\int[\\kappa/2 |\\grac u|^ + g(u)]\n"
               << "where the gradient term is the interfacial energy density with \\kappa quantifying the energy of the interface,\n"
               << "and g(u) is the bulk energy density\n"
               << "    g(u) = \\theta L(u) + \\theta_c W(u),\n"
               << "L(u) is the (optional, in this model) logarithmic term corresponding to the entropy of the mixture:\n"
               << "    L(u) = (\\theta/2)[(1+u)\\ln((1+u)/2) + (1-u)\\ln((1-u)/2)],\n"
               << "where \\theta is related to the Boltzmann factor k_B T - a proxy for the absolute temperature T.\n"
               << "L can be optionally approximated ('truncated') using a quadratic or a cubic polynomial on [-1,1]\n"
               << "W(u) is the (optional, in this model) potential promoting demixing.  It can take the form of \n"
               << "a 'double well' potential\n"
               << "    W(u) = \\theta_c (u^4/4 - u^2/2),\n"
               << "         or \n"
               << "a 'double obstacle' potential\n"
               << "    W(u) = (\\theta_c/2)(1-u^2),\n"
               << "where \\theta_c is the critical 'temperature'.\n"
               << "Finally, mobility M can be constant of 'degenerate', by which we mean that M is varying with u and \n"
               << "vanishing (degenerating) whenever u reaches critical values +1 or -1:\n"
               << "    M(u) = 1.0\n"
               << "      or\n"
               << "    M(u) = (1.0 - u^2)\n"
               << "Degenerate mobility should generally be used only in conjunction with logarithmic free energy terms.\n\n"
               << "The equation is solved on a periodic domain (in 1D, 2D or 3D)\n"
               << "using a Galerkin formulation with C^1 elements approximating the H^2_{per} function space.\n\n"
               << "\n-----------\n"
               << "COMPILING: "
               << "\n-----------\n"
               << "Compile as follows (assuming libmesh has been built): \n"
               << "METHOD=<method> make \n"
               << "where <method> is dbg or opt.\n"
               << "\n-----------\n"
               << "HELP:        "
               << "\n-----------\n"
               << "Print this help message:\n"
               << argv[0] << " --help\n"
               << "\n-----------\n"
               << "RUNNING:     "
               << "\n-----------\n"
               << "Run in serial with PETSc-3.4 and above as follows:\n"
               << "\n"
               << argv[0] << "\n"
               << "               [--verbose] dim=<1|2|3> N=<number_of_linear_elements> \n"
               << "               kappa=<kappa_value> growth=<yes|no> degenerate=<yes|no> [--cahn-hillard]                                           \n"
               << "               [--netforce]  energy=<double_well|double_obstacle|log_double_well|log_double_obstacle>  log_truncation_order=<2|3> \n"
               << "               theta=<theta_value> theta_c=<theta_c_value>                                                                        \n"
               << "               initial_state=<ball|rod|strip> initial_center='x [y [z]]' initial_width=<width>                                    \n"
               << "               min_time=<initial_time> max_time=<final_time> dt=<timestep_size> crank_nicholson_weight=<between_0_and_1>          \n"
               << "               output_base=<base_filename> output_dt=<output_timestep_size> [-snes_type vinewtonrsls | vinewtonssls]              \n"
               << "(with PETSc-3.3 use -snes_type virs | viss).\n"
               << argv[0] << " --verbose \n"
               << "is a pretty good start.\n"
               << "\nModeling a 1D system with 2D or 3D (for a strip the second and third components of the center are immaterial):\n"
               << argv[0]<< " --verbose dim=1 N=1024 initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-6\n"
               << argv[0]<< " --verbose dim=2 N=64   initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
               << argv[0]<< " --verbose dim=3 N=32   initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
               << "\n"
               << "Modeling a 2D system with 3D (for a rod the third component of the center is immaterial) \n"
               << argv[0]<< " --verbose dim=2 N=64   initial_state=rod initial_center='0.5 0.5' initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
               << argv[0]<< " --verbose dim=3 N=32   initial_state=rod initial_center='0.5 0.5' initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
               << "\n"
               << "A 3D system with an initial ball in the center\n"
               << argv[0] << " --verbose dim=3 N=32   initial_state=ball initial_center='0.5 0.5 0.5' initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
               << "\n"
               << "Add -snes_type vinewtonrsls to run the variational inequality version that ensures the solution is between -1.0 and 1.0 at all times.\n"
               << "If using PETSc-3.3, run with -snes_type virs.\n"
               << "\n"
               << std::endl;
}
