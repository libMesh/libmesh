// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// <h1>Vector Finite Elements Example 9 - Hybridizable Discontinuous Galerkin Navier Stokes</h1>
// \author Alexander Lindsay
// \date 2023

// Basic utilities.
#include "libmesh/string_to_enum.h"

// The solver packages supported by libMesh.
#include "libmesh/enum_solver_package.h"

// The mesh object and mesh generation and modification utilities.
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"

// Matrix and vector types.
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/numeric_vector.h"

// The finite element object and the geometric element type.
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"

// Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"
#include "libmesh/enum_quadrature_type.h"

// The dof map, which handles degree of freedom indexing.
#include "libmesh/dof_map.h"

// The system of equations.
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/static_condensation.h"

// The exact solution and error computation.
#include "libmesh/exact_solution.h"
#include "libmesh/enum_norm_type.h"
#include "exact_soln.h"

// I/O utilities.
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"

// The HDGProblem application context
#include "hdg_problem.h"

#if defined(LIBMESH_HAVE_EIGEN_DENSE) && defined(LIBMESH_HAVE_PETSC)

using namespace libMesh;

int
main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init(argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Parse the input file.
  GetPot infile("vector_fe_ex9.in");

  // But allow the command line to override it.
  infile.parse_command_line(argc, argv);

  // Read in parameters from the command line and the input file.
  const unsigned int dimension = 2;
  const unsigned int grid_size = infile("grid_size", 2);
  const bool mms = infile("mms", true);
  const Real nu = infile("nu", 1.);
  const bool cavity = infile("cavity", false);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build.
  libmesh_example_requires(dimension <= LIBMESH_DIM, dimension << "D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the cube [-1,1]^D. To accomodate first order side hierarchics, we must
  // use TRI6/7 elements
  const std::string elem_str = infile("element_type", std::string("TRI6"));

  libmesh_error_msg_if(elem_str != "TRI6" && elem_str != "TRI7",
                       "You selected "
                           << elem_str
                           << " but this example must be run with TRI6, TRI7, QUAD8, or QUAD9 in 2d"
                           << " or with TET14, or HEX27 in 3d.");

  if (mms && !cavity)
    MeshTools::Generation::build_square(
        mesh, grid_size, grid_size, 0., 2., -1, 1., Utility::string_to_enum<ElemType>(elem_str));
  else if (cavity)
    MeshTools::Generation::build_square(
        mesh, grid_size, grid_size, -1, 1, -1, 1., Utility::string_to_enum<ElemType>(elem_str));
  else
    MeshTools::Generation::build_square(mesh,
                                        5 * grid_size,
                                        grid_size,
                                        0.,
                                        10.,
                                        -1,
                                        1.,
                                        Utility::string_to_enum<ElemType>(elem_str));

  // Create an equation systems object.
  EquationSystems equation_systems(mesh);

  // Declare the system "Mixed" and its variables.
  auto & system = equation_systems.add_system<NonlinearImplicitSystem>("HDG");

  // Adds the velocity variables and their gradients
  system.add_variable("qu", FIRST, L2_LAGRANGE_VEC);
  system.add_variable("qv", FIRST, L2_LAGRANGE_VEC);
  system.add_variable("vel_x", FIRST, L2_LAGRANGE);
  system.add_variable("vel_y", FIRST, L2_LAGRANGE);

  // Add our Lagrange multiplier to the implicit system
  system.add_variable("lm_u", FIRST, SIDE_HIERARCHIC);
  system.add_variable("lm_v", FIRST, SIDE_HIERARCHIC);
  const auto p_num = system.add_variable("pressure", FIRST, L2_LAGRANGE);
  if (cavity)
    system.add_variable("global_lm", FIRST, SCALAR);

  const FEType vector_fe_type(FIRST, L2_LAGRANGE_VEC);
  const FEType scalar_fe_type(FIRST, L2_LAGRANGE);
  const FEType lm_fe_type(FIRST, SIDE_HIERARCHIC);

  auto & sc = system.get_dof_map().add_static_condensation();
  sc.dont_condense_vars({p_num});

  HDGProblem hdg(nu, cavity);
  hdg.mesh = &mesh;
  hdg.system = &system;
  hdg.dof_map = &system.get_dof_map();
  hdg.vector_fe = FEVectorBase::build(dimension, vector_fe_type);
  hdg.scalar_fe = FEBase::build(dimension, scalar_fe_type);
  hdg.qrule = QBase::build(QGAUSS, dimension, FIFTH);
  hdg.qface = QBase::build(QGAUSS, dimension - 1, FIFTH);
  hdg.vector_fe_face = FEVectorBase::build(dimension, vector_fe_type);
  hdg.scalar_fe_face = FEBase::build(dimension, scalar_fe_type);
  hdg.lm_fe_face = FEBase::build(dimension, lm_fe_type);
  hdg.mms = mms;
  hdg.sc = &sc;

  system.nonlinear_solver->residual_object = &hdg;
  system.nonlinear_solver->jacobian_object = &hdg;
  system.nonlinear_solver->attach_preconditioner(&sc);

  hdg.init();

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Solve the implicit system for the Lagrange multiplier
  system.solve();

  if (mms)
  {
    //
    // Now we will compute our solution approximation errors
    //

    equation_systems.parameters.set<const ExactSoln *>("vel_x_exact_sol") = &hdg.u_true_soln;
    equation_systems.parameters.set<const ExactSoln *>("vel_y_exact_sol") = &hdg.v_true_soln;
    equation_systems.parameters.set<const ExactSoln *>("pressure_exact_sol") = &hdg.p_true_soln;
    ExactSolution exact_sol(equation_systems);
    exact_sol.attach_exact_value(&compute_error);

    // Compute the error.
    exact_sol.compute_error("HDG", "vel_x");
    exact_sol.compute_error("HDG", "vel_y");
    exact_sol.compute_error("HDG", "pressure");

    // Print out the error values.
    libMesh::out << "L2 error for u is: " << exact_sol.l2_error("HDG", "vel_x") << std::endl;
    libMesh::out << "L2 error for v is: " << exact_sol.l2_error("HDG", "vel_y") << std::endl;
    libMesh::out << "L2 error for pressure is: " << exact_sol.l2_error("HDG", "pressure")
                 << std::endl;
  }

#ifdef LIBMESH_HAVE_EXODUS_API

  // We write the file in the ExodusII format.
  ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}

#else

int
main()
{
  return 0;
}

#endif
