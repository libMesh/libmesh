// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic
// This file is part of rbOOmit.



// <h1>Reduced Basis: Example 5 - Reduced Basis Cantilever</h1>
// 3D cantilever beam using the Reduced Basis Method
// \author David Knezevic
// \author Kyunghoon "K" Lee
// \date 2012
//
// We consider four parameters in this problem:
//  x_scaling: scales the length of the cantilever
//  load_Fx: the traction in the x-direction on the right boundary of the cantilever
//  load_Fy: the traction in the y-direction on the right boundary of the cantilever
//  load_Fz: the traction in the z-direction on the right boundary of the cantilever

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/rb_data_serialization.h"
#include "libmesh/rb_data_deserialization.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/enum_solver_type.h"

// local includes
#include "rb_classes.h"
#include "assembly.h"

// boundary IDs
#define BOUNDARY_ID_MIN_Z 0
#define BOUNDARY_ID_MIN_Y 1
#define BOUNDARY_ID_MAX_X 2
#define BOUNDARY_ID_MAX_Y 3
#define BOUNDARY_ID_MIN_X 4
#define BOUNDARY_ID_MAX_Z 5
#define NODE_BOUNDARY_ID 10

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Define a function to scale the mesh according to the parameter.
void scale_mesh_and_plot(EquationSystems & es,
                         const RBParameters & mu,
                         const std::string & filename);

// Post-process the solution to compute stresses
void compute_stresses(EquationSystems & es);

// The main program.
int main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

#if !defined(LIBMESH_HAVE_XDR)
  // We need XDR support to write out reduced bases
  libmesh_example_requires(false, "--enable-xdr");
#elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
  // XDR binary support requires double precision
  libmesh_example_requires(false, "--disable-singleprecision");
#elif defined(LIBMESH_DEFAULT_TRIPLE_PRECISION)
  // I have no idea why long double isn't working here... [RHS]
  libmesh_example_requires(false, "double precision");
#endif

  // Eigen can take forever to solve the offline mode portion of this
  // example
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc or --enable-laspack");

  // Trilinos reports "true residual is too large" in the offline mode
  // portion of this example
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc or --enable-laspack");

  // This example only works if libMesh was compiled for 3D
  const unsigned int dim = 3;
  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  std::string parameters_filename = "reduced_basis_ex5.in";
  GetPot infile(parameters_filename);

  unsigned int n_elem_x = infile("n_elem_x", 0);
  unsigned int n_elem_y = infile("n_elem_y", 0);
  unsigned int n_elem_z = infile("n_elem_z", 0);
  Real x_size           = infile("x_size", 0.);
  Real y_size           = infile("y_size", 0.);
  Real z_size           = infile("z_size", 0.);

  bool store_basis_functions = infile("store_basis_functions", true);

  // Read the "online_mode" flag from the command line
  GetPot command_line(argc, argv);
  int online_mode = 0;
  if (command_line.search(1, "-online_mode"))
    online_mode = command_line.next(online_mode);


  Mesh mesh (init.comm(), dim);
  MeshTools::Generation::build_cube (mesh,
                                     n_elem_x,
                                     n_elem_y,
                                     n_elem_z,
                                     0., x_size,
                                     0., y_size,
                                     0., z_size,
                                     HEX8);

  // Let's add a node boundary condition so that we can impose a point
  // load.
  // Each processor should know about each boundary condition it can
  // see, so we loop over all elements, not just local elements.
  for (const auto & elem : mesh.element_ptr_range())
    {
      unsigned int side_max_x = 0, side_max_y = 0, side_max_z = 0;
      bool found_side_max_x = false, found_side_max_y = false, found_side_max_z = false;
      for (auto side : elem->side_index_range())
        {
          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_X))
            {
              side_max_x = side;
              found_side_max_x = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_Y))
            {
              side_max_y = side;
              found_side_max_y = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_Z))
            {
              side_max_z = side;
              found_side_max_z = true;
            }
        }

      // If elem has sides on boundaries
      // BOUNDARY_ID_MAX_X, BOUNDARY_ID_MAX_Y, BOUNDARY_ID_MAX_Z
      // then let's set a node boundary condition
      if (found_side_max_x && found_side_max_y && found_side_max_z)
        {
          for (auto n : elem->node_index_range())
            {
              if (elem->is_node_on_side(n, side_max_x) &&
                  elem->is_node_on_side(n, side_max_y) &&
                  elem->is_node_on_side(n, side_max_z))
                {
                  mesh.get_boundary_info().add_node(elem->node_ptr(n), NODE_BOUNDARY_ID);
                }
            }
        }
    }

  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems(mesh);

  // We override RBConstruction with ElasticityRBConstruction in order to
  // specialize a few functions for this particular problem.
  ElasticityRBConstruction & rb_con =
    equation_systems.add_system<ElasticityRBConstruction>("RBElasticity");

  // Also, initialize an ExplicitSystem to store stresses
  ExplicitSystem & stress_system =
    equation_systems.add_system<ExplicitSystem> ("StressSystem");
  stress_system.add_variable("sigma_00", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_01", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_02", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_10", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_11", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_12", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_20", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_21", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_22", CONSTANT, MONOMIAL);
  stress_system.add_variable("vonMises", CONSTANT, MONOMIAL);

  // Initialize the data structures for the equation system.
  equation_systems.init ();
  equation_systems.print_info();

  // Build a new RBEvaluation object which will be used to perform
  // Reduced Basis calculations. This is required in both the
  // "Offline" and "Online" stages.
  ElasticityRBEvaluation rb_eval(mesh.comm());

  // We need to give the RBConstruction object a pointer to
  // our RBEvaluation object
  rb_con.set_rb_evaluation(rb_eval);

  if (!online_mode) // Perform the Offline stage of the RB method
    {
      // Use CG solver.  This echoes the Petsc command line options set
      // in run.sh, but also works for non-PETSc linear solvers.
      rb_con.get_linear_solver()->set_solver_type(CG);

      // Read in the data that defines this problem from the specified text file
      rb_con.process_parameters_file(parameters_filename);

      // Print out info that describes the current setup of rb_con
      rb_con.print_info();

      // Prepare rb_con for the Construction stage of the RB method.
      // This sets up the necessary data structures and performs
      // initial assembly of the "truth" affine expansion of the PDE.
      rb_con.initialize_rb_construction();

      // Compute the reduced basis space by computing "snapshots", i.e.
      // "truth" solves, at well-chosen parameter values and employing
      // these snapshots as basis functions.
      rb_con.train_reduced_basis();

      // Write out the data that will subsequently be required for the Evaluation stage
#if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataSerialization::RBEvaluationSerialization rb_eval_writer(rb_con.get_rb_evaluation());
      rb_eval_writer.write_to_file("rb_eval.bin");
#else
      rb_con.get_rb_evaluation().legacy_write_offline_data_to_files();
#endif

      // If requested, write out the RB basis functions for visualization purposes
      if (store_basis_functions)
        {
          // Write out the basis functions
          rb_con.get_rb_evaluation().write_out_basis_functions(rb_con);
        }
    }
  else // Perform the Online stage of the RB method
    {
      // Read in the reduced basis data
#if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataDeserialization::RBEvaluationDeserialization rb_eval_reader(rb_eval);
      rb_eval_reader.read_from_file("rb_eval.bin", /*read_error_bound_data*/ true);
#else
      rb_eval.legacy_read_offline_data_from_files();
#endif

      // Initialize online parameters
      Real online_x_scaling = infile("online_x_scaling", 0.);
      Real online_load_Fx   = infile("online_load_Fx",   0.);
      Real online_load_Fy   = infile("online_load_Fy",   0.);
      Real online_load_Fz   = infile("online_load_Fz",   0.);
      Real online_point_load_Fx   = infile("online_point_load_Fx",   0.);
      Real online_point_load_Fy   = infile("online_point_load_Fy",   0.);
      Real online_point_load_Fz   = infile("online_point_load_Fz",   0.);
      RBParameters online_mu;
      online_mu.set_value("x_scaling", online_x_scaling);
      online_mu.set_value("load_Fx",   online_load_Fx);
      online_mu.set_value("load_Fy",   online_load_Fy);
      online_mu.set_value("load_Fz",   online_load_Fz);
      online_mu.set_value("point_load_Fx",   online_point_load_Fx);
      online_mu.set_value("point_load_Fy",   online_point_load_Fy);
      online_mu.set_value("point_load_Fz",   online_point_load_Fz);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();

      // Now do the Online solve using the precomputed reduced basis
      rb_eval.rb_solve(rb_eval.get_n_basis_functions());

      if (store_basis_functions)
        {
          // Read in the basis functions
          rb_eval.read_in_basis_functions(rb_con);

          // Plot the solution
          rb_con.load_rb_solution();

          const RBParameters & rb_eval_params = rb_eval.get_parameters();
          scale_mesh_and_plot(equation_systems, rb_eval_params, "RB_sol.e");
        }
    }

  return 0;
}

void scale_mesh_and_plot(EquationSystems & es,
                         const RBParameters & mu,
                         const std::string & filename)
{
  // Loop over the mesh nodes and move them!
  MeshBase & mesh = es.get_mesh();

  for (auto & node : mesh.node_ptr_range())
    (*node)(0) *= mu.get_value("x_scaling");

  // Post-process the solution to compute the stresses
  compute_stresses(es);

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_equation_systems (filename, es);
#endif

  // Loop over the mesh nodes and move them!
  for (auto & node : mesh.node_ptr_range())
    (*node)(0) /= mu.get_value("x_scaling");
}

void compute_stresses(EquationSystems & es)
{
  LOG_SCOPE("compute_stresses()", "main");

  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  ElasticityRBConstruction & system = es.get_system<ElasticityRBConstruction>("RBElasticity");

  unsigned int displacement_vars[3];
  displacement_vars[0] = system.variable_number ("u");
  displacement_vars[1] = system.variable_number ("v");
  displacement_vars[2] = system.variable_number ("w");
  const unsigned int u_var = system.variable_number ("u");

  const DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem & stress_system = es.get_system<ExplicitSystem>("StressSystem");
  const DofMap & stress_dof_map = stress_system.get_dof_map();
  unsigned int sigma_vars[3][3];
  sigma_vars[0][0] = stress_system.variable_number ("sigma_00");
  sigma_vars[0][1] = stress_system.variable_number ("sigma_01");
  sigma_vars[0][2] = stress_system.variable_number ("sigma_02");
  sigma_vars[1][0] = stress_system.variable_number ("sigma_10");
  sigma_vars[1][1] = stress_system.variable_number ("sigma_11");
  sigma_vars[1][2] = stress_system.variable_number ("sigma_12");
  sigma_vars[2][0] = stress_system.variable_number ("sigma_20");
  sigma_vars[2][1] = stress_system.variable_number ("sigma_21");
  sigma_vars[2][2] = stress_system.variable_number ("sigma_22");
  unsigned int vonMises_var = stress_system.variable_number ("vonMises");

  // Storage for the stress dof indices on each element
  std::vector<std::vector<dof_id_type>> dof_indices_var(system.n_vars());
  std::vector<dof_id_type> stress_dof_indices_var;

  // To store the stress tensor on each element
  DenseMatrix<Number> elem_sigma;

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      for (unsigned int var=0; var<3; var++)
        dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);

      fe->reinit (elem);

      elem_sigma.resize(3, 3);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        for (unsigned int C_i=0; C_i<3; C_i++)
          for (unsigned int C_j=0; C_j<3; C_j++)
            for (unsigned int C_k=0; C_k<3; C_k++)
              {
                const unsigned int n_var_dofs = dof_indices_var[C_k].size();

                // Get the gradient at this quadrature point
                Gradient displacement_gradient;
                for (unsigned int l=0; l<n_var_dofs; l++)
                  displacement_gradient.add_scaled(dphi[l][qp], system.current_solution(dof_indices_var[C_k][l]));

                for (unsigned int C_l=0; C_l<3; C_l++)
                  elem_sigma(C_i,C_j) += JxW[qp] * (elasticity_tensor(C_i, C_j, C_k, C_l) * displacement_gradient(C_l));
              }

      // Get the average stresses by dividing by the element volume
      elem_sigma.scale(1./elem->volume());

      // load elem_sigma data into stress_system
      for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
          {
            stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[i][j]);

            // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
            // one dof index per variable
            dof_id_type dof_index = stress_dof_indices_var[0];

            if ((stress_system.solution->first_local_index() <= dof_index) &&
                (dof_index < stress_system.solution->last_local_index()))
              {
                stress_system.solution->set(dof_index, elem_sigma(i,j));
              }

          }

      // Also, the von Mises stress
      Number vonMises_value = std::sqrt(0.5*(pow(elem_sigma(0,0) - elem_sigma(1,1),2.) +
                                             pow(elem_sigma(1,1) - elem_sigma(2,2),2.) +
                                             pow(elem_sigma(2,2) - elem_sigma(0,0),2.) +
                                             6.*(pow(elem_sigma(0,1),2.) + pow(elem_sigma(1,2),2.) + pow(elem_sigma(2,0),2.))
                                             ));
      stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);
      dof_id_type dof_index = stress_dof_indices_var[0];
      if ((stress_system.solution->first_local_index() <= dof_index) &&
          (dof_index < stress_system.solution->last_local_index()))
        {
          stress_system.solution->set(dof_index, vonMises_value);
        }

    }

  // Should call close and update when we set vector entries directly
  stress_system.solution->close();
  stress_system.update();
}
