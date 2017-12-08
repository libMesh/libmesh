// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// \author Robert Weidlich
// \date Copyright 2012

#include "libmesh/boundary_info.h"
#include "libmesh/diff_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/steady_solver.h"
#include "libmesh/transient_system.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

#include "nonlinear_neohooke_cc.h"
#include "solid_system.h"

// Solaris Studio has no NAN
#ifdef __SUNPRO_CC
#define NAN (1.0/0.0)
#endif

// Bring in everything from the libMesh namespace
using namespace libMesh;

SolidSystem::SolidSystem(EquationSystems & es,
                         const std::string & name_in,
                         const unsigned int number_in) :
  FEMSystem(es, name_in, number_in)
{

  // Add a time solver. We are just looking at a steady state problem here.
  this->time_solver = libmesh_make_unique<SteadySolver>(*this);
}

void SolidSystem::save_initial_mesh()
{
  System & aux_sys = this->get_equation_systems().get_system("auxiliary");

  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Loop over all nodes and copy the location from the current system to
  // the auxiliary system.
  for (const auto & node : this->get_mesh().local_node_ptr_range())
    for (unsigned int d = 0; d < dim; ++d)
      {
        unsigned int source_dof = node->dof_number(this->number(), var[d], 0);
        unsigned int dest_dof = node->dof_number(aux_sys.number(), undefo_var[d], 0);
        Number value = this->current_local_solution->el(source_dof);
        aux_sys.current_local_solution->set(dest_dof, value);
      }
}

void SolidSystem::init_data()
{
  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Get the default order of the used elements. Assumption:
  // Just one type of elements in the mesh.
  Order order = (*(this->get_mesh().elements_begin()))->default_order();

  // Add the node positions as primary variables.
  var[0] = this->add_variable("x", order);
  var[1] = this->add_variable("y", order);
  if (dim == 3)
    var[2] = this->add_variable("z", order);
  else
    var[2] = var[1];

  // Add variables for storing the initial mesh to an auxiliary system.
  System & aux_sys = this->get_equation_systems().get_system("auxiliary");
  undefo_var[0] = aux_sys.add_variable("undefo_x", order);
  undefo_var[1] = aux_sys.add_variable("undefo_y", order);
  undefo_var[2] = aux_sys.add_variable("undefo_z", order);

  // Set the time stepping options
  this->deltat = args("schedule/dt", 0.2);

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  //// Tell the system to march velocity forward in time, but
  //// leave p as a constraint only
  //this->time_evolving(var[0]);
  //this->time_evolving(var[1]);
  //if (dim == 3)
  //this->time_evolving(var[2]);

  // Tell the system which variables are containing the node positions
  set_mesh_system(this);

  this->set_mesh_x_var(var[0]);
  this->set_mesh_y_var(var[1]);
  if (dim == 3)
    this->set_mesh_z_var(var[2]);

  // Fill the variables with the position of the nodes
  this->mesh_position_get();

  System::reinit();

  // Set some options for the DiffSolver
  DiffSolver & solver = *(this->time_solver->diff_solver().get());
  solver.quiet = args("solver/quiet", false);
  solver.max_nonlinear_iterations = args("solver/nonlinear/max_nonlinear_iterations", 100);
  solver.relative_step_tolerance = args("solver/nonlinear/relative_step_tolerance", 1.e-3);
  solver.relative_residual_tolerance = args("solver/nonlinear/relative_residual_tolerance", 1.e-8);
  solver.absolute_residual_tolerance = args("solver/nonlinear/absolute_residual_tolerance", 1.e-8);
  solver.verbose = !args("solver/quiet", false);

  ((NewtonSolver &) solver).require_residual_reduction = args("solver/nonlinear/require_reduction", false);

  // And the linear solver options
  solver.max_linear_iterations = args("max_linear_iterations", 50000);
  solver.initial_linear_tolerance = args("initial_linear_tolerance", 1.e-3);
}

void SolidSystem::update()
{
  System::update();
  this->mesh_position_set();
}

void SolidSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // Pre-request all the data needed
  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);

  elem_fe->get_JxW();
  elem_fe->get_phi();
  elem_fe->get_dphi();
  elem_fe->get_xyz();

  FEBase * side_fe = libmesh_nullptr;
  c.get_side_fe(0, side_fe);

  side_fe->get_JxW();
  side_fe->get_phi();
  side_fe->get_xyz();
}

/**
 * Compute contribution from internal forces in elem_residual and contribution from
 * linearized internal forces (stiffness matrix) in elem_jacobian.
 */
bool SolidSystem::element_time_derivative(bool request_jacobian,
                                          DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * elem_fe = libmesh_nullptr;
  c.get_element_fe(0, elem_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // Element basis functions
  const std::vector<std::vector<RealGradient>> & dphi = elem_fe->get_dphi();

  // Dimension of the mesh
  const unsigned int dim = this->get_mesh().mesh_dimension();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.get_dof_indices(var[0]).size();
  libmesh_assert(n_u_dofs == c.get_dof_indices(var[1]).size());
  if (dim == 3)
    libmesh_assert(n_u_dofs ==  c.get_dof_indices(var[2]).size());

  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Some matrices and vectors for storing the results of the constitutive
  // law
  DenseMatrix<Real> stiff;
  DenseVector<Real> res;
  VectorValue<Gradient> grad_u;

  // Instantiate the constitutive law
  // Just calculate jacobian contribution when we need to
  NonlinearNeoHookeCurrentConfig material(dphi, args, request_jacobian);

  // Get a reference to the auxiliary system
  TransientExplicitSystem & aux_system = this->get_equation_systems().get_system<
    TransientExplicitSystem>("auxiliary");
  std::vector<dof_id_type> undefo_index;

  // Assume symmetry of local stiffness matrices
  bool use_symmetry = args("assembly/use_symmetry", false);

  // Now we will build the element Jacobian and residual.
  // This must be calculated at each quadrature point by
  // summing the solution degree-of-freedom values by
  // the appropriate weight functions.
  // This class just takes care of the assembly. The matrix of
  // the jacobian and the residual vector are provided by the
  // constitutive formulation.

  for (unsigned int qp = 0; qp != n_qpoints; qp++)
    {
      // Compute the displacement gradient
      grad_u(0) = grad_u(1) = grad_u(2) = 0;
      for (unsigned int d = 0; d < dim; ++d)
        {
          std::vector<Number> u_undefo;
          aux_system.get_dof_map().dof_indices(&c.get_elem(), undefo_index, undefo_var[d]);
          aux_system.current_local_solution->get(undefo_index, u_undefo);
          for (unsigned int l = 0; l != n_u_dofs; l++)
            grad_u(d).add_scaled(dphi[l][qp], u_undefo[l]); // u_current(l)); // -
        }

      // initialize the constitutive formulation with the current displacement
      // gradient
      material.init_for_qp(grad_u, qp);

      // Acquire, scale and assemble residual and stiffness
      for (unsigned int i = 0; i < n_u_dofs; i++)
        {
          res.resize(dim);
          material.get_residual(res, i);
          res.scale(JxW[qp]);
          for (unsigned int ii = 0; ii < dim; ++ii)
            (c.get_elem_residual(ii))(i) += res(ii);

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert(c.elem_solution_derivative == 1.0);
              for (unsigned int j = (use_symmetry ? i : 0); j < n_u_dofs; j++)
                {
                  material.get_linearized_stiffness(stiff, i, j);
                  stiff.scale(JxW[qp]);
                  for (unsigned int ii = 0; ii < dim; ++ii)
                    {
                      for (unsigned int jj = 0; jj < dim; ++jj)
                        {
                          (c.get_elem_jacobian(ii,jj))(i, j) += stiff(ii, jj);
                          if (use_symmetry && i != j)
                            (c.get_elem_jacobian(ii,jj))(j, i) += stiff(jj, ii);
                        }
                    }
                }
            }
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

bool SolidSystem::side_time_derivative(bool request_jacobian,
                                       DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  // Apply displacement boundary conditions with penalty method

  // Get the current load step
  Real ratio = this->get_equation_systems().parameters.get<Real>("progress") + 0.001;

  // The BC are stored in the simulation parameters as array containing sequences of
  // four numbers: Id of the side for the displacements and three values describing the
  // displacement. E.g.: bc/displacement = '5 nan nan -1.0'. This will move all nodes of
  // side 5 about 1.0 units down the z-axis while leaving all other directions unrestricted

  // Get number of BCs to enforce
  unsigned int num_bc = args.vector_variable_size("bc/displacement");
  if (num_bc % 4 != 0)
    libmesh_error_msg("ERROR, Odd number of values in displacement boundary condition.");
  num_bc /= 4;

  // Loop over all BCs
  for (unsigned int nbc = 0; nbc < num_bc; nbc++)
    {
      // Get IDs of the side for this BC
      short int positive_boundary_id = args("bc/displacement", 1, nbc * 4);

      // The current side may not be on the boundary to be restricted
      if (!this->get_mesh().get_boundary_info().has_boundary_id
          (&c.get_elem(),c.get_side(),positive_boundary_id))
        continue;

      // Read values from configuration file
      Point diff_value;
      for (unsigned int d = 0; d < c.get_dim(); ++d)
        diff_value(d) = args("bc/displacement", NAN, nbc * 4 + 1 + d);

      // Scale according to current load step
      diff_value *= ratio;

      Real penalty_number = args("bc/displacement_penalty", 1e7);

      FEBase * fe = libmesh_nullptr;
      c.get_side_fe(0, fe);

      const std::vector<std::vector<Real>> & phi = fe->get_phi();
      const std::vector<Real> & JxW = fe->get_JxW();
      const std::vector<Point> & coords = fe->get_xyz();

      unsigned int n_x_dofs = c.get_dof_indices(this->var[0]).size();

      // get mappings for dofs for auxiliary system for original mesh positions
      const System & auxsys = this->get_equation_systems().get_system("auxiliary");
      const DofMap & auxmap = auxsys.get_dof_map();
      std::vector<dof_id_type> undefo_dofs[3];
      for (unsigned int d = 0; d < c.get_dim(); ++d)
        auxmap.dof_indices(&c.get_elem(), undefo_dofs[d], undefo_var[d]);

      for (unsigned int qp = 0; qp < c.get_side_qrule().n_points(); ++qp)
        {
          // calculate coordinates of qp on undeformed mesh
          Point orig_point;
          for (unsigned int i = 0; i < n_x_dofs; ++i)
            {
              for (unsigned int d = 0; d < c.get_dim(); ++d)
                {
                  Number orig_val = auxsys.current_solution(undefo_dofs[d][i]);

#if LIBMESH_USE_COMPLEX_NUMBERS
                  orig_point(d) += phi[i][qp] * orig_val.real();
#else
                  orig_point(d) += phi[i][qp] * orig_val;
#endif
                }
            }

          // Calculate displacement to be enforced.
          Point diff = coords[qp] - orig_point - diff_value;

          // Assemble
          for (unsigned int i = 0; i < n_x_dofs; ++i)
            {
              for (unsigned int d1 = 0; d1 < c.get_dim(); ++d1)
                {
                  if (libmesh_isnan(diff(d1)))
                    continue;
                  Real val = JxW[qp] * phi[i][qp] * diff(d1) * penalty_number;
                  (c.get_elem_residual(var[d1]))(i) += val;
                }
              if (request_jacobian)
                {
                  for (unsigned int j = 0; j < n_x_dofs; ++j)
                    {
                      for (unsigned int d1 = 0; d1 < c.get_dim(); ++d1)
                        {
                          if (libmesh_isnan(diff(d1)))
                            continue;
                          Real val = JxW[qp] * phi[i][qp] * phi[j][qp] * penalty_number;
                          (c.get_elem_jacobian(var[d1],var[d1]))(i, j) += val;
                        }
                    }
                }
            }
        }
    }

  return request_jacobian;
}
