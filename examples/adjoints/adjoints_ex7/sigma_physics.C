// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/getpot.h"

#include "sigma_physics.h"

#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/zero_function.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "heatsystem.h"

void SigmaPhysics::init_data (System & sys)
{
  T_var = dynamic_cast<FEMSystem&>(sys).add_variable ("T", static_cast<Order>(1), LAGRANGE);

  GetPot infile("heat.in");
  _k = infile("k", 1.0);
  _analytic_jacobians = infile("analytic_jacobians", true);

  // The temperature is evolving, with a first order time derivative
  this->time_evolving(T_var, 1);

  const boundary_id_type all_ids[4] = {0, 1, 2, 3};
  std::set<boundary_id_type> all_bdys(all_ids, all_ids+(4));

  std::vector<unsigned int> T_only(1, T_var);

  ZeroFunction<Number> zero;

  ConstFunction<Number> one(1.0);

  // Dirichlet boundary conditions for the primal.
  sys.get_dof_map().add_dirichlet_boundary(DirichletBoundary (all_bdys, T_only, &one));

  // Set adjoint boundary conditions.
  sys.get_dof_map().add_adjoint_dirichlet_boundary(DirichletBoundary (all_bdys, T_only, &zero), 0);
  sys.get_dof_map().add_adjoint_dirichlet_boundary(DirichletBoundary (all_bdys, T_only, &zero), 1);

}

void SigmaPhysics::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * elem_fe = nullptr;
  c.get_element_fe(T_var, elem_fe);

  // Now make sure we have requested all the data
  // we need to build the linear system.
  elem_fe->get_JxW();
  elem_fe->get_dphi();
  elem_fe->get_phi();
  elem_fe->get_xyz();

  // Don't waste time on side computations for T
  FEBase * side_fe = nullptr;
  c.get_side_fe(T_var, side_fe);
  side_fe->get_nothing();

  // We'll have a more automatic solution to preparing adjoint
  // solutions for time integration, eventually...
  if (c.is_adjoint())
    {
      // A reference to the system context is built with
      const System & sys = c.get_system();

      // Get a pointer to the adjoint solution vector
      NumericVector<Number> & adjoint_solution0 =
        const_cast<System &>(sys).get_adjoint_solution(0);

      // Add this adjoint solution to the vectors that diff context should localize
      c.add_localized_vector(adjoint_solution0, sys);

      // Get a pointer to the adjoint solution vector
      NumericVector<Number> & adjoint_solution1 =
        const_cast<System &>(sys).get_adjoint_solution(1);

      // Add this adjoint solution to the vectors that diff context should localize
      c.add_localized_vector(adjoint_solution1, sys);
    }

}

bool SigmaPhysics::element_time_derivative (bool request_jacobian,
                                          DiffContext & context)
{
  bool compute_jacobian = request_jacobian && _analytic_jacobians;

  FEMContext & c = cast_ref<FEMContext &>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase * elem_fe = nullptr;
  c.get_element_fe(T_var, elem_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> & JxW = elem_fe->get_JxW();

  // Element basis functions
  const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = elem_fe->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.n_dof_indices(T_var);

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> & K = c.get_elem_jacobian(T_var, T_var);
  DenseSubVector<Number> & F = c.get_elem_residual(T_var);

  // Quadrature point locations
  const std::vector<Point > & q_point = elem_fe->get_xyz();

  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Conductivity
  Real sigma = 1.0;

  // Forcing function
  Real f = 1.0;

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution gradient at the Newton iterate
      Gradient grad_T = c.interior_gradient(T_var, qp);

      // Location of the current qp
      const Real x = q_point[qp](0);

      // Spatially varying conductivity
      sigma = 0.001 + x;

      for (unsigned int i=0; i != n_u_dofs; i++)
      {
F(i) += JxW[qp] * ( ( -sigma * (grad_T * dphi[i][qp]) ) + (f * phi[i][qp]) );
      }

      if (compute_jacobian)
        for (unsigned int i=0; i != n_u_dofs; i++)
          for (unsigned int j=0; j != n_u_dofs; ++j)
            K(i,j) += JxW[qp] * -sigma * (dphi[i][qp] * dphi[j][qp]);
    } // end of the quadrature point qp-loop

  return compute_jacobian;
}

std::unique_ptr<DifferentiablePhysics> SigmaPhysics::clone_physics()
{
  return libmesh_make_unique<SigmaPhysics>(*this);
}
