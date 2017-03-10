// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/first_order_unsteady_solver.h"
#include "libmesh/diff_system.h"
#include "libmesh/quadrature.h"

namespace libMesh
{

void FirstOrderUnsteadySolver::prepare_accel(DiffContext & context)
{
  context.get_elem_solution_accel() = context.get_elem_solution_rate();

  context.elem_solution_accel_derivative = context.get_elem_solution_rate_derivative();
}

bool FirstOrderUnsteadySolver::compute_second_order_eqns(bool compute_jacobian, DiffContext & c)
{
  FEMContext & context = cast_ref<FEMContext &>(c);

  unsigned int n_qpoints = context.get_element_qrule().n_points();

  for (unsigned int var = 0; var != context.n_vars(); ++var)
    {
      if (!this->_system.is_second_order_var(var))
        continue;

      unsigned int dot_var = this->_system.get_second_order_dot_var(var);

      // We're assuming that the FE space for var and dot_var are the same
      libmesh_assert( context.get_system().variable(var).type() ==
                      context.get_system().variable(dot_var).type() );

      FEBase * elem_fe = libmesh_nullptr;
      context.get_element_fe( var, elem_fe );

      const std::vector<Real> & JxW = elem_fe->get_JxW();

      const std::vector<std::vector<Real> > & phi = elem_fe->get_phi();

      const unsigned int n_dofs = cast_int<unsigned int>
        (context.get_dof_indices(dot_var).size());

      DenseSubVector<Number> & Fu = context.get_elem_residual(var);
      DenseSubMatrix<Number> & Kuu = context.get_elem_jacobian( var, var );
      DenseSubMatrix<Number> & Kuv = context.get_elem_jacobian( var, dot_var );

      for (unsigned int qp = 0; qp != n_qpoints; ++qp)
        {
          Number udot, v;
          context.interior_rate(var, qp, udot);
          context.interior_value(dot_var, qp, v);

          for (unsigned int i = 0; i < n_dofs; i++)
            {
              Fu(i) += JxW[qp]*(udot-v)*phi[i][qp];

              if (compute_jacobian)
                {
                  Number rate_factor = JxW[qp]*context.get_elem_solution_rate_derivative()*phi[i][qp];
                  Number soln_factor = JxW[qp]*context.get_elem_solution_derivative()*phi[i][qp];

                  Kuu(i,i) += rate_factor*phi[i][qp];
                  Kuv(i,i) -= soln_factor*phi[i][qp];

                  for (unsigned int j = i+1; j < n_dofs; j++)
                    {
                      Kuu(i,j) += rate_factor*phi[j][qp];
                      Kuu(j,i) += rate_factor*phi[j][qp];

                      Kuv(i,j) -= soln_factor*phi[j][qp];
                      Kuv(j,i) -= soln_factor*phi[j][qp];
                    }
                }
            }
        }
    }

  return compute_jacobian;
}

} // namespace libMesh
