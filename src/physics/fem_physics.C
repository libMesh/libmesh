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



#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/time_solver.h"
#include "libmesh/unsteady_solver.h" // For eulerian_residual


namespace libMesh
{

bool FEMPhysics::eulerian_residual (bool request_jacobian,
                                    DiffContext & /*c*/)
{
  // Only calculate a mesh movement residual if it's necessary
  if (!_mesh_sys)
    return request_jacobian;

  libmesh_not_implemented();

#if 0
  FEMContext & context = cast_ref<FEMContext &>(c);

  // This function only supports fully coupled mesh motion for now
  libmesh_assert_equal_to (_mesh_sys, this);

  unsigned int n_qpoints = (context.get_element_qrule())->n_points();

  const unsigned int n_x_dofs = (_mesh_x_var == libMesh::invalid_uint) ?
    0 : context.dof_indices_var[_mesh_x_var].size();
  const unsigned int n_y_dofs = (_mesh_y_var == libMesh::invalid_uint) ?
    0 : context.dof_indices_var[_mesh_y_var].size();
  const unsigned int n_z_dofs = (_mesh_z_var == libMesh::invalid_uint) ?
    0 : context.dof_indices_var[_mesh_z_var].size();

  const unsigned int mesh_xyz_var = n_x_dofs ? _mesh_x_var :
    (n_y_dofs ? _mesh_y_var :
     (n_z_dofs ? _mesh_z_var :
      libMesh::invalid_uint));

  // If we're our own _mesh_sys, we'd better be in charge of
  // at least one coordinate, and we'd better have the same
  // FE type for all coordinates we are in charge of
  libmesh_assert_not_equal_to (mesh_xyz_var, libMesh::invalid_uint);
  libmesh_assert(!n_x_dofs || context.element_fe_var[_mesh_x_var] ==
                 context.element_fe_var[mesh_xyz_var]);
  libmesh_assert(!n_y_dofs || context.element_fe_var[_mesh_y_var] ==
                 context.element_fe_var[mesh_xyz_var]);
  libmesh_assert(!n_z_dofs || context.element_fe_var[_mesh_z_var] ==
                 context.element_fe_var[mesh_xyz_var]);

  const std::vector<std::vector<Real>>     & psi =
    context.element_fe_var[mesh_xyz_var]->get_phi();

  for (unsigned int var = 0; var != context.n_vars(); ++var)
    {
      // Mesh motion only affects time-evolving variables
      if (this->is_time_evolving(var))
        continue;

      // The mesh coordinate variables themselves are Lagrangian,
      // not Eulerian, and no convective term is desired.
      if (/*_mesh_sys == this && */
          (var == _mesh_x_var ||
           var == _mesh_y_var ||
           var == _mesh_z_var))
        continue;

      // Some of this code currently relies on the assumption that
      // we can pull mesh coordinate data from our own system
      if (_mesh_sys != this)
        libmesh_not_implemented();

      // This residual should only be called by unsteady solvers:
      // if the mesh is steady, there's no mesh convection term!
      UnsteadySolver * unsteady;
      if (this->time_solver->is_steady())
        return request_jacobian;
      else
        unsteady = cast_ptr<UnsteadySolver*>(this->time_solver.get());

      const std::vector<Real> & JxW =
        context.element_fe_var[var]->get_JxW();

      const std::vector<std::vector<Real>> & phi =
        context.element_fe_var[var]->get_phi();

      const std::vector<std::vector<RealGradient>> & dphi =
        context.element_fe_var[var]->get_dphi();

      const unsigned int n_u_dofs = context.dof_indices_var[var].size();

      DenseSubVector<Number> & Fu = *context.elem_subresiduals[var];
      DenseSubMatrix<Number> & Kuu = *context.elem_subjacobians[var][var];

      DenseSubMatrix<Number> * Kux = n_x_dofs ?
        context.elem_subjacobians[var][_mesh_x_var] : nullptr;
      DenseSubMatrix<Number> * Kuy = n_y_dofs ?
        context.elem_subjacobians[var][_mesh_y_var] : nullptr;
      DenseSubMatrix<Number> * Kuz = n_z_dofs ?
        context.elem_subjacobians[var][_mesh_z_var] : nullptr;

      std::vector<Real> delta_x(n_x_dofs, 0.);
      std::vector<Real> delta_y(n_y_dofs, 0.);
      std::vector<Real> delta_z(n_z_dofs, 0.);

      for (unsigned int i = 0; i != n_x_dofs; ++i)
        {
          unsigned int j = context.dof_indices_var[_mesh_x_var][i];
          delta_x[i] = libmesh_real(this->current_solution(j)) -
            libmesh_real(unsteady->old_nonlinear_solution(j));
        }

      for (unsigned int i = 0; i != n_y_dofs; ++i)
        {
          unsigned int j = context.dof_indices_var[_mesh_y_var][i];
          delta_y[i] = libmesh_real(this->current_solution(j)) -
            libmesh_real(unsteady->old_nonlinear_solution(j));
        }

      for (unsigned int i = 0; i != n_z_dofs; ++i)
        {
          unsigned int j = context.dof_indices_var[_mesh_z_var][i];
          delta_z[i] = libmesh_real(this->current_solution(j)) -
            libmesh_real(unsteady->old_nonlinear_solution(j));
        }

      for (unsigned int qp = 0; qp != n_qpoints; ++qp)
        {
          Gradient grad_u = context.interior_gradient(var, qp);
          RealGradient convection(0.);

          for (unsigned int i = 0; i != n_x_dofs; ++i)
            convection(0) += delta_x[i] * psi[i][qp];
          for (unsigned int i = 0; i != n_y_dofs; ++i)
            convection(1) += delta_y[i] * psi[i][qp];
          for (unsigned int i = 0; i != n_z_dofs; ++i)
            convection(2) += delta_z[i] * psi[i][qp];

          for (unsigned int i = 0; i != n_u_dofs; ++i)
            {
              Number JxWxPhiI = JxW[qp] * phi[i][qp];
              Fu(i) += (convection * grad_u) * JxWxPhiI;
              if (request_jacobian)
                {
                  Number JxWxPhiI = JxW[qp] * phi[i][qp];
                  for (unsigned int j = 0; j != n_u_dofs; ++j)
                    Kuu(i,j) += JxWxPhiI * (convection * dphi[j][qp]);

                  Number JxWxPhiIoverDT = JxWxPhiI/this->deltat;

                  Number JxWxPhiIxDUDXoverDT = JxWxPhiIoverDT * grad_u(0);
                  for (unsigned int j = 0; j != n_x_dofs; ++j)
                    (*Kux)(i,j) += JxWxPhiIxDUDXoverDT * psi[j][qp];

                  Number JxWxPhiIxDUDYoverDT = JxWxPhiIoverDT * grad_u(1);
                  for (unsigned int j = 0; j != n_y_dofs; ++j)
                    (*Kuy)(i,j) += JxWxPhiIxDUDYoverDT * psi[j][qp];

                  Number JxWxPhiIxDUDZoverDT = JxWxPhiIoverDT * grad_u(2);
                  for (unsigned int j = 0; j != n_z_dofs; ++j)
                    (*Kuz)(i,j) += JxWxPhiIxDUDZoverDT * psi[j][qp];
                }
            }
        }
    }
#endif // 0

  return request_jacobian;
}



bool FEMPhysics::mass_residual (bool request_jacobian,
                                DiffContext & c)
{
  FEMContext & context = cast_ref<FEMContext &>(c);

  unsigned int n_qpoints = context.get_element_qrule().n_points();

  for (unsigned int var = 0; var != context.n_vars(); ++var)
    {
      if (!this->is_time_evolving(var))
        continue;

      FEBase * elem_fe = nullptr;
      context.get_element_fe( var, elem_fe );

      const std::vector<Real> & JxW = elem_fe->get_JxW();

      const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();

      const unsigned int n_dofs = cast_int<unsigned int>
        (context.get_dof_indices(var).size());

      DenseSubVector<Number> & Fu = context.get_elem_residual(var);
      DenseSubMatrix<Number> & Kuu = context.get_elem_jacobian( var, var );

      for (unsigned int qp = 0; qp != n_qpoints; ++qp)
        {
          Number uprime;
          context.interior_rate(var, qp, uprime);
          const Number JxWxU = JxW[qp] * uprime;
          for (unsigned int i = 0; i != n_dofs; ++i)
            {
              Fu(i) -= JxWxU * phi[i][qp];
              if (request_jacobian && context.elem_solution_rate_derivative)
                {
                  const Number JxWxPhiIxDeriv = JxW[qp] * phi[i][qp] *
                    context.elem_solution_rate_derivative;
                  Kuu(i,i) -= JxWxPhiIxDeriv * phi[i][qp];
                  for (unsigned int j = i+1; j < n_dofs; ++j)
                    {
                      const Number Kij = JxWxPhiIxDeriv * phi[j][qp];
                      Kuu(i,j) -= Kij;
                      Kuu(j,i) -= Kij;
                    }
                }
            }
        }
    }

  return request_jacobian;
}

} // namespace libMesh
