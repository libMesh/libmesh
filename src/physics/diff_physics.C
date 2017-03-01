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


#include "libmesh/diff_context.h"
#include "libmesh/diff_physics.h"
#include "libmesh/system.h"

namespace libMesh
{

DifferentiablePhysics::~DifferentiablePhysics()
{
  DifferentiablePhysics::clear_physics();
}



void DifferentiablePhysics::clear_physics ()
{
  _time_evolving.resize(0);
}



void DifferentiablePhysics::init_physics (const System & sys)
{
  // give us flags for every variable that might be time evolving
  _time_evolving.resize(sys.n_vars(), false);
}

void DifferentiablePhysics::time_evolving (unsigned int var,
                                           unsigned int order)
{
  if (order != 1 && order != 2)
    libmesh_error_msg("Input order must be 1 or 2!");

  if (_time_evolving.size() <= var)
    _time_evolving.resize(var+1, 0);

  _time_evolving[var] = order;

  if( order == 1 )
    _first_order_vars.insert(var);
  else
    _second_order_vars.insert(var);
}

bool DifferentiablePhysics::nonlocal_mass_residual(bool request_jacobian,
                                                   DiffContext & c)
{
  FEMContext & context = cast_ref<FEMContext &>(c);

  for (unsigned int var = 0; var != context.n_vars(); ++var)
    {
      if (!this->is_time_evolving(var))
        continue;

      if (c.get_system().variable(var).type().family != SCALAR)
        continue;

      const std::vector<dof_id_type> & dof_indices =
        context.get_dof_indices(var);

      const unsigned int n_dofs = cast_int<unsigned int>
        (dof_indices.size());

      DenseSubVector<Number> & Fs = context.get_elem_residual(var);
      DenseSubMatrix<Number> & Kss = context.get_elem_jacobian( var, var );

      const libMesh::DenseSubVector<libMesh::Number> & Us =
        context.get_elem_solution(var);

      for (unsigned int i=0; i != n_dofs; ++i)
        {
          Fs(i) -= Us(i);

          if (request_jacobian)
            Kss(i,i) -= context.elem_solution_rate_derivative;
        }
    }

  return request_jacobian;
}



bool DifferentiablePhysics::_eulerian_time_deriv (bool request_jacobian,
                                                  DiffContext & context)
{
  // For any problem we need time derivative terms
  request_jacobian =
    this->element_time_derivative(request_jacobian, context);

  // For a moving mesh problem we may need the pseudoconvection term too
  return this->eulerian_residual(request_jacobian, context) &&
    request_jacobian;
}

} // namespace libMesh
