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
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"

namespace libMesh
{

void FirstOrderUnsteadySolver::init()
{
  // First call parent init
  UnsteadySolver::init();

  // Now check for second order variables and add their velocities to the System.
  const std::set<unsigned int> & second_order_vars = this->get_second_order_vars();
  if( !second_order_vars.empty() )
    {
      for( std::set<unsigned int>::const_iterator var_it = second_order_vars.begin();
           var_it != second_order_vars.end(); ++var_it )
        {
          const Variable & var = _system.variable(*var_it);
          std::string new_var_name = std::string("dot_")+var.name();

          unsigned int v_var_idx;

          if( var.active_subdomains().empty() )
            v_var_idx = _system.add_variable( new_var_name, var.type() );
          else
            v_var_idx = _system.add_variable( new_var_name, var.type(), &var.active_subdomains() );

          _second_order_dot_vars.insert( std::pair<unsigned int,unsigned int>(*var_it,v_var_idx) );

          // The new velocities are time evolving variables of first order
          this->system().time_evolving( v_var_idx, 1 );

          // And if there are any boundary conditions set on the second order
          // variable, we also need to set it on its velocity variable.
          this->add_dot_var_dirichlet_bcs( *var_it, v_var_idx );
        }
    }
}

unsigned int FirstOrderUnsteadySolver::get_second_order_dot_var(unsigned int var) const
{
  libmesh_assert(this->is_second_order_var(var));
  libmesh_assert(_second_order_dot_vars.find(var) != _second_order_dot_vars.end());
  return _second_order_dot_vars.find(var)->second;
}

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
      if (!this->is_second_order_var(var))
        continue;

      unsigned int dot_var = this->get_second_order_dot_var(var);

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


void FirstOrderUnsteadySolver::add_dot_var_dirichlet_bcs( unsigned int var_idx,
                                                          unsigned int dot_var_idx )
{
  // We're assuming that there could be a lot more variables than
  // boundary conditions, so we search each of the boundary conditions
  // for this variable rather than looping over boundary conditions
  // in a separate loop and searching through all the variables.
  const DirichletBoundaries * all_dbcs =
    this->system().get_dof_map().get_dirichlet_boundaries();

  if( all_dbcs )
    {
      // We need to cache the DBCs to be added so that we add them
      // after looping over the existing DBCs. Otherwise, we're polluting
      // the thing we're looping over.
      std::vector<DirichletBoundary*> new_dbcs;

      DirichletBoundaries::const_iterator dbc_it = all_dbcs->begin();
      for( ; dbc_it != all_dbcs->end(); ++dbc_it )
        {
          libmesh_assert(*dbc_it);
          DirichletBoundary & dbc = *(*dbc_it);

          // Look for second order variable in the current
          // DirichletBoundary object
          std::vector<unsigned int>::const_iterator dbc_var_it =
            std::find( dbc.variables.begin(), dbc.variables.end(), var_idx );

          // If we found it, then we also need to add it's corresponding
          // "dot" variable to a DirichletBoundary
          std::vector<unsigned int> vars_to_add;
          if( dbc_var_it != dbc.variables.end() )
            vars_to_add.push_back(dot_var_idx);

          if( !vars_to_add.empty() )
            {
              DirichletBoundary * new_dbc;

              if( dbc.f )
                {
                  if( dbc.g )
                    new_dbc = new DirichletBoundary(dbc.b,
                                                    vars_to_add,
                                                    dbc.f.get(),
                                                    dbc.g.get());

                  else
                    new_dbc = new DirichletBoundary(dbc.b,
                                                    vars_to_add,
                                                    dbc.f.get());
                }
              else if( dbc.f_fem )
                {
                  if( dbc.g_fem )
                    new_dbc = new DirichletBoundary(dbc.b,
                                                    vars_to_add,
                                                    this->system(),
                                                    dbc.f_fem.get(),
                                                    dbc.g_fem.get());

                  else
                    new_dbc = new DirichletBoundary(dbc.b,
                                                    vars_to_add,
                                                    this->system(),
                                                    dbc.f_fem.get());
                }
              else
                libmesh_error_msg("Could not find valid boundary function!");

              new_dbcs.push_back(new_dbc);
            }
        }

      // Now add the new DBCs for the "dot" vars to the DofMap
      std::vector<DirichletBoundary*>::iterator new_dbc_it =
        new_dbcs.begin();

      for( ; new_dbc_it != new_dbcs.end(); ++new_dbc_it )
        {
          const DirichletBoundary & dbc = *(*new_dbc_it);
          this->system().get_dof_map().add_dirichlet_boundary(dbc);
          delete *new_dbc_it;
        }

    } // if(all_dbcs)
}

} // namespace libMesh
