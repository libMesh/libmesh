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

#include "libmesh/second_order_unsteady_solver.h"

#include "libmesh/diff_system.h"
#include "libmesh/dof_map.h"

namespace libMesh
{
SecondOrderUnsteadySolver::SecondOrderUnsteadySolver (sys_type & s)
  : UnsteadySolver(s),
    _old_local_solution_rate(NumericVector<Number>::build(s.comm())),
    _old_local_solution_accel(NumericVector<Number>::build(s.comm()))
{}

SecondOrderUnsteadySolver::~SecondOrderUnsteadySolver ()
{}

void SecondOrderUnsteadySolver::init ()
{
  UnsteadySolver::init();

  _system.add_vector("_old_solution_rate");
  _system.add_vector("_old_solution_accel");
}

void SecondOrderUnsteadySolver::init_data()
{
  UnsteadySolver::init_data();

#ifdef LIBMESH_ENABLE_GHOSTED
  _old_local_solution_rate->init (_system.n_dofs(), _system.n_local_dofs(),
                                  _system.get_dof_map().get_send_list(), false,
                                  GHOSTED);

  _old_local_solution_accel->init (_system.n_dofs(), _system.n_local_dofs(),
                                   _system.get_dof_map().get_send_list(), false,
                                   GHOSTED);
#else
  _old_local_solution_rate->init (_system.n_dofs(), false, SERIAL);
  _old_local_solution_accel->init (_system.n_dofs(), false, SERIAL);
#endif
}

void SecondOrderUnsteadySolver::reinit ()
{
  UnsteadySolver::reinit();

#ifdef LIBMESH_ENABLE_GHOSTED
  _old_local_solution_rate->init (_system.n_dofs(), _system.n_local_dofs(),
                                  _system.get_dof_map().get_send_list(), false,
                                  GHOSTED);

  _old_local_solution_accel->init (_system.n_dofs(), _system.n_local_dofs(),
                                   _system.get_dof_map().get_send_list(), false,
                                   GHOSTED);
#else
  _old_local_solution_rate->init (_system.n_dofs(), false, SERIAL);
  _old_local_solution_accel->init (_system.n_dofs(), false, SERIAL);
#endif

  // localize the old solutions
  NumericVector<Number> & old_solution_rate =
    _system.get_vector("_old_solution_rate");

  NumericVector<Number> & old_solution_accel =
    _system.get_vector("_old_solution_accel");

  old_solution_rate.localize
    (*_old_local_solution_rate,
     _system.get_dof_map().get_send_list());

  old_solution_accel.localize
    (*_old_local_solution_accel,
     _system.get_dof_map().get_send_list());
}

void SecondOrderUnsteadySolver::retrieve_timestep()
{
  libmesh_not_implemented();
}

void SecondOrderUnsteadySolver::project_initial_rate(FunctionBase<Number> * f,
                                                     FunctionBase<Gradient> * g)
{
  NumericVector<Number> & old_solution_rate =
    _system.get_vector("_old_solution_rate");

  _system.project_vector( old_solution_rate, f, g );
}

Number SecondOrderUnsteadySolver::old_solution_rate(const dof_id_type global_dof_number)
  const
{
  libmesh_assert_less (global_dof_number, _system.get_dof_map().n_dofs());
  libmesh_assert_less (global_dof_number, _old_local_solution_rate->size());

  return (*_old_local_solution_rate)(global_dof_number);
}

Number SecondOrderUnsteadySolver::old_solution_accel(const dof_id_type global_dof_number)
  const
{
  libmesh_assert_less (global_dof_number, _system.get_dof_map().n_dofs());
  libmesh_assert_less (global_dof_number, _old_local_solution_accel->size());

  return (*_old_local_solution_accel)(global_dof_number);
}

} // end namespace libMesh
