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



#include "libmesh/meshfunction_solution_transfer.h"

#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_function.h"
#include "libmesh/node.h"

namespace libMesh
{

MeshFunctionSolutionTransfer::MeshFunctionSolutionTransfer(const libMesh::Parallel::Communicator & comm_in) :
  SolutionTransfer(comm_in)
{}

MeshFunctionSolutionTransfer::~MeshFunctionSolutionTransfer()
{}

void
MeshFunctionSolutionTransfer::transfer(const Variable & from_var,
                                       const Variable & to_var)
{
  // This only works when transferring to a Lagrange variable
  libmesh_assert(to_var.type().family == LAGRANGE);

  unsigned int to_var_num = to_var.number();

  System * from_sys = from_var.system();
  System * to_sys = to_var.system();

  // Only works with a serialized mesh to transfer from!
  libmesh_assert(from_sys->get_mesh().is_serial());

  unsigned int to_sys_num = to_sys->number();

  EquationSystems & from_es = from_sys->get_equation_systems();

  //Create a serialized version of the solution vector
  std::unique_ptr<NumericVector<Number>> serialized_solution = NumericVector<Number>::build(from_sys->get_mesh().comm());
  serialized_solution->init(from_sys->n_dofs(), false, SERIAL);

  // Need to pull down a full copy of this vector on every processor
  // so we can get values in parallel
  from_sys->solution->localize(*serialized_solution);

  MeshFunction from_func(from_es, *serialized_solution, from_sys->get_dof_map(), to_var_num);
  from_func.init();

  // Now loop over the nodes of the 'To' mesh setting values for each variable.
  for (const auto & node : to_sys->get_mesh().local_node_ptr_range())
    to_sys->solution->set(node->dof_number(to_sys_num, to_var_num, 0), from_func(*node)); // 0 is for the value component

  to_sys->solution->close();
  to_sys->update();
}

} // namespace libMesh
