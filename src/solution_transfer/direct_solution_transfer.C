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



// C++ includes
#include "libmesh/direct_solution_transfer.h"

#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"

namespace libMesh {

DirectSolutionTransfer::DirectSolutionTransfer(const libMesh::Parallel::Communicator & comm_in) :
  SolutionTransfer(comm_in)
{}

DirectSolutionTransfer::~DirectSolutionTransfer()
{}

void
DirectSolutionTransfer::transfer(const Variable & from_var,
                                 const Variable & to_var)
{
  libmesh_experimental();

  System * from_sys = from_var.system();
  System * to_sys = to_var.system();

  // Just a couple of (not completely thorough)
  libmesh_assert(from_sys->get_equation_systems().get_mesh().n_nodes() == from_sys->get_equation_systems().get_mesh().n_nodes());
  libmesh_assert(from_var.type() == to_var.type());

  // get dof indices for source variable
  unsigned int from_vn = from_var.number();
  std::set<dof_id_type> from_var_indices;
  from_sys->local_dof_indices(from_vn, from_var_indices);

  // get dof indices for dest variable
  unsigned int to_vn = to_var.number();
  std::set<dof_id_type> to_var_indices;
  to_sys->local_dof_indices(to_vn, to_var_indices);

  // copy the values from from solution vector to to solution vector
  std::set<dof_id_type>::iterator from_it = from_var_indices.begin();
  std::set<dof_id_type>::iterator from_it_end = from_var_indices.end();
  std::set<dof_id_type>::iterator to_it = to_var_indices.begin();

  NumericVector<Number> & from_solution = *from_sys->solution;

  for (; from_it != from_it_end; ++from_it, ++to_it)
    to_sys->solution->set(*to_it, from_solution(*from_it));

  to_sys->solution->close();
  to_sys->update();
}

} // namespace libMesh
