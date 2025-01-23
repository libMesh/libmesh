// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/meshfree_solution_transfer.h"

#include "libmesh/mesh.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/threads.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/function_base.h"
#include "libmesh/node.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward Declarations
template <typename T>
class DenseVector;

void
MeshfreeSolutionTransfer::transfer(const Variable & from_var,
                                   const Variable & to_var)
{
  libmesh_experimental();

  System * from_sys = from_var.system();
  System * to_sys = to_var.system();

  EquationSystems & from_es = from_sys->get_equation_systems();

  MeshBase & from_mesh = from_es.get_mesh();

  InverseDistanceInterpolation<LIBMESH_DIM> idi
    (from_mesh.comm(), 4, 2);

  std::vector<Point>  & src_pts  (idi.get_source_points());
  std::vector<Number> & src_vals (idi.get_source_vals());

  std::vector<std::string> field_vars;
  field_vars.push_back(from_var.name());
  idi.set_field_variables(field_vars);

  // We now will loop over every node in the source mesh
  // and add it to a source point list, along with the solution
  for (const auto & node : from_mesh.local_node_ptr_range())
    {
      src_pts.push_back(*node);
      src_vals.push_back((*from_sys->solution)(node->dof_number(from_sys->number(),from_var.number(),0)));
    }

  // We have only set local values - prepare for use by gathering remote data
  idi.prepare_for_use();

  // Create a MeshfreeInterpolationFunction that uses our
  // InverseDistanceInterpolation object.  Since each
  // MeshfreeInterpolationFunction shares the same
  // InverseDistanceInterpolation object in a threaded environment we
  // must also provide a locking mechanism.
  Threads::spin_mutex mutex;
  MeshfreeInterpolationFunction mif(idi, mutex);

  // project the solution
  to_sys->project_solution(&mif);
}

} // namespace libMesh
