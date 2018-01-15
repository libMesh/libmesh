// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Helper function for doing the projection
class MeshlessInterpolationFunction : public FunctionBase<Number>
{
public:
  MeshlessInterpolationFunction (const MeshfreeInterpolation & mfi,
                                 Threads::spin_mutex & mutex) :
    _mfi(mfi),
    _mutex(mutex)
  {}

  void init () {}
  void clear () {}

  virtual std::unique_ptr<FunctionBase<Number>> clone () const
  {
    return libmesh_make_unique<MeshlessInterpolationFunction>(_mfi, _mutex);
  }

  Number operator() (const Point & p,
                     const Real /*time*/)
  {
    _pts.clear();
    _pts.push_back(p);
    _vals.resize(1);

    Threads::spin_mutex::scoped_lock lock(_mutex);

    _mfi.interpolate_field_data(_mfi.field_variables(), _pts, _vals);

    return _vals.front();
  }


  void operator() (const Point & p,
                   const Real time,
                   DenseVector<Number> & output)
  {
    output.resize(1);
    output(0) = (*this)(p,time);
    return;
  }

private:
  const MeshfreeInterpolation & _mfi;
  mutable std::vector<Point> _pts;
  mutable std::vector<Number> _vals;
  Threads::spin_mutex & _mutex;
};

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

  // Create a MeshlessInterpolationFunction that uses our
  // InverseDistanceInterpolation object.  Since each
  // MeshlessInterpolationFunction shares the same
  // InverseDistanceInterpolation object in a threaded environment we
  // must also provide a locking mechanism.
  Threads::spin_mutex mutex;
  MeshlessInterpolationFunction mif(idi, mutex);

  // project the solution
  to_sys->project_solution(&mif);
}

} // namespace libMesh
