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

// Local Includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/periodic_boundary.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{

PeriodicBoundary::PeriodicBoundary() :
  PeriodicBoundaryBase(),
  translation_vector()
{
}



PeriodicBoundary::PeriodicBoundary(const PeriodicBoundary & o, TransformationType t) :
  PeriodicBoundaryBase(o),
  translation_vector(o.translation_vector)
{
  if (t == INVERSE)
    {
      std::swap(myboundary, pairedboundary);
      translation_vector *= -1.0;
    }
}



PeriodicBoundary::PeriodicBoundary(const RealVectorValue & vector) :
  PeriodicBoundaryBase(),
  translation_vector(vector)
{
}



Point PeriodicBoundary::get_corresponding_pos(const Point & pt) const
{
  return pt + translation_vector;
}



std::unique_ptr<PeriodicBoundaryBase> PeriodicBoundary::clone(TransformationType t) const
{
  return libmesh_make_unique<PeriodicBoundary>(*this, t);
}


} // namespace libMesh


#endif // LIBMESH_ENABLE_PERIODIC
