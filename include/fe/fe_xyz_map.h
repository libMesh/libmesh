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



#ifndef LIBMESH_FE_XYZ_MAP_H
#define LIBMESH_FE_XYZ_MAP_H

#include "libmesh/fe_map.h"
#include "libmesh/libmesh_logging.h"

namespace libMesh
{

// Forward declarations
class Elem;

class FEXYZMap : public FEMap
{
public:

  FEXYZMap()
    : FEMap()
  {
    // All FEXYZ objects are going to be querying xyz coordinates
    calculate_xyz = true;
  }

  virtual ~FEXYZMap(){}

  /**
   * Special implementation for XYZ finite elements
   */
  virtual void compute_face_map(int dim,
                                const std::vector<Real> & qw,
                                const Elem * side) libmesh_override;

}; // class FEXYZMap
} // namespace libMesh

#endif // LIBMESH_FE_XYZ_MAP_H
