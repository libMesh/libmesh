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



#ifndef LIBMESH_LOCATION_MAPS_H
#define LIBMESH_LOCATION_MAPS_H

#include "libmesh/libmesh_config.h"

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"

// C++ Includes
#include <unordered_map>
#include <vector>

namespace libMesh
{

// Forward Declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;
template <typename> class NodeTempl;
typedef NodeTempl<Real> Node;


/**
 * Data structures that enable location-based lookups
 * The key is a hash of the Point location.
 * For efficiency we will use a hashed multimap if it is
 * available, otherwise a regular multimap.
 *
 * \author Roy Stogner
 * \date 2008
 * \brief std::map-like data structure using hashed Points for keys.
 */
template <typename T>
class LocationMap
{
  typedef std::unordered_multimap<unsigned int, T *> map_type;
public:
  void init(MeshBase &);

  void clear() { _map.clear(); }

  void insert(T &);

  bool empty() const { return _map.empty(); }

  T * find(const Point &,
           const Real tol = TOLERANCE);

  Point point_of(const T &) const;

protected:
  unsigned int key(const Point &);

  void fill(MeshBase &);

private:
  map_type          _map;
  std::vector<Real> _lower_bound;
  std::vector<Real> _upper_bound;
};

} // namespace libMesh


#endif // LIBMESH_LOCATION_MAPS_H
