// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/location_maps.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/node.h"
#include "libmesh/parallel.h"

// C++ Includes
#include <limits>
#include <utility>

namespace
{
using libMesh::Real;

// 10 bits per coordinate, to work with 32+ bit machines
const unsigned int chunkmax = 1024;
const Real chunkfloat = 1024.0;
}



namespace libMesh
{

template <typename T>
void LocationMap<T>::init(MeshBase & mesh)
{
  // This function must be run on all processors at once
  // for non-serial meshes
  if (!mesh.is_serial())
    libmesh_parallel_only(mesh.comm());

  LOG_SCOPE("init()", "LocationMap");

  // Clear the old map
  _map.clear();

  // Cache a bounding box
  _bbox = MeshTools::create_nodal_bounding_box(mesh);

  this->fill(mesh);
}



template <typename T>
void LocationMap<T>::insert(T & t)
{
  this->_map.emplace(this->key(this->point_of(t)), &t);
}



template <>
Point LocationMap<Node>::point_of(const Node & node) const
{
  return node;
}



template <>
Point LocationMap<Elem>::point_of(const Elem & elem) const
{
  return elem.vertex_average();
}



template <typename T>
T * LocationMap<T>::find(const Point & p,
                         const Real tol)
{
  LOG_SCOPE("find()", "LocationMap");

  // Look for a likely key in the multimap
  unsigned int pointkey = this->key(p);

  // Look for the exact key first
  for (const auto & pr : as_range(_map.equal_range(pointkey)))
    if (p.absolute_fuzzy_equals(this->point_of(*(pr.second)), tol))
      return pr.second;

  // Look for neighboring bins' keys next
  for (int xoffset = -1; xoffset != 2; ++xoffset)
    {
      for (int yoffset = -1; yoffset != 2; ++yoffset)
        {
          for (int zoffset = -1; zoffset != 2; ++zoffset)
            {
              auto key_pos = _map.equal_range(pointkey +
                                              xoffset*chunkmax*chunkmax +
                                              yoffset*chunkmax +
                                              zoffset);
              for (const auto & pr : as_range(key_pos))
                if (p.absolute_fuzzy_equals(this->point_of(*(pr.second)), tol))
                  return pr.second;
            }
        }
    }

  return nullptr;
}



template <typename T>
unsigned int LocationMap<T>::key(const Point & p)
{
  Real xscaled = 0., yscaled = 0., zscaled = 0.;

  const Real deltax = _bbox.max()(0) - _bbox.min()(0);

  if (std::abs(deltax) > TOLERANCE)
    xscaled = (p(0) - _bbox.min()(0))/deltax;

  // Only check y-coords if libmesh is compiled with LIBMESH_DIM>1
#if LIBMESH_DIM > 1
  const Real deltay = _bbox.max()(1) - _bbox.min()(1);

  if (std::abs(deltay) > TOLERANCE)
    yscaled = (p(1) - _bbox.min()(1))/deltay;
#endif

  // Only check z-coords if libmesh is compiled with LIBMESH_DIM>2
#if LIBMESH_DIM > 2
  const Real deltaz = _bbox.max()(2) - _bbox.min()(2);

  if (std::abs(deltaz) > TOLERANCE)
    zscaled = (p(2) - _bbox.min()(2))/deltaz;
#endif

  unsigned int n0 = static_cast<unsigned int> (chunkfloat * xscaled),
    n1 = static_cast<unsigned int> (chunkfloat * yscaled),
    n2 = static_cast<unsigned int> (chunkfloat * zscaled);

  return chunkmax*chunkmax*n0 + chunkmax*n1 + n2;
}



template <>
void LocationMap<Node>::fill(MeshBase & mesh)
{
  // Populate the nodes map
  for (auto & node : mesh.node_ptr_range())
    this->insert(*node);
}



template <>
void LocationMap<Elem>::fill(MeshBase & mesh)
{
  // Populate the elem map
  for (auto & elem : mesh.active_element_ptr_range())
    this->insert(*elem);
}



template class LIBMESH_EXPORT LocationMap<Elem>;
template class LIBMESH_EXPORT LocationMap<Node>;

} // namespace libMesh
