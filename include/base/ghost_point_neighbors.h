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



#ifndef LIBMESH_GHOST_POINT_NEIGHBORS_H
#define LIBMESH_GHOST_POINT_NEIGHBORS_H

// Local Includes
#include "libmesh/ghosting_functor.h"

namespace libMesh
{
template <typename> class ElemTempl;
template <typename> class RemoteElemTempl;
template <typename> class MeshBaseTempl;
template <typename> class NodeTempl;

/**
 * This class implements the default geometry ghosting in libMesh:
 * point neighbors and interior_parent elements are ghosted.
 *
 * \author Roy H. Stogner
 * \date 2016
 */
template <typename RealType = Real>
class GhostPointNeighborsTempl : public GhostingFunctorTempl<RealType>
{
public:
  typedef GhostPointNeighborsTempl<RealType> GhostPointNeighbors;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef ElemTempl<RealType> Elem;
  typedef RemoteElemTempl<RealType> RemoteElem;
  typedef NodeTempl<RealType> Node;
  using map_type = std::unordered_map<const ElemTempl<RealType> *, const CouplingMatrix*>;

  /**
   * Constructor.
   */
  GhostPointNeighborsTempl(const MeshBase & mesh) : _mesh(mesh) {}

  /**
   * For the specified range of active elements, find their point
   * neighbors and interior_parent elements, ignoring those on
   * processor p.
   */
  virtual void operator() (const typename MeshBase::const_element_iterator & range_begin,
                           const typename MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements);

private:

  const MeshBase & _mesh;
};

typedef GhostPointNeighborsTempl<Real> GhostPointNeighbors;

} // namespace libMesh

#endif // LIBMESH_GHOST_POINT_NEIGHBORS_H
