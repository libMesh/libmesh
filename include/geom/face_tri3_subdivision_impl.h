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

#ifndef LIBMESH_FACE_TRI3_SUBDIVISION_IMPL_H
#define LIBMESH_FACE_TRI3_SUBDIVISION_IMPL_H

// Local includes
#include "libmesh/face_tri3_subdivision.h"
#include "libmesh/mesh_subdivision_support.h"
#include "libmesh/enum_order.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tri3 subdivision class member functions
template <typename RealType>
Tri3SubdivisionTempl<RealType>::Tri3SubdivisionTempl(Elem * p) : Tri3(p), _subdivision_updated(true)
{
  if (p)
    {
      libmesh_assert_equal_to(p->type(), TRI3SUBDIVISION);
      Tri3Subdivision * sd_elem = static_cast<Tri3Subdivision *>(p);
      _is_ghost = sd_elem->is_ghost();

      if (!_is_ghost)
        {
          _ordered_nodes[0] = sd_elem->get_ordered_node(0);
          _ordered_nodes[1] = sd_elem->get_ordered_node(1);
          _ordered_nodes[2] = sd_elem->get_ordered_node(2);
        }
    }
}



template <typename RealType>
Order Tri3SubdivisionTempl<RealType>::default_order() const
{
  return FOURTH;
}



template <typename RealType>
void Tri3SubdivisionTempl<RealType>::prepare_subdivision_properties()
{
  /*
   * Find the index of the irregular vertex, if any.
   * The current implementation can only handle triangles with
   * no more than one irregular vertex. That is, a vertex with
   * valence != 6.
   */
  unsigned int irregular_idx = 0;
  for (unsigned int i = 0; i < 3; ++i)
    {
      if (this->node_ptr(i)->valence() != 6)
        {
          irregular_idx = i;
          if (this->node_ptr(MeshTools::Subdivision::next[i])->valence() != 6 || this->node_ptr(MeshTools::Subdivision::prev[i])->valence() != 6)
            libmesh_error_msg("Error: The mesh contains elements with more than one irregular vertex!");
        }
    }

  /*
   * Rotate ordered vertices such that ordered_nodes[0] is the
   * irregular vertex. Doing this once in advance lets the evaluation
   * of subdivision interpolation be much more efficient afterward.
   */
  switch (irregular_idx)
    {
    case 0:
      _ordered_nodes[0] = this->node_ptr(0);
      _ordered_nodes[1] = this->node_ptr(1);
      _ordered_nodes[2] = this->node_ptr(2);
      break;
    case 1:
      _ordered_nodes[0] = this->node_ptr(1);
      _ordered_nodes[1] = this->node_ptr(2);
      _ordered_nodes[2] = this->node_ptr(0);
      break;
    case 2:
      _ordered_nodes[0] = this->node_ptr(2);
      _ordered_nodes[1] = this->node_ptr(0);
      _ordered_nodes[2] = this->node_ptr(1);
      break;
    default:
      libmesh_error_msg("Unrecognized irregular_idx = " << irregular_idx);
    }

  _subdivision_updated = true;
}


template <typename RealType>
unsigned int Tri3SubdivisionTempl<RealType>::local_node_number(unsigned int node_id) const
{
  return (this->_nodes[0]->id() == node_id) ? 0 : ( (this->_nodes[1]->id() == node_id) ? 1 : ( (this->_nodes[2]->id() == node_id) ? 2 : 3 ) );
}


template <typename RealType>
unsigned int Tri3SubdivisionTempl<RealType>::get_ordered_valence(unsigned int node_id) const
{
  libmesh_assert_less(node_id, this->n_neighbors());
  libmesh_assert(_subdivision_updated);
  return get_ordered_node(node_id)->valence();
}


template <typename RealType>
NodeTempl<RealType> * Tri3SubdivisionTempl<RealType>::get_ordered_node(unsigned int node_id) const
{
  libmesh_assert_less(node_id, 3);
  libmesh_assert(_subdivision_updated);
  return _ordered_nodes[node_id];
}

} // namespace libMesh

#endif // LIBMESH_FACE_TRI3_SUBDIVISION_IMPL_H
