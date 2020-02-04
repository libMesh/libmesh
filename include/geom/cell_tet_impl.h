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

#ifndef LIBMESH_CELL_TET_IMPL_H
#define LIBMESH_CELL_TET_IMPL_H

// Local includes
#include "libmesh/cell_tet.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/enum_elem_quality.h"

namespace libMesh
{

// ------------------------------------------------------------
// Tet class member functions
template <typename RealType>
dof_id_type TetTempl<RealType>::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  return this->compute_key(this->node_id(Tet4::side_nodes_map[s][0]),
                           this->node_id(Tet4::side_nodes_map[s][1]),
                           this->node_id(Tet4::side_nodes_map[s][2]));
}



template <typename RealType>
unsigned int TetTempl<RealType>::which_node_am_i(unsigned int side,
                                  unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, 3);

  return Tet4::side_nodes_map[side][side_node];
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> TetTempl<RealType>::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face = libmesh_make_unique<Tri3>();

  for (auto n : face->node_index_range())
    face->set_node(n) = this->node_ptr(Tet4::side_nodes_map[i][n]);

  return face;
}



template <typename RealType>
void TetTempl<RealType>::side_ptr (std::unique_ptr<Elem> & side,
                    const unsigned int i)
{
  this->template simple_side_ptr<Tet,Tet4>(side, i, TRI3);
}



template <typename RealType>
void TetTempl<RealType>::select_diagonal (const Diagonal diag) const
{
  libmesh_assert_equal_to (_diagonal_selection, INVALID_DIAG);
  _diagonal_selection = diag;
}





#ifdef LIBMESH_ENABLE_AMR


template <typename RealType>
bool TetTempl<RealType>::is_child_on_side_helper(const unsigned int c,
                                  const unsigned int s,
                                  const unsigned int checked_nodes[][3]) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  // For the 4 vertices, child c touches vertex c, so we can return
  // true if that vertex is on side s
  for (unsigned int i = 0; i != 3; ++i)
    if (Tet4::side_nodes_map[s][i] == c)
      return true;

  // If we are a "vertex child" and we didn't already return true,
  // we must not be on the side in question
  if (c < 4)
    return false;

  // For the 4 non-vertex children, the child ordering depends on the
  // diagonal selection.  We'll let the embedding matrix figure that
  // out: if this child has three nodes that don't depend on the
  // position of the node_facing_side[s], then we're on side s.  Which
  // three nodes those are depends on the subclass, so their responsibility
  // is to call this function with the proper check_nodes array
  const unsigned int node_facing_side[4] = {3, 2, 0, 1};
  const unsigned int n = node_facing_side[s];

  // Add up the absolute values of the entries of the embedding matrix for the
  // nodes opposite node n.  If it is equal to zero, then the child in question is
  // on side s, so return true.
  Real embedding_sum = 0.;
  for (unsigned i=0; i<3; ++i)
    embedding_sum += std::abs(this->embedding_matrix(c, checked_nodes[n][i], n));

  return ( std::abs(embedding_sum) < 1.e-3 );
}

#else

template <typename RealType>
bool TetTempl<RealType>::is_child_on_side_helper(const unsigned int /*c*/,
                                  const unsigned int /*s*/,
                                  const unsigned int /*checked_nodes*/[][3]) const
{
  libmesh_not_implemented();
  return false;
}

#endif //LIBMESH_ENABLE_AMR




template <typename RealType>
void TetTempl<RealType>::choose_diagonal() const
{
  // Check for uninitialized diagonal selection
  if (this->_diagonal_selection==INVALID_DIAG)
    {
      Real diag_01_23 = (this->point(0) + this->point(1) - this->point(2) - this->point(3)).norm_sq();
      Real diag_02_13 = (this->point(0) - this->point(1) + this->point(2) - this->point(3)).norm_sq();
      Real diag_03_12 = (this->point(0) - this->point(1) - this->point(2) + this->point(3)).norm_sq();

      this->_diagonal_selection=DIAG_02_13;

      if (diag_01_23 < diag_02_13 || diag_03_12 < diag_02_13)
        {
          if (diag_01_23 < diag_03_12)
            this->_diagonal_selection=DIAG_01_23;

          else
            this->_diagonal_selection=DIAG_03_12;
        }
    }
}



template <typename RealType>
bool TetTempl<RealType>::is_edge_on_side(const unsigned int e,
                          const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (this->is_node_on_side(Tet4::edge_nodes_map[e][0],s) &&
          this->is_node_on_side(Tet4::edge_nodes_map[e][1],s));
}



template <typename RealType>
Real TetTempl<RealType>::quality(const ElemQuality q) const
{
  return Elem::quality(q); // Not implemented
}




template <typename RealType>
std::pair<Real, Real> TetTempl<RealType>::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {

    case ASPECT_RATIO_BETA:
    case ASPECT_RATIO_GAMMA:
      bounds.first  = 1.;
      bounds.second = 3.;
      break;

    case SIZE:
    case SHAPE:
      bounds.first  = 0.2;
      bounds.second = 1.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 3.;
      break;

    case DISTORTION:
      bounds.first  = 0.6;
      bounds.second = 1.;
      break;

    case JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.414;
      break;

    default:
      libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}

} // namespace libMesh

#endif // LIBMESH_CELL_TET_IMPL_H
