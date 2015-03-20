// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
#include "libmesh/cell_tet.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/face_tri3.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tet class member functions
dof_id_type Tet::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:
      return
        this->compute_key (this->node(0),
                           this->node(2),
                           this->node(1));

    case 1:
      return
        this->compute_key (this->node(0),
                           this->node(1),
                           this->node(3));

    case 2:
      return
        this->compute_key (this->node(1),
                           this->node(2),
                           this->node(3));

    case 3:
      return
        this->compute_key (this->node(2),
                           this->node(0),
                           this->node(3));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
}



UniquePtr<Elem> Tet::side (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_sides());

  Elem* face = new Tri3;

  switch (i)
    {
    case 0:
      {
        face->set_node(0) = this->get_node(0);
        face->set_node(1) = this->get_node(2);
        face->set_node(2) = this->get_node(1);
        break;
      }
    case 1:
      {
        face->set_node(0) = this->get_node(0);
        face->set_node(1) = this->get_node(1);
        face->set_node(2) = this->get_node(3);
        break;
      }
    case 2:
      {
        face->set_node(0) = this->get_node(1);
        face->set_node(1) = this->get_node(2);
        face->set_node(2) = this->get_node(3);
        break;
      }
    case 3:
      {
        face->set_node(0) = this->get_node(2);
        face->set_node(1) = this->get_node(0);
        face->set_node(2) = this->get_node(3);
        break;
      }
    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  return UniquePtr<Elem>(face);
}


void Tet::select_diagonal (const Diagonal diag) const
{
  libmesh_assert_equal_to (_diagonal_selection, INVALID_DIAG);
  _diagonal_selection = diag;
}





#ifdef LIBMESH_ENABLE_AMR


bool Tet::is_child_on_side_helper(const unsigned int c,
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

bool Tet::is_child_on_side_helper(const unsigned int /*c*/,
                                  const unsigned int /*s*/,
                                  const unsigned int /*checked_nodes*/[][3]) const
{
  libmesh_not_implemented();
  return false;
}

#endif //LIBMESH_ENABLE_AMR




void Tet::choose_diagonal() const
{
  // Check for uninitialized diagonal selection
  if (this->_diagonal_selection==INVALID_DIAG)
    {
      Real diag_01_23 = (this->point(0)+this->point(1)-this->point(2)-this->point(3)).size_sq();
      Real diag_02_13 = (this->point(0)-this->point(1)+this->point(2)-this->point(3)).size_sq();
      Real diag_03_12 = (this->point(0)-this->point(1)-this->point(2)+this->point(3)).size_sq();

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



bool Tet::is_edge_on_side(const unsigned int e,
                          const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (is_node_on_side(Tet4::edge_nodes_map[e][0],s) &&
          is_node_on_side(Tet4::edge_nodes_map[e][1],s));
}



Real Tet::quality(const ElemQuality q) const
{
  return Elem::quality(q); // Not implemented
}




std::pair<Real, Real> Tet::qual_bounds (const ElemQuality q) const
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
