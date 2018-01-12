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


// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes cont'd
#include "libmesh/face_inf_quad.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/face_inf_quad4.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfQuad class static member initializations


// We need to require C++11...
const Real InfQuad::_master_points[6][3] =
  {
    {-1, 0},
    {1, 0},
    {-1, 1},
    {1, 1},
    {0, 0},
    {0, 1}
  };


// ------------------------------------------------------------
// InfQuad class member functions
dof_id_type InfQuad::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  // The order of the node ids does not matter, they are sorted by the
  // compute_key() function.
  return this->compute_key(this->node_id(InfQuad4::side_nodes_map[s][0]),
                           this->node_id(InfQuad4::side_nodes_map[s][1]));
}



unsigned int InfQuad::which_node_am_i(unsigned int side,
                                      unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, 2);

  return InfQuad4::side_nodes_map[side][side_node];
}



std::unique_ptr<Elem> InfQuad::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  // Return value
  std::unique_ptr<Elem> edge;

  switch (i)
    {
    case 0: // base face
      {
        edge = libmesh_make_unique<Edge2>();
        break;
      }

    case 1: // adjacent to another infinite element
    case 2: // adjacent to another infinite element
      {
        edge = libmesh_make_unique<InfEdge2>();
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (unsigned n=0; n<edge->n_nodes(); ++n)
    edge->set_node(n) = this->node_ptr(InfQuad4::side_nodes_map[i][n]);

  return edge;
}



bool InfQuad::is_child_on_side(const unsigned int c,
                               const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  return (s == 0 || s == c+1);
}



Real InfQuad::quality (const ElemQuality) const
{
  return 0.; // Not implemented
}




std::pair<Real, Real> InfQuad::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {

    case ASPECT_RATIO:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case SKEW:
      bounds.first  = 0.;
      bounds.second = 0.5;
      break;

    case TAPER:
      bounds.first  = 0.;
      bounds.second = 0.7;
      break;

    case WARP:
      bounds.first  = 0.9;
      bounds.second = 1.;
      break;

    case STRETCH:
      bounds.first  = 0.25;
      bounds.second = 1.;
      break;

    case MIN_ANGLE:
      bounds.first  = 45.;
      bounds.second = 90.;
      break;

    case MAX_ANGLE:
      bounds.first  = 90.;
      bounds.second = 135.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.;
      break;

    case SHEAR:
    case SHAPE:
    case SIZE:
      bounds.first  = 0.3;
      bounds.second = 1.;
      break;

    case DISTORTION:
      bounds.first  = 0.6;
      bounds.second = 1.;
      break;

    default:
      libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}

} // namespace libMesh



#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
