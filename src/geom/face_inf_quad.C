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


// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes cont'd
#include "libmesh/face_inf_quad.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_inf_edge2.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfQuad class member functions
dof_id_type InfQuad::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());


  switch (s)
    {
    case 0:

      return
        this->compute_key (this->node(0),
                           this->node(1));

    case 1:

      return
        this->compute_key (this->node(1),
                           this->node(3));

    case 2:

      return
        this->compute_key (this->node(0),
                           this->node(2));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
}



UniquePtr<Elem> InfQuad::side (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_sides());

  // To be returned wrapped in an UniquePtr
  Elem* edge = NULL;

  switch (i)
    {
    case 0:
      {
        // base face
        edge = new Edge2;
        edge->set_node(0) = this->get_node(0);
        edge->set_node(1) = this->get_node(1);
        break;
      }

    case 1:
      {
        // adjacent to another infinite element
        edge = new InfEdge2;
        edge->set_node(0) = this->get_node(1);
        edge->set_node(1) = this->get_node(3);
        break;
      }

    case 2:
      {
        // adjacent to another infinite element
        edge = new InfEdge2;
        edge->set_node(0) = this->get_node(0); // be aware of swapped nodes,
        edge->set_node(1) = this->get_node(2); // compared to conventional side numbering
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  return UniquePtr<Elem>(edge);
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
