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
// include <algorithm>

// Local includes cont'd
#include "libmesh/cell_inf_prism.h"
#include "libmesh/cell_inf_prism6.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_inf_quad4.h"

namespace libMesh
{




// ------------------------------------------------------------
// InfPrism class member functions
dof_id_type InfPrism::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:  // the triangular face at z=-1, base face

      return
        this->compute_key (this->node(0),
                           this->node(2),
                           this->node(1));

    case 1:  // the quad face at y=0

      return
        this->compute_key (this->node(0),
                           this->node(1),
                           this->node(4),
                           this->node(3));

    case 2:  // the other quad face

      return
        this->compute_key (this->node(1),
                           this->node(2),
                           this->node(5),
                           this->node(4));

    case 3: // the quad face at x=0

      return
        this->compute_key (this->node(2),
                           this->node(0),
                           this->node(3),
                           this->node(5));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
}



AutoPtr<Elem> InfPrism::side (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_sides());

  Elem* face = NULL;

  switch (i)
    {
    case 0:  // the triangular face at z=-1, base face
      {
        face = new Tri3;

        // Note that for this face element, the normal points inward
        face->set_node(0) = this->get_node(0);
        face->set_node(1) = this->get_node(1);
        face->set_node(2) = this->get_node(2);

        break;
      }

    case 1:  // the quad face at y=0
      {
        face = new InfQuad4;

        face->set_node(0) = this->get_node(0);
        face->set_node(1) = this->get_node(1);
        face->set_node(2) = this->get_node(3);
        face->set_node(3) = this->get_node(4);

        break;
      }

    case 2:  // the other quad face
      {
        face = new InfQuad4;

        face->set_node(0) = this->get_node(1);
        face->set_node(1) = this->get_node(2);
        face->set_node(2) = this->get_node(4);
        face->set_node(3) = this->get_node(5);

        break;
      }

    case 3: // the quad face at x=0
      {
        face = new InfQuad4;

        face->set_node(0) = this->get_node(2);
        face->set_node(1) = this->get_node(0);
        face->set_node(2) = this->get_node(5);
        face->set_node(3) = this->get_node(3);

        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  return AutoPtr<Elem>(face);
}



bool InfPrism::is_child_on_side(const unsigned int c,
                                const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  return (s == 0 || c+1 == s || c == s%3);
}



bool InfPrism::is_edge_on_side (const unsigned int e,
                                const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (is_node_on_side(InfPrism6::edge_nodes_map[e][0],s) &&
          is_node_on_side(InfPrism6::edge_nodes_map[e][1],s));
}



} // namespace libMesh



#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
