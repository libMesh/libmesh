// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "cell_tet.h"
#include "cell_tet4.h"
#include "face_tri3.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tet class member functions
unsigned int Tet::key (const unsigned int s) const
{
  libmesh_assert (s < this->n_sides());

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
    }

  // We'll never get here.
  libmesh_error();
  return 0;
}



AutoPtr<Elem> Tet::side (const unsigned int i) const
{
  libmesh_assert (i < this->n_sides());


  
  Elem* face = new Tri3;
  
  switch (i)
    {
    case 0:
      {
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(1);

        AutoPtr<Elem> ap_face(face);
	return ap_face;
      }
    case 1:
      {
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(3);

        AutoPtr<Elem> ap_face(face);
	return ap_face;
      }
    case 2:
      {
	face->set_node(0) = this->get_node(1);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(3);

        AutoPtr<Elem> ap_face(face);
	return ap_face;
      }
    case 3:
      {
	face->set_node(0) = this->get_node(2);
	face->set_node(1) = this->get_node(0);
	face->set_node(2) = this->get_node(3);
	
        AutoPtr<Elem> ap_face(face);
	return ap_face;
      }
    default:
      {
	libmesh_error();
      }
    }

  // We'll never get here.
  libmesh_error();
  AutoPtr<Elem> ap_face(face);
  return ap_face;
}



bool Tet::is_child_on_side(const unsigned int c,
                           const unsigned int s) const
{
  libmesh_assert (c < this->n_children());
  libmesh_assert (s < this->n_sides());

  // For the 4 vertices, child c touches vertex c, so we can return
  // true if that vertex is on side s
  for (unsigned int i = 0; i != 3; ++i)
    if (Tet4::side_nodes_map[s][i] == c)
      return true;

  // For the 4 non-vertex children, the child ordering depends on the
  // diagonal selection.  We'll let the embedding matrix figure that
  // out: if this child has three nodes that don't depend on the
  // position of the node_facing_side[s], then we're on side s.

  const unsigned int node_facing_side[4] = {3, 2, 0, 1};
  const unsigned int n = node_facing_side[s];

  unsigned int independent_nodes = 0;

  for (unsigned int nc = 0; nc != 3; ++nc)
    {
      independent_nodes++;  // Hey, we're independent so far!
      
#ifdef LIBMESH_ENABLE_AMR
      if (this->embedding_matrix(c,nc,n) != 0.)
        {
          independent_nodes--;  // No, wait, we're not
          continue;
        }
#endif //LIBMESH_ENABLE_AMR
    }

  // No subtet of an octahedron touches a side at all nodes
  libmesh_assert(independent_nodes != 4);
  // Every subtet of an octahedron touches each side at at least one
  // node
  libmesh_assert(independent_nodes != 0);

  return (independent_nodes == 3);
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
