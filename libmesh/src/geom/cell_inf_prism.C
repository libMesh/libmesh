// $Id: cell_inf_prism.C,v 1.3 2003-07-12 16:33:18 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS


// C++ includes
// include <algorithm>

// Local includes cont'd
#include "cell_inf_prism.h"
#include "face_tri3.h"
#include "face_inf_quad4.h"




// ------------------------------------------------------------
// InfPrism class member functions
unsigned int InfPrism::key (const unsigned int s) const
{
  assert (s < this->n_sides());
  
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
    }

  // We'll never get here.
  error();
  return 0;
}



AutoPtr<Elem> InfPrism::side (const unsigned int i) const
{
  assert (i < this->n_sides());
  
  switch (i)
    {
    case 0:  // the triangular face at z=-1, base face
      {
	AutoPtr<Elem> face(new Tri3);

	// Note that for this face element, the normal points inward
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(2);

	return face;
      }

    case 1:  // the quad face at y=0
      {
	AutoPtr<Elem> face(new InfQuad4);
	
	face->set_node(0) = this->get_node(0);
	face->set_node(1) = this->get_node(1);
	face->set_node(2) = this->get_node(3);
	face->set_node(3) = this->get_node(4);
	
	return face;
      }

    case 2:  // the other quad face
      {
	AutoPtr<Elem> face(new InfQuad4);

	face->set_node(0) = this->get_node(1);
	face->set_node(1) = this->get_node(2);
	face->set_node(2) = this->get_node(4);
	face->set_node(3) = this->get_node(5);

	return face;
      }

    case 3: // the quad face at x=0
      {
	AutoPtr<Elem> face(new InfQuad4);

	face->set_node(0) = this->get_node(2);
	face->set_node(1) = this->get_node(0);
	face->set_node(2) = this->get_node(5);
	face->set_node(3) = this->get_node(3);
	
	return face;
      }

    default:
      {
	error(); 
	AutoPtr<Elem> ap(NULL);  return ap;
      }
    }

  // We'll never get here.
  error();
  AutoPtr<Elem> ap(NULL);  return ap;
}




#ifdef ENABLE_AMR

const unsigned int InfPrism::_side_children_matrix[4][4] =
{
  {0, 1,  2,  3}, // 4 side-0 children
  {0, 1, 42, 42}, // 2 side-1 children
  {1, 2, 42, 42}, // 2 side-2 children
  {0, 2, 42, 42}, // 2 side-3 children
};

#endif



#endif // ifdef ENABLE_INFINITE_ELEMENTS
