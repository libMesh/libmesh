// $Id: mesh_base_modification.C,v 1.18 2004-11-12 22:36:09 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include <algorithm> // for std::count
#include <map>
#include <math.h>    // for fabs

// Local includes
#include "mesh_base.h"
#include "face_tri3.h"
#include "face_tri6.h"


// ------------------------------------------------------------
// MeshBase class member functions for mesh modification
void MeshBase::all_tri ()
{
  // assert (this->mesh_dimension() == 2);
	  
  std::vector<Elem*> new_elements;
  new_elements.reserve (2*this->n_active_elem());

//   active_elem_iterator el (this->elements_begin());
//   active_elem_iterator end(this->elements_end());

  element_iterator       el  = this->active_elements_begin();
  const element_iterator end = this->active_elements_end(); 

  for (; el!=end; ++el)
    {
      if ((*el)->type() == QUAD4)
	{
	  Elem* tri0 = new Tri3;
	  Elem* tri1 = new Tri3;
	  
	  // Check for possible edge swap
	  if (((*el)->point(0) - (*el)->point(2)).size() <
	      ((*el)->point(1) - (*el)->point(3)).size())
	    {	      
	      tri0->set_node(0) = (*el)->get_node(0);
	      tri0->set_node(1) = (*el)->get_node(1);
	      tri0->set_node(2) = (*el)->get_node(2);
	      
	      tri1->set_node(0) = (*el)->get_node(0);
	      tri1->set_node(1) = (*el)->get_node(2);
	      tri1->set_node(2) = (*el)->get_node(3);
	    }

	  else
	    {
	      tri0->set_node(0) = (*el)->get_node(0);
	      tri0->set_node(1) = (*el)->get_node(1);
	      tri0->set_node(2) = (*el)->get_node(3);
	      
	      tri1->set_node(0) = (*el)->get_node(1);
	      tri1->set_node(1) = (*el)->get_node(2);
	      tri1->set_node(2) = (*el)->get_node(3);
	    }
	  
	  new_elements.push_back(tri0);
	  new_elements.push_back(tri1);
	  
	  delete *el;
	}
      
      else if ((*el)->type() == QUAD8)
	{
	  Elem* tri0 = new Tri6;
	  Elem* tri1 = new Tri6;
	  
	  Node* new_node = add_point((node((*el)->node(0)) +
				      node((*el)->node(1)) +
				      node((*el)->node(2)) +
				      node((*el)->node(3)))*.25
				     );
	  
	  // Check for possible edge swap
	  if (((*el)->point(0) - (*el)->point(2)).size() <
	      ((*el)->point(1) - (*el)->point(3)).size())
	    {	      
	      tri0->set_node(0) = (*el)->get_node(0);
	      tri0->set_node(1) = (*el)->get_node(1);
	      tri0->set_node(2) = (*el)->get_node(2);
	      tri0->set_node(3) = (*el)->get_node(4);
	      tri0->set_node(4) = (*el)->get_node(5);
	      tri0->set_node(5) = new_node;
	      
	      tri1->set_node(0) = (*el)->get_node(0);
	      tri1->set_node(1) = (*el)->get_node(2);
	      tri1->set_node(2) = (*el)->get_node(3);
	      tri1->set_node(3) = new_node;
	      tri1->set_node(4) = (*el)->get_node(6);
	      tri1->set_node(5) = (*el)->get_node(7);

	    }
	  
	  else
	    {
	      tri0->set_node(0) = (*el)->get_node(3);
	      tri0->set_node(1) = (*el)->get_node(0);
	      tri0->set_node(2) = (*el)->get_node(1);
	      tri0->set_node(3) = (*el)->get_node(7);
	      tri0->set_node(4) = (*el)->get_node(4);
	      tri0->set_node(5) = new_node;
	      
	      tri1->set_node(0) = (*el)->get_node(1);
	      tri1->set_node(1) = (*el)->get_node(2);
	      tri1->set_node(2) = (*el)->get_node(3);
	      tri1->set_node(3) = (*el)->get_node(5);
	      tri1->set_node(4) = (*el)->get_node(6);
	      tri1->set_node(5) = new_node;
	    }
	  
	  new_elements.push_back(tri0);
	  new_elements.push_back(tri1);
	  
	  delete *el;
	}
      
      else if ((*el)->type() == QUAD9)
	{
	  Elem* tri0 = new Tri6;
	  Elem* tri1 = new Tri6;

	  // Check for possible edge swap
	  if (((*el)->point(0) - (*el)->point(2)).size() <
	      ((*el)->point(1) - (*el)->point(3)).size())
	    {	      
	      tri0->set_node(0) = (*el)->get_node(0);
	      tri0->set_node(1) = (*el)->get_node(1);
	      tri0->set_node(2) = (*el)->get_node(2);
	      tri0->set_node(3) = (*el)->get_node(4);
	      tri0->set_node(4) = (*el)->get_node(5);
	      tri0->set_node(5) = (*el)->get_node(8);
	      
	      tri1->set_node(0) = (*el)->get_node(0);
	      tri1->set_node(1) = (*el)->get_node(2);
	      tri1->set_node(2) = (*el)->get_node(3);
	      tri1->set_node(3) = (*el)->get_node(8);
	      tri1->set_node(4) = (*el)->get_node(6);
	      tri1->set_node(5) = (*el)->get_node(7);
	    }

	  else
	    {
	      tri0->set_node(0) = (*el)->get_node(0);
	      tri0->set_node(1) = (*el)->get_node(1);
	      tri0->set_node(2) = (*el)->get_node(3);
	      tri0->set_node(3) = (*el)->get_node(4);
	      tri0->set_node(4) = (*el)->get_node(8);
	      tri0->set_node(5) = (*el)->get_node(7);
	      
	      tri1->set_node(0) = (*el)->get_node(1);
	      tri1->set_node(1) = (*el)->get_node(2);
	      tri1->set_node(2) = (*el)->get_node(3);
	      tri1->set_node(3) = (*el)->get_node(5);
	      tri1->set_node(4) = (*el)->get_node(6);
	      tri1->set_node(5) = (*el)->get_node(8);
	    }
	  
	  new_elements.push_back(tri0);
	  new_elements.push_back(tri1);

	  delete *el;
	}
      else
	new_elements.push_back(*el);
    }
  
  _elements = new_elements;

  this->prepare_for_use();
}






