// $Id: face_inf_quad4.C,v 1.11 2003-02-26 04:43:14 jwpeterson Exp $

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
#include "mesh_common.h"

#ifdef ENABLE_INFINITE_ELEMENTS


// Local includes cont'd
#include "mesh.h"
#include "face_inf_quad4.h"
#include "edge_edge2.h"
#include "edge_inf_edge2.h"




// ------------------------------------------------------------
// InfQuad4 class static member initialization
#ifdef ENABLE_AMR

const float  InfQuad4::embedding_matrix[2][4][4] =
{
  // embedding matrix for child 0
  {
    // 0    1    2    3rd parent node
    {1.0, 0.0, 0.0, 0.0}, // 0th child node
    {0.5, 0.5, 0.0, 0.0}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.0, 0.0, 1.0}  // 3
  },

  // embedding matrix for child 1
  {
    // 0    1    2    3
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 1.0, 0.0, 0.0}, // 1
    {0.0, 0.0, 1.0, 0.0}, // 2
    {0.0, 0.0, 0.5, 0.5}  // 3
  }
};



const unsigned int InfQuad4::side_children_matrix[4][3] =
{
  // note different storage scheme
  {2,   0, 1}, // 2 side-0 children
  {1,   1,42}, // 1 side-1 children
  {2,   0, 1}, // 2 side-2 children
  {1,   0,42}  // 1 side-3 children
};

#endif



// ------------------------------------------------------------
// InfQuad4 class member functions
AutoPtr<Elem> InfQuad4::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  switch (i)
    {
    case 0:
      {
	// base face
	Edge2* edge = new Edge2;

	edge->set_node(0) = this->get_node(0);
	edge->set_node(1) = this->get_node(1);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 1:
      {
	// adjacent to another infinite element
        InfEdge2* edge = new InfEdge2;

	edge->set_node(0) = this->get_node(1);
	edge->set_node(1) = this->get_node(2);	

	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 2:
      {
	// supposed to lie at infinity
	std::cerr << "Side represents the exterior. No face." << std::endl;
	error(); 
     }
    case 3:
      {
	// adjacent to another infinite element	
	InfEdge2* edge = new InfEdge2;

	edge->set_node(0) = this->get_node(0); // be aware of swapped nodes,
	edge->set_node(1) = this->get_node(3); // compared to conventional side numbering

	AutoPtr<Elem> ap(edge);  return ap;
      }
    default:
      {
	error();
      }
    }


  // We will never get here...  Look at the code above.
  error();
  AutoPtr<Elem> ap(NULL);  return ap;
}







const std::vector<unsigned int> InfQuad4::tecplot_connectivity(const unsigned int sf) const
{
  assert (sf < this->n_sub_elem());
  
  std::vector<unsigned int> conn(4);

  conn[0] = this->node(0)+1;
  conn[1] = this->node(1)+1;
  conn[2] = this->node(2)+1;
  conn[3] = this->node(3)+1;

  return conn;
}



#ifdef ENABLE_AMR

void InfQuad4::refine(Mesh& mesh)
{
  assert (this->refinement_flag() == Elem::REFINE);
  assert (this->active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[this->n_children()];

    for (unsigned int c=0; c<this->n_children(); c++)
      {
	_children[c] = new InfQuad4(this);
	_children[c]->set_refinement_flag() = Elem::JUST_REFINED;
      }
  }

  // Compute new nodal locations
  // and asssign nodes to children
  {
    std::vector<std::vector<Point> >  p(this->n_children());
    
    for (unsigned int c=0; c<this->n_children(); c++)
      p[c].resize(this->child(c)->n_nodes());
    

    // compute new nodal locations
    for (unsigned int c=0; c<this->n_children(); c++)
      for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
	for (unsigned int n=0; n<this->n_nodes(); n++)
	  if (embedding_matrix[c][nc][n] != 0.)
	    p[c][nc].add_scaled (this->point(n), static_cast<Real>(embedding_matrix[c][nc][n]));
    
    
    // assign nodes to children & add them to the mesh
    for (unsigned int c=0; c<this->n_children(); c++)
      {
	for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
	  _children[c]->set_node(nc) = mesh.mesh_refinement.add_point(p[c][nc]);

	mesh.add_elem(this->child(c), mesh.mesh_refinement.new_element_number());
      }
  }


  
  // Possibly add boundary information
  {
    for (unsigned int s=0; s<this->n_sides(); s++)
      if (this->neighbor(s) == NULL)
	{
	  const short int id = mesh.boundary_info.boundary_id(this, s);
	
	  if (id != mesh.boundary_info.invalid_id)
	    // the upper limit for sc is stored in the 0th column
	    for (unsigned int sc=1; sc <=side_children_matrix[s][0]; sc++)
	      mesh.boundary_info.add_side(this->child(side_children_matrix[s][sc]), s, id);


	}
  }

  
  // Un-set my refinement flag now
  this->set_refinement_flag() = Elem::DO_NOTHING;
}



#endif


#endif // #ifdef ENABLE_INFINITE_ELEMENTS
