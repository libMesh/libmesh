// $Id: face_inf_quad6.C,v 1.13 2003-03-03 02:15:58 benkirk Exp $

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
#include "mesh_base.h"
#include "face_inf_quad6.h"
#include "edge_edge3.h"
#include "edge_inf_edge2.h"




// ------------------------------------------------------------
// InfQuad6 class static member initialization
#ifdef ENABLE_AMR

const float InfQuad6::_embedding_matrix[2][6][6] =
{
  // embedding matrix for child 0
  {
    //     0       1       2       3       4       5th parent node
    {    1.0,    0.0,    0.0,    0.0,    0.0,    0.0 }, // 0th child node
    {    0.0,    0.0,    0.0,    0.0,    1.0,    0.0 }, // 1
    {    0.0,    0.0,    0.0,    0.0,    0.0,    1.0 }, // 2
    {    0.0,    0.0,    0.0,    1.0,    0.0,    0.0 }, // 3
    {  0.375, -0.125,    0.0,    0.0,   0.75,    0.0 }, // 4
    {    0.0,    0.0, -0.125,  0.375,    0.0,   0.75 }  // 5
  },

  // embedding matrix for child 1
  {
    //     0       1       2       3       4       5th parent node
    {    0.0,    0.0,    0.0,    0.0,    1.0,    0.0 }, // 0th child node
    {    0.0,    1.0,    0.0,    0.0,    0.0,    0.0 }, // 1
    {    0.0,    0.0,    1.0,    0.0,    0.0,    0.0 }, // 2
    {    0.0,    0.0,    0.0,    0.0,    0.0,    1.0 }, // 3
    { -0.125,  0.375,    0.0,    0.0,   0.75,    0.0 }, // 4
    {    0.0,    0.0,  0.375, -0.125,    0.0,   0.75 }  // 5
  }
};



const unsigned int InfQuad6::_side_children_matrix[4][3] =
{
  // note different storage scheme
  {2,   0, 1}, // 2 side-0 children
  {1,   1,42}, // 1 side-1 children
  {2,   0, 1}, // 2 side-2 children
  {1,   0,42}  // 1 side-3 children
};

#endif



// ------------------------------------------------------------
// InfQuad6 class member functions
AutoPtr<Elem> InfQuad6::build_side (const unsigned int i) const
{
  assert (i < this->n_sides());

  
  switch (i)
    {
    case 0:
      {
	Edge3* edge = new Edge3;

	edge->set_node(0) = this->get_node(0);
	edge->set_node(1) = this->get_node(1);
	edge->set_node(2) = this->get_node(4);
	
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



const std::vector<unsigned int> InfQuad6::tecplot_connectivity(const unsigned int sf) const
{
  assert (_nodes != NULL);
  assert (sf < this->n_sub_elem());

  std::vector<unsigned int> conn(4);

  switch(sf)
    {
    case 0:
      // linear sub-quad 0
      conn[0] = this->node(0)+1;
      conn[1] = this->node(4)+1;
      conn[2] = this->node(5)+1;
      conn[3] = this->node(3)+1;

      return conn;

    case 1:
      // linear sub-quad 1
      conn[0] = this->node(4)+1;
      conn[1] = this->node(1)+1;
      conn[2] = this->node(2)+1;
      conn[3] = this->node(5)+1;

      return conn;

    default:
      error();

    }

  error();
  
  return conn;
}



#ifdef ENABLE_AMR

void InfQuad6::refine (MeshBase& mesh)
{
  assert (this->refinement_flag() == Elem::REFINE);
  assert (this->active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[this->n_children()];

    for (unsigned int c=0; c<this->n_children(); c++)
      {
	_children[c] = new InfQuad6(this);
	_children[c]->set_refinement_flag(Elem::JUST_REFINED);
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
	  if (_embedding_matrix[c][nc][n] != 0.)
	    p[c][nc].add_scaled (this->point(n), static_cast<Real>(_embedding_matrix[c][nc][n]));
    
    
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
	    for (unsigned int sc=0; sc < 2; sc++)
	      mesh.boundary_info.add_side(this->child(_side_children_matrix[s][sc]), s, id);
	}
  }

  
  // Un-set my refinement flag now
  this->set_refinement_flag(Elem::DO_NOTHING);
}



#endif // #ifdef ENABLE_AMR


#endif // #ifdef ENABLE_INFINITE_ELEMENTS
