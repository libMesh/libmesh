// $Id: face_tri6.C,v 1.8 2003-02-06 23:02:56 benkirk Exp $

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

// C++ includes

// Local includes
#include "mesh_config.h"
#include "edge_edge3.h"
#include "face_tri6.h"
#include "mesh.h"




// ------------------------------------------------------------
// Tri6 class static member initializations
#ifdef ENABLE_AMR

const Real Tri6::embedding_matrix[4][6][6] =
{
  // embedding matrix for child 0
  {
    //  0      1      2    3    4    5
    { 1.0,   0.0,   0.0, 0.0, 0.0, 0.0}, // 0
    { 0.0,   0.0,   0.0, 1.0, 0.0, 0.0}, // 1
    { 0.0,   0.0,   0.0, 0.0, 0.0, 1.0}, // 2
    {.375, -.125,   0.0, .75, 0.0, 0.0}, // 3
    { 0.0, -.125, -.125, 0.5, .25, 0.5}, // 4
    {.375,   0.0, -.125, 0.0, 0.0, .75}  // 5
  },

  // embedding matrix for child 1
  {
    //  0      1      2    3    4    5
    {  0.0,  0.0,   0.0, 1.0, 0.0, 0.0}, // 0
    {  0.0,  1.0,   0.0, 0.0, 0.0, 0.0}, // 1
    {  0.0,  0.0,   0.0, 0.0, 1.0, 0.0}, // 2
    {-.125, .375,   0.0, .75, 0.0, 0.0}, // 3
    {  0.0, .375, -.125, 0.0, .75, 0.0}, // 4
    {-.125,  0.0, -.125, 0.5, 0.5, .25}  // 5
  },

  // embedding matrix for child 2
  {
    //  0       1     2    3    4    5
    {  0.0,   0.0,  0.0, 0.0, 0.0, 1.0}, // 0
    {  0.0,   0.0,  0.0, 0.0, 1.0, 0.0}, // 1
    {  0.0,   0.0,  1.0, 0.0, 0.0, 0.0}, // 2
    {-.125, -.125,  0.0, .25, 0.5, 0.5}, // 3
    {  0.0, -.125, .375, 0.0, .75, 0.0}, // 4
    {-.125,   0.0, .375, 0.0, 0.0, .75}  // 5
  },

  // embedding matrix for child 3
  {
    //  0       1      2    3    4    5
    {  0.0,   0.0,   0.0, 1.0, 0.0, 0.0}, // 0
    {  0.0,   0.0,   0.0, 0.0, 1.0, 0.0}, // 1
    {  0.0,   0.0,   0.0, 0.0, 0.0, 1.0}, // 2
    {-.125,   0.0, -.125, 0.5, 0.5, .25}, // 3
    {-.125, -.125,   0.0, .25, 0.5, 0.5}, // 4
    {  0.0, -.125, -.125, 0.5, .25, 0.5}  // 5
  }
};



const unsigned int Tri6::side_children_matrix[3][2] =
{
  {0, 1}, // side-0 children
  {1, 2}, // side-1 children
  {0, 2}  // side-2 children
};

#endif



// ------------------------------------------------------------
// Tri6 class member functions
AutoPtr<Elem> Tri6::build_side (const unsigned int i) const
{
  assert (i < n_sides());

  
  Edge3* edge = new Edge3;

  switch (i)
    {
    case 0:
      {
	edge->set_node(0) = get_node(0);
	edge->set_node(1) = get_node(1);
	edge->set_node(2) = get_node(3);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 1:
      {
	edge->set_node(0) = get_node(1);
	edge->set_node(1) = get_node(2);
	edge->set_node(2) = get_node(4);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 2:
      {
	edge->set_node(0) = get_node(2);
	edge->set_node(1) = get_node(0);
	edge->set_node(2) = get_node(5);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    default:
      {
	error();
      }
    };

  
  // We will never get here...  Look at the code above.
  error();
  AutoPtr<Elem> ap(NULL);  return ap;
};



const std::vector<unsigned int> Tri6::tecplot_connectivity(const unsigned int sf) const
{
  assert (_nodes != NULL);
  assert (sf < n_sub_elem());

  std::vector<unsigned int> conn(4);

  switch(sf)
    {
    case 0:
      // linear sub-triangle 0
      conn[0] = node(0)+1;
      conn[1] = node(3)+1;
      conn[2] = node(5)+1;
      conn[3] = node(5)+1;

      return conn;

    case 1:
      // linear sub-triangle 1
      conn[0] = node(3)+1;
      conn[1] = node(1)+1;
      conn[2] = node(4)+1;
      conn[3] = node(4)+1;

      return conn;

    case 2:
      // linear sub-triangle 2
      conn[0] = node(5)+1;
      conn[1] = node(4)+1;
      conn[2] = node(2)+1;
      conn[3] = node(2)+1;

      return conn;

    case 3:
      // linear sub-triangle 3
      conn[0] = node(3)+1;
      conn[1] = node(4)+1;
      conn[2] = node(5)+1;
      conn[3] = node(5)+1;

      return conn;

    default:
      error();
    };

  error();
  
  return conn;
};



void Tri6::vtk_connectivity(const unsigned int sf,
			    std::vector<unsigned int> *conn) const
{
  assert (_nodes != NULL);
  assert (sf < n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(3);
  
  switch(sf)
    {
    case 0:
      // linear sub-triangle 0
      (*conn)[0] = node(0);
      (*conn)[1] = node(3);
      (*conn)[2] = node(5);

      return;

    case 1:
      // linear sub-triangle 1
      (*conn)[0] = node(3);
      (*conn)[1] = node(1);
      (*conn)[2] = node(4);

      return;

    case 2:
      // linear sub-triangle 2
      (*conn)[0] = node(5);
      (*conn)[1] = node(4);
      (*conn)[2] = node(2);

      return;

    case 3:
      // linear sub-triangle 3
      (*conn)[0] = node(3);
      (*conn)[1] = node(4);
      (*conn)[2] = node(5);

      return;

    default:
      error();
    };

  error();
  
  return;
};



#ifdef ENABLE_AMR

void Tri6::refine(Mesh& mesh)
{
  assert (refinement_flag() == Elem::REFINE);
  assert (active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[n_children()];

    for (unsigned int c=0; c<n_children(); c++)
      {
	_children[c] = Elem::build (type());
	_children[c]->set_refinement_flag() = Elem::JUST_REFINED;
      };
  };

  // Compute new nodal locations
  // and asssign nodes to children
  {
    std::vector<std::vector<Point> >  p(n_children());
    
    for (unsigned int c=0; c<n_children(); c++)
      p[c].resize(child(c)->n_nodes());
    

    // compute new nodal locations
    for (unsigned int c=0; c<n_children(); c++)
      for (unsigned int nc=0; nc<child(c)->n_nodes(); nc++)
	for (unsigned int n=0; n<n_nodes(); n++)
	  if (embedding_matrix[c][nc][n] != 0.)
	    p[c][nc].add_scaled (point(n), embedding_matrix[c][nc][n]);
    
    
    // Assign nodes to children & add them to the mesh
    for (unsigned int c=0; c<n_children(); c++)
      {
	for (unsigned int nc=0; nc<child(c)->n_nodes(); nc++)
	  _children[c]->set_node(nc) = mesh.mesh_refinement.add_point(p[c][nc]);

	mesh.add_elem(child(c), mesh.mesh_refinement.new_element_number());
      };
  };


  
  // Possibly add boundary information
  {
    for (unsigned int s=0; s<n_sides(); s++)
      if (neighbor(s) == NULL)
	{
	  const short int id = mesh.boundary_info.boundary_id(this, s);
	
	  if (id != mesh.boundary_info.invalid_id)
	    for (unsigned int sc=0; sc<2; sc++)
	      mesh.boundary_info.add_side(child(side_children_matrix[s][sc]), s, id);
	};
  };

  
  // Un-set my refinement flag now
  set_refinement_flag() = Elem::DO_NOTHING;
};



#endif




