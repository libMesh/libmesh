// $Id: face_quad9.C,v 1.9 2003-02-13 22:56:12 benkirk Exp $

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
#include "face_quad9.h"
#include "mesh.h"




// ------------------------------------------------------------
// Quad9 class static member initializations
#ifdef ENABLE_AMR

const Real Quad9::embedding_matrix[4][9][9] =
{
  // embedding matrix for child 0
  {
    //         0           1           2           3           4           5           6           7           8
    {    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
    {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 3
    {   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
    {    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 5
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,   0.750000 }, // 6
    {   0.375000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 7
    {   0.140625, -0.0468750,  0.0156250, -0.0468750,   0.281250, -0.0937500, -0.0937500,   0.281250,   0.562500 }  // 8
  },

  // embedding matrix for child 1
  {
    //         0           1           2           3           4           5           6           7           8
    {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
    {    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 3
    {  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
    {    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 5
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,   0.750000 }, // 6
    {    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 7
    { -0.0468750,   0.140625, -0.0468750,  0.0156250,   0.281250,   0.281250, -0.0937500, -0.0937500,   0.562500 }  // 8
  },

  // embedding matrix for child 2
  {
    //         0           1           2           3           4           5           6           7           8
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 0
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 1
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,   0.750000 }, // 4
    {    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 5
    {    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 6
    {  -0.125000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 7
    { -0.0468750,  0.0156250, -0.0468750,   0.140625, -0.0937500, -0.0937500,   0.281250,   0.281250,   0.562500 }  // 8
  },

  // embedding matrix for child 3
  {
    //         0           1           2           3           4           5           6           7           8
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 0
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 1
    {    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 3
    {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,   0.750000 }, // 4
    {    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 5
    {    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 6
    {    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 7
    {  0.0156250, -0.0468750,   0.140625, -0.0468750, -0.0937500,   0.281250,   0.281250, -0.0937500,   0.562500 }  // 8
  }
};



const unsigned int Quad9::side_children_matrix[4][2] =
{
  {0, 1}, // side-0 children
  {1, 3}, // side-1 children
  {2, 3}, // side-2 children
  {0, 2}  // side-3 children
};

#endif



// ------------------------------------------------------------
// Quad9 class member functions
AutoPtr<Elem> Quad9::build_side (const unsigned int i) const
{
  assert (i < n_sides());

  
  Edge3* edge = new Edge3;

  switch (i)
    {
    case 0:
      {
	edge->set_node(0) = get_node(0);
	edge->set_node(1) = get_node(1);
	edge->set_node(2) = get_node(4);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 1:
      {
	edge->set_node(0) = get_node(1);
	edge->set_node(1) = get_node(2);
	edge->set_node(2) = get_node(5);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 2:
      {
	edge->set_node(0) = get_node(2);
	edge->set_node(1) = get_node(3);
	edge->set_node(2) = get_node(6);
	
	AutoPtr<Elem> ap(edge);  return ap;
      }
    case 3:
      {
	edge->set_node(0) = get_node(3);
	edge->set_node(1) = get_node(0);
	edge->set_node(2) = get_node(7);
	
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



const std::vector<unsigned int> Quad9::tecplot_connectivity(const unsigned int sf) const
{
  assert (_nodes != NULL);
  assert (sf < n_sub_elem());

  std::vector<unsigned int> conn(4);

  switch(sf)
    {
    case 0:
      // linear sub-quad 0
      conn[0] = node(0)+1;
      conn[1] = node(4)+1;
      conn[2] = node(8)+1;
      conn[3] = node(7)+1;

      return conn;

    case 1:
      // linear sub-quad 1
      conn[0] = node(4)+1;
      conn[1] = node(1)+1;
      conn[2] = node(5)+1;
      conn[3] = node(8)+1;

      return conn;

    case 2:
      // linear sub-quad 2
      conn[0] = node(7)+1;
      conn[1] = node(8)+1;
      conn[2] = node(6)+1;
      conn[3] = node(3)+1;

      return conn;

    case 3:
      // linear sub-quad 3
      conn[0] = node(8)+1;
      conn[1] = node(5)+1;
      conn[2] = node(2)+1;
      conn[3] = node(6)+1;

      return conn;

    default:
      error();
    }

  error();
  
  return conn;
}






void Quad9::vtk_connectivity(const unsigned int sf,
			     std::vector<unsigned int> *conn) const
{
  assert (_nodes != NULL);
  assert (sf < n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(4);

  switch(sf)
    {
    case 0:
      // linear sub-quad 0
      (*conn)[0] = node(0);
      (*conn)[1] = node(4);
      (*conn)[2] = node(8);
      (*conn)[3] = node(7);

      return;

    case 1:
      // linear sub-quad 1
      (*conn)[0] = node(4);
      (*conn)[1] = node(1);
      (*conn)[2] = node(5);
      (*conn)[3] = node(8);

      return;

    case 2:
      // linear sub-quad 2
      (*conn)[0] = node(7);
      (*conn)[1] = node(8);
      (*conn)[2] = node(6);
      (*conn)[3] = node(3);

      return;

    case 3:
      // linear sub-quad 3
      (*conn)[0] = node(8);
      (*conn)[1] = node(5);
      (*conn)[2] = node(2);
      (*conn)[3] = node(6);

      return;

    default:
      error();
    }

  error();
  
  return;
}



#ifdef ENABLE_AMR

void Quad9::refine(Mesh& mesh)
{
  assert (refinement_flag() == Elem::REFINE);
  assert (active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[n_children()];

    for (unsigned int c=0; c<n_children(); c++)
      {
	_children[c] = new Quad9(this);
	_children[c]->set_refinement_flag() = Elem::JUST_REFINED;
      }
  }

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
    
    
    // assign nodes to children & add them to the mesh
    for (unsigned int c=0; c<n_children(); c++)
      {
	for (unsigned int nc=0; nc<child(c)->n_nodes(); nc++)
	  _children[c]->set_node(nc) = mesh.mesh_refinement.add_point(p[c][nc]);

	mesh.add_elem(child(c), mesh.mesh_refinement.new_element_number());
      }
  }


  
  // Possibly add boundary information
  {
    for (unsigned int s=0; s<n_sides(); s++)
      if (neighbor(s) == NULL)
	{
	  const short int id = mesh.boundary_info.boundary_id(this, s);
	
	  if (id != mesh.boundary_info.invalid_id)
	    for (unsigned int sc=0; sc<2; sc++)
	      mesh.boundary_info.add_side(child(side_children_matrix[s][sc]), s, id);
	}
  }

  
  // Un-set my refinement flag now
  set_refinement_flag() = Elem::DO_NOTHING;
}



#endif

