// $Id: face_quad4.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include "edge_edge2.h"
#include "face_quad4.h"
#include "mesh.h"




// ------------------------------------------------------------
// Quad class static member initialization
#ifdef ENABLE_AMR

const real Quad4::embedding_matrix[4][4][4] =
{
  // embedding matrix for child 0
  {
    // 0    1    2    3
    {1.0, 0.0, 0.0, 0.0}, // 0
    {0.5, 0.5, 0.0, 0.0}, // 1
    {.25, .25, .25, .25}, // 2
    {0.5, 0.0, 0.0, 0.5}  // 3
  },

  // embedding matrix for child 1
  {
    // 0    1    2    3
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 1.0, 0.0, 0.0}, // 1
    {0.0, 0.5, 0.5, 0.0}, // 2
    {.25, .25, .25, .25}  // 3
  },

  // embedding matrix for child 2
  {
    // 0    1    2    3
    {0.5, 0.0, 0.0, 0.5}, // 0
    {.25, .25, .25, .25}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.0, 0.0, 1.0}  // 3
  },

  // embedding matrix for child 3
  {
    // 0    1    2    3
    {.25, .25, .25, .25}, // 0
    {0.0, 0.5, 0.5, 0.0}, // 1
    {0.0, 0.0, 1.0, 0.0}, // 2
    {0.0, 0.0, 0.5, 0.5}  // 3
  }
};



const unsigned int Quad4::side_children_matrix[4][2] =
{
  {0, 1}, // side-0 children
  {1, 3}, // side-1 children
  {2, 3}, // side-2 children
  {0, 2}  // side-3 children
};

#endif





// ------------------------------------------------------------
// Quad4 class member functions
std::auto_ptr<Elem> Quad4::build_side (const unsigned int i) const
{
  assert (i < n_sides());

  switch (i)
    {
    case 0:
      {
	Edge2* edge = new Edge2;

	edge->node(0) = node(0);
	edge->node(1) = node(1);
	
	std::auto_ptr<Elem> ap(edge);  return ap;
      }
    case 1:
      {
	Edge2* edge = new Edge2;

	edge->node(0) = node(1);
	edge->node(1) = node(2);
	
	std::auto_ptr<Elem> ap(edge);  return ap;
      }
    case 2:
      {
	Edge2* edge = new Edge2;

	edge->node(0) = node(2);
	edge->node(1) = node(3);
	
	std::auto_ptr<Elem> ap(edge);  return ap;
      }
    case 3:
      {
	Edge2* edge = new Edge2;

	edge->node(0) = node(3);
	edge->node(1) = node(0);
	
	std::auto_ptr<Elem> ap(edge);  return ap;
      }
    default:
      {
	error();
      }
    };


  // We will never get here...  Look at the code above.
  std::auto_ptr<Elem> ap(NULL);  return ap;
};



const std::vector<unsigned int> Quad4::tecplot_connectivity(const unsigned int sf) const
{
  assert (sf < n_sub_elem());
  
  std::vector<unsigned int> conn(4);

  conn[0] = node(0)+1;
  conn[1] = node(1)+1;
  conn[2] = node(2)+1;
  conn[3] = node(3)+1;

  return conn;
};







void Quad4::vtk_connectivity(const unsigned int sf,
			     std::vector<unsigned int> *conn) const
{
  assert (!_nodes.empty());
  assert (sf < n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(4);

  (*conn)[0] = node(0);
  (*conn)[1] = node(1);
  (*conn)[2] = node(2);
  (*conn)[3] = node(3);

  return;
};






unsigned int Quad4::vtk_element_type (const unsigned int sf) const
{
  return 9;
};



#ifdef ENABLE_AMR

void Quad4::refine(Mesh& mesh)
{
  assert (refinement_flag() == Elem::REFINE);
  assert (active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[n_children()];

    for (unsigned int c=0; c<n_children(); c++)
      {
	_children[c] = new Quad4(this);
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
	    p[c][nc] += mesh.vertex(node(n))*embedding_matrix[c][nc][n];
    
    
    // assign nodes to children & add them to the mesh
    for (unsigned int c=0; c<n_children(); c++)
      {
	for (unsigned int nc=0; nc<child(c)->n_nodes(); nc++)
	  _children[c]->node(nc) = mesh.mesh_refinement.add_node(p[c][nc]);

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



void Quad4::coarsen()
{
  assert (refinement_flag() == Elem::COARSEN);
  assert (!active());
  
  delete [] _children;

  _children = NULL;

  set_refinement_flag() = Elem::DO_NOTHING;
};

#endif
