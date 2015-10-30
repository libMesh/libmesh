// $Id: face_tri3.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include "face_tri3.h"
#include "mesh.h"




// ------------------------------------------------------------
// Tri3 class static member initializations
#ifdef ENABLE_AMR

const real Tri3::embedding_matrix[4][3][3] =
{
  // embedding matrix for child 0
  {
    // 0    1    2  
    {1.0, 0.0, 0.0}, // 0
    {0.5, 0.5, 0.0}, // 1
    {0.5, 0.0, 0.5}  // 2
  },

  // embedding matrix for child 1
  {
    // 0    1    2  
    {0.5, 0.5, 0.0}, // 0
    {0.0, 1.0, 0.0}, // 1
    {0.0, 0.5, 0.5}  // 2
  },

  // embedding matrix for child 2
  {
    // 0    1    2  
    {0.5, 0.0, 0.5}, // 0
    {0.0, 0.5, 0.5}, // 1
    {0.0, 0.0, 1.0}  // 2
  },

  // embedding matrix for child 3
  {
    // 0    1    2  
    {0.5, 0.5, 0.0}, // 0
    {0.0, 0.5, 0.5}, // 1
    {0.5, 0.0, 0.5}  // 2
  }
};



const unsigned int Tri3::side_children_matrix[3][2] =
{
  {0, 1}, // side-0 children
  {1, 2}, // side-1 children
  {0, 2}  // side-2 children
};

#endif



// ------------------------------------------------------------
// Tri3 class member functions
std::auto_ptr<Elem> Tri3::build_side (const unsigned int i) const
{
  assert (i < n_sides());

  
  Edge2* edge = new Edge2;

  switch (i)
    {
    case 0:
      {
	edge->node(0) = node(0);
	edge->node(1) = node(1);
	
	std::auto_ptr<Elem> ap(edge);  return ap;
      }
    case 1:
      {
	edge->node(0) = node(1);
	edge->node(1) = node(2);
	
	std::auto_ptr<Elem> ap(edge);  return ap;
      }
    case 2:
      {
	edge->node(0) = node(2);
	edge->node(1) = node(0);
	
	std::auto_ptr<Elem> ap(edge);  return ap;
      }
    default:
      {
	error();
      }
    };

  
  // We will never get here...  Look at the code above.
  error();
  std::auto_ptr<Elem> ap(NULL);  return ap;
};





const std::vector<unsigned int> Tri3::tecplot_connectivity(const unsigned int sf) const
{
  assert (sf <n_sub_elem());
  
  std::vector<unsigned int> conn(4);
  
  conn[0] = node(0)+1;
  conn[1] = node(1)+1;
  conn[2] = node(2)+1;
  conn[3] = node(2)+1;

  return conn;
};





void Tri3::vtk_connectivity(const unsigned int sf,
			    std::vector<unsigned int> *conn) const
{
  assert (!_nodes.empty());
  assert (sf < n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(3);

  (*conn)[0] = node(0);
  (*conn)[1] = node(1);
  (*conn)[2] = node(2);

  return;
};





unsigned int Tri3::vtk_element_type (const unsigned int sc) const
{
  return 5;
};



#ifdef ENABLE_AMR

void Tri3::refine(Mesh& mesh)
{
  assert (refinement_flag() == Elem::REFINE);
  assert (active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[n_children()];

    for (unsigned int c=0; c<n_children(); c++)
      {
	_children[c] = new Tri3(this);
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




void Tri3::coarsen()
{
  assert (refinement_flag() == Elem::COARSEN);
  assert (!active());
  
  delete [] _children;

  _children = NULL;

  set_refinement_flag() = Elem::DO_NOTHING;
};

#endif 


