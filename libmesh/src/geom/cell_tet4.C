// $Id: cell_tet4.C,v 1.9 2003-02-06 23:02:55 benkirk Exp $

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
#include "mesh.h"
#include "cell_tet4.h"
#include "face_tri3.h"



// ------------------------------------------------------------
// Tet4 class member functions
AutoPtr<Elem> Tet4::build_side (const unsigned int i) const
{
  assert (i < n_sides());

  
  
  AutoPtr<Elem> face(new Tri3);

  switch (i)
    {
    case 0:
      {
	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(2);
	face->set_node(2) = get_node(1);

	return face;
      }
    case 1:
      {
	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(1);
	face->set_node(2) = get_node(3);

	return face;
      }
    case 2:
      {
	face->set_node(0) = get_node(1);
	face->set_node(1) = get_node(2);
	face->set_node(2) = get_node(3);

	return face;
      }
    case 3:
      {
	face->set_node(0) = get_node(2);
	face->set_node(1) = get_node(0);
	face->set_node(2) = get_node(3);
	
	return face;
      }
    default:
      {
	error();
      }
    };

  // We'll never get here.
  error();  
  return face;
};



const std::vector<unsigned int> Tet4::tecplot_connectivity(const unsigned int sc) const
{
  assert (_nodes != NULL);
  assert (sc < n_sub_elem());

  std::vector<unsigned int> conn(8);
  
  conn[0] = node(0)+1;
  conn[1] = node(1)+1;
  conn[2] = node(2)+1;
  conn[3] = node(2)+1;
  conn[4] = node(3)+1;
  conn[5] = node(3)+1;
  conn[6] = node(3)+1;
  conn[7] = node(3)+1;

  return conn;
};



void Tet4::vtk_connectivity(const unsigned int sc,
			    std::vector<unsigned int> *conn) const
{
  assert (_nodes != NULL);
  assert (sc < n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(4);

  (*conn)[0] = node(0);
  (*conn)[1] = node(1);
  (*conn)[2] = node(2);
  (*conn)[3] = node(3);

  return;
};



#ifdef ENABLE_AMR

const Real Tet4::embedding_matrix[8][4][4] =
{
  // embedding matrix for child 0
  {
    // 0    1    2    3  
    {1.0, 0.0, 0.0, 0.0}, // 0
    {0.5, 0.5, 0.0, 0.0}, // 1
    {0.5, 0.0, 0.5, 0.0}, // 2
    {0.5, 0.0, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 1
  {
    // 0    1    2    3  
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 1.0, 0.0, 0.0}, // 1
    {0.0, 0.5, 0.5, 0.0}, // 2
    {0.0, 0.5, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 2
  {
    // 0    1    2    3  
    {0.5, 0.0, 0.5, 0.0}, // 0
    {0.0, 0.5, 0.5, 0.0}, // 1
    {0.0, 0.0, 1.0, 0.0}, // 2
    {0.0, 0.0, 0.5, 0.5}  // 3
  },
  
  // embedding matrix for child 3
  {
    // 0    1    2    3  
    {0.5, 0.0, 0.0, 0.5}, // 0
    {0.0, 0.5, 0.0, 0.5}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.0, 0.0, 1.0}  // 3
  },
  
  // embedding matrix for child 4
  {
    // 0    1    2    3  
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 0.5, 0.0, 0.5}, // 1
    {0.5, 0.0, 0.5, 0.0}, // 2
    {0.5, 0.0, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 5
  {
    // 0    1    2    3  
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 0.5, 0.5, 0.0}, // 1
    {0.5, 0.0, 0.5, 0.0}, // 2
    {0.0, 0.5, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 6
  {
    // 0    1    2    3  
    {0.5, 0.0, 0.5, 0.0}, // 0
    {0.0, 0.5, 0.5, 0.0}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.5, 0.0, 0.5}  // 3
  },
  
  // embedding matrix for child 7
  {
    // 0    1    2    3  
    {0.5, 0.0, 0.5, 0.0}, // 0
    {0.0, 0.5, 0.0, 0.5}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.5, 0.0, 0.0, 0.5}  // 3
  }  
};



const unsigned int Tet4::side_children_matrix[4][4] =
{
  {0, 1, 2, 5}, // side-0 children
  {0, 1, 3, 4}, // side-1 children
  {1, 2, 3, 6}, // side-2 children
  {0, 2, 3, 7}  // side-3 children
};


void Tet4::refine(Mesh& mesh)
{
  assert (refinement_flag() == Elem::REFINE);
  assert (active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[n_children()];

    for (unsigned int c=0; c<n_children(); c++)
      {
	_children[c] = new Tet4(this);
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
    
    
    // assign nodes to children & add them to the mesh
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
	    for (unsigned int sc=0; sc<4; sc++)
	      mesh.boundary_info.add_side(child(side_children_matrix[s][sc]), s, id);
	};
  };

  
  // Un-set my refinement flag now
  set_refinement_flag() = Elem::DO_NOTHING;
};


#endif
