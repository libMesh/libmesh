// $Id: cell_inf_prism12.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include "cell.h"

#ifdef ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes cont'd
#include "mesh.h"

// Temporary includes
#include "cell_inf_prism12.h"
#include "face_tri6.h"
#include "face_inf_quad6.h"




// ------------------------------------------------------------
// InfPrism12 class member functions
std::auto_ptr<Elem> InfPrism12::build_side (const unsigned int i) const
{
  assert (i < n_sides());
  assert (_nodes.size() == n_nodes());

  
  switch (i)
    {
    case 0:  // the triangular face at z=-1, base face
      {
	Tri6*  face = new Tri6;

	face->node(0) = node(0);
	face->node(1) = node(2);
	face->node(2) = node(1);
	face->node(3) = node(8);
	face->node(4) = node(7);
	face->node(5) = node(6);

	std::auto_ptr<Elem> ap(face);  return ap;
      }

    case 1:  // the quad face at y=0
      {
	InfQuad6* face = new InfQuad6;
	
	face->node(0) = node(0);
	face->node(1) = node(1);
	face->node(2) = node(4);
	face->node(3) = node(3);
	face->node(4) = node(6);
	face->node(5) = node(9);
	
	std::auto_ptr<Elem> ap(face);  return ap;
      }

    case 2:  // the other quad face
      {
	InfQuad6* face = new InfQuad6;

	face->node(0) = node(1);
	face->node(1) = node(2);
	face->node(2) = node(5);
	face->node(3) = node(4);
	face->node(4) = node(7);
	face->node(5) = node(10);

	std::auto_ptr<Elem> ap(face);  return ap;
      }

    case 3: // the quad face at x=0
      {
	InfQuad6* face = new InfQuad6;

	face->node(0) = node(2);
	face->node(1) = node(0);
	face->node(2) = node(3);
	face->node(3) = node(5);
	face->node(4) = node(8);
	face->node(5) = node(11);
	
	std::auto_ptr<Elem> ap(face);  return ap;
      }

    case 4: // the triangular face at z=1
      {
        std::cerr << "No face 4 in case of infinite elements!" << std::endl;
        error();
	std::auto_ptr<Elem> ap(NULL);  return ap;
      }

    default:
      {
	error();
	std::auto_ptr<Elem> ap(NULL);  return ap;
      }
    };

  // We'll never get here.
  error();
  std::auto_ptr<Elem> ap(NULL);  return ap;
};



const std::vector<unsigned int> InfPrism12::tecplot_connectivity(const unsigned int sc) const
{
  assert (!_nodes.empty());
  assert (sc < n_sub_elem());

  std::vector<unsigned int> conn(8);

  switch (sc)
    {
    case 0:

      // guess this is a collapsed hex8
      conn[0] = node(0)+1;
      conn[1] = node(6)+1;
      conn[2] = node(8)+1;
      conn[3] = node(8)+1;
      conn[4] = node(3)+1;
      conn[5] = node(9)+1;
      conn[6] = node(11)+1;
      conn[7] = node(11)+1;

      return conn;

    case 1:

      conn[0] = node(6)+1;
      conn[1] = node(7)+1;
      conn[2] = node(8)+1;
      conn[3] = node(8)+1;
      conn[4] = node(9)+1;
      conn[5] = node(10)+1;
      conn[6] = node(11)+1;
      conn[7] = node(11)+1;

      return conn;

    case 2:

      conn[0] = node(6)+1;
      conn[1] = node(1)+1;
      conn[2] = node(7)+1;
      conn[3] = node(7)+1;
      conn[4] = node(9)+1;
      conn[5] = node(4)+1;
      conn[6] = node(10)+1;
      conn[7] = node(10)+1;

      return conn;

    case 3:

      conn[0] = node(8)+1;
      conn[1] = node(7)+1;
      conn[2] = node(2)+1;
      conn[3] = node(2)+1;
      conn[4] = node(11)+1;
      conn[5] = node(10)+1;
      conn[6] = node(5)+1;
      conn[7] = node(5)+1;

      return conn;

    default:
      error();
      
    };

  error();
  return conn;
};



void InfPrism12::write_tecplot_connectivity(std::ostream &out) const
{
  assert (out);
  assert (!_nodes.empty());

  for (unsigned int sc=0; sc<n_sub_elem(); sc++)
    {
      std::vector<unsigned int> conn = tecplot_connectivity(sc);

      for (unsigned int i=0; i<conn.size(); i++)
	out << conn[i] << " ";

      out << std::endl;
    };
};



#ifdef ENABLE_AMR

const real InfPrism12::embedding_matrix[4][12][12] =
{
  // embedding matrix for child 0
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11 th parent Node
    {         1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 5
    {       0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 6
    {         0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.5,       0.25,        0.5,        0.0,        0.0,        0.0}, // 7
    {       0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 8
    {         0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 9
    {         0.0,        0.0,        0.0,        0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.5,       0.25,        0.5}, // 10
    {         0.0,        0.0,        0.0,      0.375,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75}  // 11
  },

  // embedding matrix for child 1
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11 th parent Node
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 5
    {      -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0,        0.0}, // 6
    {         0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 7
    {      -0.125,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.5,        0.5,       0.25,        0.0,        0.0,        0.0}, // 8
    {         0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0}, // 9
    {         0.0,        0.0,        0.0,        0.0,      0.375,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 10
    {         0.0,        0.0,        0.0,     -0.125,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.5,        0.5,       0.25}  // 11
  },

  // embedding matrix for child 2
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11 th parent Node
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 5
    {      -0.125,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.25,        0.5,        0.5,        0.0,        0.0,        0.0}, // 6
    {         0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0,        0.0}, // 7
    {      -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0,        0.0,        0.0}, // 8
    {         0.0,        0.0,        0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.25,        0.5,        0.5}, // 9
    {         0.0,        0.0,        0.0,        0.0,     -0.125,      0.375,        0.0,        0.0,        0.0,        0.0,       0.75,        0.0}, // 10
    {         0.0,        0.0,        0.0,     -0.125,        0.0,      0.375,        0.0,        0.0,        0.0,        0.0,        0.0,       0.75}  // 11
  },

  // embedding matrix for child 3
  {
    //          0           1           2           3           4           5           6           7           8           9          10          11 th parent Node
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        0.0,        1.0}, // 5
    {      -0.125,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.5,        0.5,       0.25,        0.0,        0.0,        0.0}, // 6
    {      -0.125,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.25,        0.5,        0.5,        0.0,        0.0,        0.0}, // 7
    {         0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.5,       0.25,        0.5,        0.0,        0.0,        0.0}, // 8
    {         0.0,        0.0,        0.0,     -0.125,        0.0,     -0.125,        0.0,        0.0,        0.0,        0.5,        0.5,       0.25}, // 9
    {         0.0,        0.0,        0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.0,       0.25,        0.5,        0.5}, // 10
    {         0.0,        0.0,        0.0,        0.0,     -0.125,     -0.125,        0.0,        0.0,        0.0,        0.5,       0.25,        0.5}  // 11
  },


};



const unsigned int InfPrism12::side_children_matrix[5][5] =
{
  {4,   0, 1, 2, 3}, // 4 side-0 children
  {2,   0, 1,42,42}, // 2 side-1 children
  {2,   1, 2,42,42}, // 2 side-2 children
  {2,   0, 2,42,42}, // 2 side-3 children
  {4,   0, 1, 2, 3}  // 4 side-4 children
};



void InfPrism12::refine(Mesh& mesh)
{
  assert (refinement_flag() == Elem::REFINE);
  assert (active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[n_children()];

    for (unsigned int c=0; c<n_children(); c++)
      {
	_children[c] = new InfPrism12(this);
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
	    // the upper limit for sc is stored in the 0th column
	    for (unsigned int sc=1; sc<=side_children_matrix[s][0]; sc++)
	      mesh.boundary_info.add_side(child(side_children_matrix[s][sc]), s, id);

	};
  };


  // Un-set my refinement flag now
  set_refinement_flag() = Elem::DO_NOTHING;
};



void InfPrism12::coarsen()
{
  assert (refinement_flag() == Elem::COARSEN);
  assert (!active());
  
  delete [] _children;

  _children = NULL;

  set_refinement_flag() = Elem::DO_NOTHING;
};


#endif

#endif
