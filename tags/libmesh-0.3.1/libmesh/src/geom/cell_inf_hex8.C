// $Id: cell_inf_hex8.C,v 1.11 2003-02-20 23:18:11 benkirk Exp $

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

// Local includes cont'd
#include "mesh.h"
#include "cell_inf_hex8.h"
#include "face_quad4.h"
#include "face_inf_quad4.h"


// ------------------------------------------------------------
// InfHex8 class member functions
AutoPtr<Elem> InfHex8::build_side (const unsigned int i) const
{
  assert (i < n_sides());



  // Think of a unit cube: (-1,1) x (-1,1)x (-1,1),
  // with the normals pointing outwards, see
  // also documentation of InfHex8 class.
  switch (i)
    {
    case 0:  // the face at z = -1
      // the base, where the infinite element couples to conventional
      // elements
      {  
        AutoPtr<Elem> face(new Quad4);

	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(1);
	face->set_node(2) = get_node(2);
	face->set_node(3) = get_node(3);
/* old code
	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(3);
	face->set_node(2) = get_node(2);
	face->set_node(3) = get_node(1);
*/

	return face;
      }
    case 1:  // the face at y = -1
      // this face connects to another infinite element
      {
	AutoPtr<Elem> face(new InfQuad4);

	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(1);
	face->set_node(2) = get_node(5);
	face->set_node(3) = get_node(4);
	
	return face;
      }
    case 2:  // the face at x = 1
      // this face connects to another infinite element
      {
	AutoPtr<Elem> face(new InfQuad4);

	face->set_node(0) = get_node(1);
	face->set_node(1) = get_node(2);
	face->set_node(2) = get_node(6);
	face->set_node(3) = get_node(5);

	return face;
      }
    case 3: // the face at y = 1
      // this face connects to another infinite element
      {
	AutoPtr<Elem> face(new InfQuad4);

	face->set_node(0) = get_node(2);
	face->set_node(1) = get_node(3);
	face->set_node(2) = get_node(7);
	face->set_node(3) = get_node(6);
	
	return face;
      }
    case 4: // the face at x = -1
      // this face connects to another infinite element
      {  
	AutoPtr<Elem> face(new InfQuad4);

	face->set_node(0) = get_node(3);
	face->set_node(1) = get_node(0);
	face->set_node(2) = get_node(4);
	face->set_node(3) = get_node(7);

	return face;
      }
    case 5: // the face at z = 1
      // disable this face, since this is supposed to lie at infinity
      {
        std::cerr << "No face 5 in case of infinite elements!" << std::endl;
        error();
	AutoPtr<Elem> ap(NULL);  return ap;
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



const std::vector<unsigned int> InfHex8::tecplot_connectivity(const unsigned int sc) const
{
  assert (_nodes != NULL);
  assert (sc < n_sub_elem());

  std::vector<unsigned int> conn(8);
  
  conn[0] = node(0)+1;
  conn[1] = node(1)+1;
  conn[2] = node(2)+1;
  conn[3] = node(3)+1;
  conn[4] = node(4)+1;
  conn[5] = node(5)+1;
  conn[6] = node(6)+1;
  conn[7] = node(7)+1;

  return conn;
}



void InfHex8::write_tecplot_connectivity(std::ostream &out) const
{
  assert (out);
  assert (_nodes != NULL);

  for (unsigned int sc=0; sc<n_sub_elem(); sc++)
    {
      std::vector<unsigned int> conn = tecplot_connectivity(sc);

      for (unsigned int i=0; i<conn.size(); i++)
	out << conn[i] << " ";

      out << std::endl;
    }
}



#ifdef ENABLE_AMR

const float  InfHex8::embedding_matrix[4][8][8] =
{
  // embedding matrix for child 0
  {
    //     0      1      2      3      4      5      6      7 th parent N.(ode)
    {    1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
    {    0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 1
    {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 2
    {    0.5,   0.0,   0.0,   0.5,   0.0,   0.0,   0.0,   0.0}, // 3
    {    0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0}, // 4
    {    0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0,   0.0}, // 5
    {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}, // 6
    {    0.0,   0.0,   0.0,   0.0,   0.5,   0.0,   0.0,   0.5}  // 7
  },

  // embedding matrix for child 1
  {
    //     0      1      2      3      4      5      6      7 th parent N.(ode)
    {    0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
    {    0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 1
    {    0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0}, // 2
    {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 3
    {    0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0,   0.0}, // 4
    {    0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0}, // 5
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0}, // 6
    {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}  // 7
  },

  // embedding matrix for child 2
  {
    //     0      1      2      3      4      5      6      7 th parent N.(ode)
    {    0.5,   0.0,   0.0,   0.5,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
    {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 1
    {    0.0,   0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0}, // 2
    {    0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0}, // 3
    {    0.0,   0.0,   0.0,   0.0,   0.5,   0.0,   0.0,   0.5}, // 4
    {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}, // 5
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5}, // 6
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0}  // 7
  },

  // embedding matrix for child 3
  {
    //     0      1      2      3      4      5      6      7 th parent N.(ode)
    {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
    {    0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0}, // 1
    {    0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 2
    {    0.0,   0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0}, // 3
    {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}, // 4
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0}, // 5
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0}, // 6
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5}  // 7
  }
};



const unsigned int InfHex8::side_children_matrix[6][5] =
{
  // note different storage scheme
  {4,   0, 1, 2, 3}, // 4 side-0 children
  {2,   0, 1,42,42}, // 2 side-1 children
  {2,   1, 3,42,42}, // 2 side-2 children
  {2,   2, 3,42,42}, // 2 side-3 children
  {2,   0, 2,42,42}, // 2 side-4 children
  {4,   0, 1, 2, 3}  // 4 side-5 children
};



void InfHex8::refine(Mesh& mesh)
{
  assert (refinement_flag() == Elem::REFINE);
  assert (active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[n_children()];

    for (unsigned int c=0; c<n_children(); c++)
      {
	_children[c] = new InfHex8(this);
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
	    p[c][nc].add_scaled (point(n), static_cast<Real>(embedding_matrix[c][nc][n]));
    
    
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
	    // the upper limit for sc is stored in the 0th column
	    for (unsigned int sc=1; sc<=side_children_matrix[s][0]; sc++)
	      mesh.boundary_info.add_side(child(side_children_matrix[s][sc]), s, id);

	}
  }


  // Un-set my refinement flag now
  set_refinement_flag() = Elem::DO_NOTHING;
}


#endif

#endif
