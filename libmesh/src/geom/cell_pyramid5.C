// $Id: cell_pyramid5.C,v 1.5 2003-01-24 17:24:43 jwpeterson Exp $

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
#include "cell_pyramid5.h"
#include "cell_tet4.h"
#include "face_tri3.h"
#include "face_quad4.h"
#include "edge_edge2.h"




// ------------------------------------------------------------
// Pyramid5 class member functions
AutoPtr<Elem> Pyramid5::build_side (const unsigned int i) const
{
  assert (i < n_sides());


  
  switch (i)
    {
    case 0:  // triangular face 1
      {
	Tri3* face = new Tri3;

	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(1);
	face->set_node(2) = get_node(4);
	
	AutoPtr<Elem> ap(face);  return ap;
      }
    case 1:  // triangular face 2
      {
	Tri3* face = new Tri3;

	face->set_node(0) = get_node(1);
	face->set_node(1) = get_node(2);
	face->set_node(2) = get_node(4);
	
	AutoPtr<Elem> ap(face);  return ap;
      }
    case 2:  // triangular face 3
      {
	Tri3* face = new Tri3;

	face->set_node(0) = get_node(2);
	face->set_node(1) = get_node(3);
	face->set_node(2) = get_node(4);
	
	AutoPtr<Elem> ap(face);  return ap;
      }
    case 3:  // triangular face 4
      {
	Tri3* face = new Tri3;

	face->set_node(0) = get_node(3);
	face->set_node(1) = get_node(0);
	face->set_node(2) = get_node(4);
	
	AutoPtr<Elem> ap(face);  return ap;
      }
    case 4:  // the quad face at z=0
      {
	Quad4* face = new Quad4;
	
	face->set_node(0) = get_node(0);
	face->set_node(1) = get_node(3);
	face->set_node(2) = get_node(2);
	face->set_node(3) = get_node(1);

	AutoPtr<Elem> ap(face);  return ap;
      }
    default:
      {
	error();
      }
    };

  // We'll never get here.
  error();

  AutoPtr<Elem> ap(NULL);  return ap;
};



const std::vector<unsigned int> Pyramid5::tecplot_connectivity(const unsigned int sc) const
{
  assert (_nodes != NULL);
  assert (sc < n_sub_elem());

  std::vector<unsigned int> conn(8);
  
  conn[0] = node(0)+1;
  conn[1] = node(1)+1;
  conn[2] = node(2)+1;
  conn[3] = node(3)+1;
  conn[4] = node(4)+1;
  conn[5] = node(4)+1;
  conn[6] = node(4)+1;
  conn[7] = node(4)+1;

  return conn;
};






void Pyramid5::vtk_connectivity(const unsigned int sc,
				std::vector<unsigned int> *conn) const
{
  assert (_nodes != NULL);
  assert (sc < n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(5);

  (*conn)[0] = node(0);
  (*conn)[1] = node(1);
  (*conn)[2] = node(2);
  (*conn)[3] = node(3);
  (*conn)[4] = node(4);

  return;
};



#ifdef ENABLE_AMR

void Pyramid5::refine(Mesh&)
{
  error();
  
  //   assert (refinement_flag() == Elem::REFINE);
  //   assert (active());
  //   assert (_children == NULL);

  //   // Create my children
  //   {
  //     _children = new Elem*[n_children()];

  //     for (unsigned int c=0; c<6; c++)
  //       {
  // 	_children[c] = new Pyramid5(this);
  //  	_children[c]->set_refinement_flag() = Elem::JUST_REFINED;
  //      }
  //     for (unsigned int c=6; c<n_children(); c++)
  //       {
  // 	_children[c] = new Tet4(this);
  //       	_children[c]->set_refinement_flag() = Elem::JUST_REFINED;
  //       };

  //   };

  //   /**
  //    * Build the nodes so we can look for them.	 
  //    */

  //   Edge2 edge;
  //   Quad4 face;

  
  
  //   edge.node(0) = node(0);
  //   edge.node(1) = node(1);
  
  //   const Point p0_1 = (mesh.point(node(0))*.5 +
  // 		      mesh.point(node(1))*.5  );
  
  //   const unsigned int n0_1 = mesh.mesh_refinement.add_point(p0_1, edge.key());
  

  //   edge.node(0) = node(1);
  //   edge.node(1) = node(2);
  
  //   const Point p1_2 = (mesh.point(node(1))*.5 +
  // 		      mesh.point(node(2))*.5  );
  
  //   const unsigned int n1_2 = mesh.mesh_refinement.add_point(p1_2, edge.key());
  

  //   edge.node(0) = node(2);
  //   edge.node(1) = node(3);
  
  //   const Point p2_3 = (mesh.point(node(2))*.5 +
  // 		      mesh.point(node(3))*.5  );
  
  //   const unsigned int n2_3 = mesh.mesh_refinement.add_point(p2_3, edge.key());
  

  //   edge.node(0) = node(0);
  //   edge.node(1) = node(3);
  
  //   const Point p0_3 = (mesh.point(node(0))*.5 +
  // 		      mesh.point(node(3))*.5  );
  
  //   const unsigned int n0_3 = mesh.mesh_refinement.add_point(p0_3, edge.key());
  

  //   edge.node(0) = node(0);
  //   edge.node(1) = node(4);
  
  //   const Point p0_4 = (mesh.point(node(0))*.5 +
  // 		      mesh.point(node(4))*.5  );
  
  //   const unsigned int n0_4 = mesh.mesh_refinement.add_point(p0_4, edge.key());
  

  //   edge.node(0) = node(1);
  //   edge.node(1) = node(4);
  
  //   const Point p1_4 = (mesh.point(node(1))*.5 +
  // 		      mesh.point(node(4))*.5  );
  
  //   const unsigned int n1_4 = mesh.mesh_refinement.add_point(p1_4, edge.key());
  

  //   edge.node(0) = node(2);
  //   edge.node(1) = node(4);
  
  //   const Point p2_4 = (mesh.point(node(2))*.5 +
  // 		      mesh.point(node(4))*.5  );

  //   const unsigned int n2_4 = mesh.mesh_refinement.add_point(p2_4, edge.key());
  

  //   edge.node(0) = node(3);
  //   edge.node(1) = node(4);
  
  //   const Point p3_4 = (mesh.point(node(3))*.5 +
  // 		      mesh.point(node(4))*.5  );
  
  //   const unsigned int n3_4 = mesh.mesh_refinement.add_point(p3_4, edge.key());
  

  //   face.node(0) = node(0);
  //   face.node(1) = node(3);
  //   face.node(2) = node(2);
  //   face.node(3) = node(1);
  
  //   const Point p0 = (mesh.point(node(0)) +
  // 		    mesh.point(node(3)) +
  // 		    mesh.point(node(2)) +
  // 		    mesh.point(node(1)))*.25;

  //   const unsigned int n0 = mesh.mesh_refinement.add_point(p0, face.key());
  
  

  //   // Tell the children about their new nodes
  //   // and add them to the mesh.
  //   {
  //     // Child 0 (Pyramid)
  //     {
  //       _children[0]->node(0) = node(0);
  //       _children[0]->node(1) = n0_1;
  //       _children[0]->node(2) = n0;
  //       _children[0]->node(3) = n0_3;
  //       _children[0]->node(4) = n0_4;
      
  //       mesh.add_elem(_children[0], mesh.mesh_refinement.new_element_number());
  //     };
    
  //     // Child 1 (Pyramid)
  //     {
  //       _children[1]->node(0) = n0_1;
  //       _children[1]->node(1) = node(1);
  //       _children[1]->node(2) = n1_2;
  //       _children[1]->node(3) = n0;
  //       _children[1]->node(4) = n1_4;

  //       mesh.add_elem(_children[1], mesh.mesh_refinement.new_element_number());
  //     };
    
  //     // Child 2 (Pyramid)
  //     {
  //       _children[2]->node(0) = n0_3;
  //       _children[2]->node(1) = n0;
  //       _children[2]->node(2) = n2_3;
  //       _children[2]->node(3) = node(3);
  //       _children[2]->node(4) = n3_4;

  //       mesh.add_elem(_children[2], mesh.mesh_refinement.new_element_number());
  //     };
    
  //     // Child 3 (Pyramid)
  //     {
  //       _children[3]->node(0) = n0;
  //       _children[3]->node(1) = n1_2;
  //       _children[3]->node(2) = node(2);
  //       _children[3]->node(3) = n2_3;
  //       _children[3]->node(4) = n2_4;
 
  //       mesh.add_elem(_children[3], mesh.mesh_refinement.new_element_number());
  //     };
    
  //     // Child 4 (Pyramid)
  //     {
  //       _children[4]->node(0) = n0_4;
  //       _children[4]->node(1) = n1_4;
  //       _children[4]->node(2) = n2_4;
  //       _children[4]->node(3) = n3_4;
  //       _children[4]->node(4) = node(4);

  //       mesh.add_elem(_children[4], mesh.mesh_refinement.new_element_number());
  //     };
    
  //     // Child 5 (Pyramid)
  //     {
  //       _children[5]->node(0) = n0_4;
  //       _children[5]->node(1) = n3_4;
  //       _children[5]->node(2) = n2_4;
  //       _children[5]->node(3) = n1_4;
  //       _children[5]->node(4) = n0;

  //       mesh.add_elem(_children[5], mesh.mesh_refinement.new_element_number());
  //     };
    
  //     // Child 6 (Tet)
  //     {
  //       _children[6]->node(0) = n0_1;
  //       _children[6]->node(1) = n0_4;
  //       _children[6]->node(2) = n1_4;
  //       _children[6]->node(3) = n0;

  //       mesh.add_elem(_children[6], mesh.mesh_refinement.new_element_number());
  //     };
    
  //     // Child 7 (Tet)
  //     {
  //       _children[7]->node(0) = n1_4;
  //       _children[7]->node(1) = n1_2;
  //       _children[7]->node(2) = n0;
  //       _children[7]->node(3) = n2_4;

  //       mesh.add_elem(_children[7], mesh.mesh_refinement.new_element_number());
  //     };
    
  //     // Child 8 (Tet)
  //     {
  //       _children[8]->node(0) = n0;
  //       _children[8]->node(1) = n2_4;
  //       _children[8]->node(2) = n2_3;
  //       _children[8]->node(3) = n3_4;

  //       mesh.add_elem(_children[8], mesh.mesh_refinement.new_element_number());
  //     };
    
  //     // Child 9 (Tet)
  //     {
  //       _children[9]->node(0) = n0_3;
  //       _children[9]->node(1) = n0;
  //       _children[9]->node(2) = n3_4;
  //       _children[9]->node(3) = n0_4;

  //       mesh.add_elem(_children[9], mesh.mesh_refinement.new_element_number());
  //     };    
  //   };


  
  //   // Possibly add boundary information
  //   {
  //     if (neighbor(0) == NULL)
  //       {
  // 	const short int id = mesh.boundary_info.boundary_id(this, 0);
	
  // 	if (id != mesh.boundary_info.invalid_id)
  // 	  {    
  // 	    mesh.boundary_info.add_side(_children[0], 0, id);
  // 	    mesh.boundary_info.add_side(_children[1], 0, id);
  // 	    mesh.boundary_info.add_side(_children[4], 0, id);
  // 	    mesh.boundary_info.add_side(_children[6], 0, id);
  // 	  }
  //       }

    

  //     if (neighbor(1) == NULL)
  //       {
  // 	const short int id = mesh.boundary_info.boundary_id(this, 1);
	
  // 	if (id != mesh.boundary_info.invalid_id)
  // 	  {    
  // 	    mesh.boundary_info.add_side(_children[1], 1, id);
  // 	    mesh.boundary_info.add_side(_children[3], 1, id);
  // 	    mesh.boundary_info.add_side(_children[4], 1, id);
  // 	    mesh.boundary_info.add_side(_children[7], 1, id);
  // 	  }
  //       }


    
  //     if (neighbor(2) == NULL)
  //       {
  // 	const short int id = mesh.boundary_info.boundary_id(this, 2);
	
  // 	if (id != mesh.boundary_info.invalid_id)
  // 	  {    
  // 	    mesh.boundary_info.add_side(_children[2], 2, id);
  // 	    mesh.boundary_info.add_side(_children[3], 2, id);
  // 	    mesh.boundary_info.add_side(_children[4], 2, id);
  // 	    mesh.boundary_info.add_side(_children[8], 2, id);
  // 	  }
  //       }


    
  //     if (neighbor(3) == NULL)
  //       {
  // 	const short int id = mesh.boundary_info.boundary_id(this, 3);
	
  // 	if (id != mesh.boundary_info.invalid_id)
  // 	  {    
  // 	    mesh.boundary_info.add_side(_children[0], 3, id);
  // 	    mesh.boundary_info.add_side(_children[2], 3, id);
  // 	    mesh.boundary_info.add_side(_children[4], 3, id);
  // 	    mesh.boundary_info.add_side(_children[9], 3, id);
  // 	  }
  //       }


    
  //     if (neighbor(4) == NULL)
  //       {
  // 	const short int id = mesh.boundary_info.boundary_id(this, 4);
	
  // 	if (id != mesh.boundary_info.invalid_id)
  // 	  {    
  // 	    mesh.boundary_info.add_side(_children[0], 4, id);
  // 	    mesh.boundary_info.add_side(_children[1], 4, id);
  // 	    mesh.boundary_info.add_side(_children[2], 4, id);
  // 	    mesh.boundary_info.add_side(_children[3], 4, id);
  // 	  }
  //       }
  //   };


  
  //   // Un-set my refinement flag now
  //   set_refinement_flag() = Elem::DO_NOTHING;
};



void Pyramid5::coarsen()
{
  assert (refinement_flag() == Elem::COARSEN);
  assert (!active());
  
  delete [] _children;

  _children = NULL;

  set_refinement_flag() = Elem::DO_NOTHING;
};


#endif
