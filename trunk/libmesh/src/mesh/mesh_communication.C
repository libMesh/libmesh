// $Id: mesh_communication.C,v 1.10 2004-11-08 00:11:05 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "libmesh_config.h"
#include "libmesh_common.h"
#include "mesh.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "mesh_communication.h"



// ------------------------------------------------------------
// MeshCommunication class members
void MeshCommunication::clear ()
{
}



void MeshCommunication::distribute (Mesh& mesh) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;
  
  this->distribute_mesh (mesh);
  this->distribute_bcs  (mesh, mesh.boundary_info);
}



void MeshCommunication::distribute_mesh (MeshBase& mesh) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;
  
#ifdef HAVE_MPI

  // Explicitly clear the mesh on all but processor 0.
  if (libMesh::processor_id() != 0)
    mesh.clear();
  
  // Get important sizes
  unsigned int n_nodes      = mesh.n_nodes();
  unsigned int n_elem       = mesh.n_elem();
  unsigned int total_weight = mesh.total_weight();

  // Broadcast the number of nodes
  MPI_Bcast (&n_nodes,      1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  // Send the number of elements
  MPI_Bcast (&n_elem,       1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  
  // Send the total_weight
  MPI_Bcast (&total_weight, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);  
  

  // First build up the pts vector which contains
  // the spatial locations of all the nodes      
  {
    std::vector<Real> pts;
	
    // If we are processor 0, we must populate this vector and
    // broadcast it to the other processors.
    if (libMesh::processor_id() == 0)
      {
	pts.reserve (3*n_nodes);
	
// 	node_iterator       it     (mesh.nodes_begin());
// 	const node_iterator it_end (mesh.nodes_end());

	MeshBase::node_iterator       it     = mesh.nodes_begin();
	const MeshBase::node_iterator it_end = mesh.nodes_end();

	for (; it != it_end; ++it)
	  {
	    assert (*it != NULL);
	    
	    const Point& p = **it;
	    
	    pts.push_back ( p(0) ); // x
	    pts.push_back ( p(1) ); // y
	    pts.push_back ( p(2) ); // z	  
	  }
      }
    else
      pts.resize (3*n_nodes);

    // Sanity check for all processors
    assert (pts.size() == (3*n_nodes));
    
    // Broadcast the pts vector
    MPI_Bcast (&pts[0], pts.size(), MPI_REAL, 0, MPI_COMM_WORLD);

    // Add the nodes we just received if we are not
    // processor 0.
    if (libMesh::processor_id() != 0)
      {
	assert (mesh.n_nodes() == 0);
	
	for (unsigned int i=0; i<pts.size(); i += 3)
	  mesh.add_point (Point(pts[i+0],
				pts[i+1],
				pts[i+2])
			  );
      }
    
    assert (mesh.n_nodes() == n_nodes);
  } // Done distributing the nodes

  
  // Now build up the elements vector which
  // contains the element types and connectivity
  {
    std::vector<unsigned int> conn;

    // If we are processor 0, we must populate this vector and
    // broadcast it to the other processors.
    if (libMesh::processor_id() == 0)
      {
	conn.reserve (n_elem + total_weight);
	
	// 	elem_iterator       it     (mesh.elements_begin());
	// 	const elem_iterator it_end (mesh.elements_end());

	MeshBase::element_iterator       it     = mesh.elements_begin();
	const MeshBase::element_iterator it_end = mesh.elements_end();

	for (; it != it_end; ++it)
	  {
	    const Elem* elem = *it;
	    
	    assert (elem != NULL);
	    
	    conn.push_back (static_cast<unsigned int>(elem->type()));
	    
	    for (unsigned int n=0; n<elem->n_nodes(); n++)
	      conn.push_back (elem->node(n));
	  }
      }
    else
      conn.resize (n_elem + total_weight);
    
    // Sanity check for all processors
    assert (conn.size() == (n_elem + total_weight));
    
    // Broadcast the element connectivity
    MPI_Bcast (&conn[0], conn.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    // Build the elements we just received if we are not
    // processor 0.
    if (libMesh::processor_id() != 0)
      {
	assert (mesh.n_elem() == 0);
	
	unsigned int cnt = 0;

	while (cnt < conn.size())
	  {
	    Elem* elem = 
	      mesh.add_elem (Elem::build(static_cast<ElemType>(conn[cnt++])).release());

	    for (unsigned int n=0; n<elem->n_nodes(); n++)
	      {
		assert (cnt < conn.size());
		
		elem->set_node(n) = mesh.node_ptr (conn[cnt++]);
	      }
	  }
      }
    
    assert (mesh.n_elem() == n_elem);
  } // Done distributing the elements

  // Print the information in the mesh for sanity.
  // mesh.print_info();
  
#else

  // no MPI but multiple processors? Huh??
  error();
  
#endif
}



void MeshCommunication::distribute_bcs (MeshBase& mesh,
					BoundaryInfo& boundary_info) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;
  
  
#ifdef HAVE_MPI

  // Explicitly clear the boundary conditions on all
  // but processor 0.
  if (libMesh::processor_id() != 0)
    boundary_info.clear();
  
  // Build up the list of elements with boundary conditions
  {
    std::vector<unsigned int>       el_id;
    std::vector<unsigned short int> side_id;
    std::vector<short int>          bc_id;

    if (libMesh::processor_id() == 0)
      boundary_info.build_side_list (el_id, side_id, bc_id);

    assert (el_id.size() == side_id.size());
    assert (el_id.size() == bc_id.size());
    
    unsigned int n_bcs = el_id.size();

    // Broadcast the number of bcs to expect from processor 0.
    MPI_Bcast (&n_bcs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    
    // Allocate space.
    el_id.resize   (n_bcs);
    side_id.resize (n_bcs);
    bc_id.resize   (n_bcs);

    // Broadcast the element identities
    MPI_Bcast (&el_id[0],   n_bcs, MPI_UNSIGNED,       0, MPI_COMM_WORLD);

    // Broadcast the side ids for those elements
    MPI_Bcast (&side_id[0], n_bcs, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

    // Broadcast the bc ids for each side
    MPI_Bcast (&bc_id[0],   n_bcs, MPI_SHORT,          0, MPI_COMM_WORLD);

    // Build the boundary_info structure if we aren't processor 0
    if (libMesh::processor_id() != 0)
      for (unsigned int e=0; e<n_bcs; e++)
	{
	  assert (el_id[e] < mesh.n_elem());
	    
	  const Elem* elem = mesh.elem(el_id[e]);

	  assert (elem != NULL);
	  assert (side_id[e] < elem->n_sides());
	    
	  boundary_info.add_side (elem, side_id[e], bc_id[e]);
	}
  }



  // Build up the list of nodes with boundary conditions
  {
    std::vector<unsigned int> node_id;
    std::vector<short int>    bc_id;

    if (libMesh::processor_id() == 0)
      boundary_info.build_node_list (node_id, bc_id);

    assert (node_id.size() == bc_id.size());
    
    unsigned int n_bcs = node_id.size();

    // Broadcast the number of bcs to expect from processor 0.
    MPI_Bcast (&n_bcs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    // Allocate space.
    node_id.resize (n_bcs);
    bc_id.resize   (n_bcs);

    // Broadcast the node ids
    MPI_Bcast (&node_id[0], n_bcs, MPI_UNSIGNED,       0, MPI_COMM_WORLD);

    // Broadcast the bc ids for each side
    MPI_Bcast (&bc_id[0],   n_bcs, MPI_SHORT,          0, MPI_COMM_WORLD);

    // Build the boundary_info structure if we aren't processor 0
    if (libMesh::processor_id() != 0)
      for (unsigned int n=0; n<n_bcs; n++)
	{
	  assert (node_id[n] < mesh.n_nodes());
	    
	  const Node* node = mesh.node_ptr (node_id[n]);

	  assert (node != NULL);
	    
	  boundary_info.add_node (node, bc_id[n]);
	}
  }
    
      

#else

  // no MPI but multiple processors? Huh??
  error();

#endif  
}
