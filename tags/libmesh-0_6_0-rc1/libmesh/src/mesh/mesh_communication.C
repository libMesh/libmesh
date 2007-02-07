// $Id: mesh_communication.C,v 1.26 2006-08-28 16:08:22 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "libmesh_logging.h"
#include "mesh_base.h"
#include "mesh_tools.h"
#include "boundary_info.h"
#include "mesh_communication.h"
#include "elem.h"
#include "sphere.h"



// ------------------------------------------------------------
// MeshCommunication class members
void MeshCommunication::clear ()
{
  _neighboring_processors.clear();
}


void MeshCommunication::find_neighboring_processors (const MeshBase& mesh)
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;
  
#ifdef HAVE_MPI
  
  _neighboring_processors.clear();

  // Get the bounding sphere for the local processor
  Sphere bounding_sphere =
    MeshTools::processor_bounding_sphere (mesh, libMesh::processor_id());

  // Just to be sure, increase its radius by 10%.  Sure would suck to
  // miss a neighboring processor!
  bounding_sphere.radius() *= 1.1;

  // Collect the bounding spheres from all processors, test for intersection
  {
    std::vector<float>
      send (4,                         0),
      recv (4*libMesh::n_processors(), 0);

    send[0] = bounding_sphere.center()(0);
    send[1] = bounding_sphere.center()(1);
    send[2] = bounding_sphere.center()(2);
    send[3] = bounding_sphere.radius();

    MPI_Allgather (&send[0], send.size(), MPI_FLOAT,
		   &recv[0], send.size(), MPI_FLOAT,
		   libMesh::COMM_WORLD);


    for (unsigned int proc=0; proc<libMesh::n_processors(); proc++)
      {
	const Point center (recv[4*proc+0],
			    recv[4*proc+1],
			    recv[4*proc+2]);
	
	const Real radius = recv[4*proc+3];

	const Sphere proc_sphere (center, radius);

	if (bounding_sphere.intersects(proc_sphere))
	  _neighboring_processors.push_back(proc);
      }

    // Print out the _neighboring_processors list
    std::cout << "Processor " << libMesh::processor_id()
	      << " intersects:" << std::endl;
    for (unsigned int p=0; p<_neighboring_processors.size(); p++)
      std::cout << " " << _neighboring_processors[p] << std::endl;
  }
  
#endif
}


void MeshCommunication::broadcast (MeshBase& mesh) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;
  
  this->broadcast_mesh (mesh);
  this->broadcast_bcs  (mesh, *(mesh.boundary_info));
}



#ifdef HAVE_MPI
void MeshCommunication::broadcast_mesh (MeshBase& mesh) const
#else // avoid spurious gcc warnings
void MeshCommunication::broadcast_mesh (MeshBase&) const
#endif
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;
  
#ifdef HAVE_MPI

  START_LOG("broadcast_mesh()","MeshCommunication");

  // Explicitly clear the mesh on all but processor 0.
  if (libMesh::processor_id() != 0)
    mesh.clear();
  
  // Get important sizes
  unsigned int n_nodes      = mesh.n_nodes();
  unsigned int n_elem       = mesh.n_elem();
  unsigned int n_levels     = MeshTools::n_levels(mesh);

  // For adaptive meshes, need to store parent and self IDs as well
  unsigned int total_weight = MeshTools::total_weight(mesh) + 2*n_elem;

  // Broadcast the sizes
  {
    std::vector<unsigned int> buf (3);
    
    buf[0] = n_nodes;
    buf[1] = n_elem;
    buf[2] = total_weight;
    
    // Broadcast
    MPI_Bcast (&buf[0], buf.size(), MPI_UNSIGNED, 0, libMesh::COMM_WORLD);

    if (libMesh::processor_id() != 0)
      {
	n_nodes      = buf[0];
	n_elem       = buf[1];
	total_weight = buf[2];
      }	
  }  

  // First build up the pts vector which contains
  // the spatial locations of all the nodes      
  {
    std::vector<Real> pts;
	
    // If we are processor 0, we must populate this vector and
    // broadcast it to the other processors.
    if (libMesh::processor_id() == 0)
      {
	pts.reserve (3*n_nodes);

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
    MPI_Bcast (&pts[0], pts.size(), MPI_REAL, 0, libMesh::COMM_WORLD);

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
    // The conn array contains the information needed to construct each element.
    // Pack all this information into one communication to avoid two latency hits
    // For each element it is of the form
    // [ etype subdomain_id self_ID parent_ID node_0 node_1 ... node_n]
    // We cannot use unsigned int because parent_ID can be negative
    std::vector<int> conn;

    // If we are processor 0, we must populate this vector and
    // broadcast it to the other processors.
    if (libMesh::processor_id() == 0)
      {
	conn.reserve (2*n_elem + total_weight);
	
	// We start from level 0. This is a bit simpler than in xdr_io.C
	// because we do not have to worry about economizing by group elements
	// of the same type. Element type is simply specified as an
	// entry in the connectivity vector, "conn".
	// By filling conn in order of levels, parents should exist before children
	// are built when we reconstruct the elements on the other processors.
	
	for(unsigned int level=0; level<=n_levels; ++level)
	  {
	    MeshBase::element_iterator it = mesh.level_elements_begin(level);
	    const MeshBase::element_iterator it_end = mesh.level_elements_end(level);
	    
	    for (; it != it_end; ++it)
	      {
		const Elem* elem = *it;
		
		assert (elem != NULL);
		
		conn.push_back (static_cast<int>(elem->type()));
		conn.push_back (static_cast<int>(elem->subdomain_id()));
		conn.push_back (elem->id());
		
		// use parent_ID of -1 to indicate a level 0 element
		if (level==0)
		  conn.push_back(-1);
		else
		  conn.push_back(elem->parent()->id());
		
		for (unsigned int n=0; n<elem->n_nodes(); n++)
		  conn.push_back (elem->node(n));		
	      }
	  }
      }
    else
      conn.resize (2*n_elem + total_weight);
    
    // Sanity check for all processors
    assert (conn.size() == (2*n_elem + total_weight));
    
    // Broadcast the element connectivity
    MPI_Bcast (&conn[0], conn.size(), MPI_INT, 0, libMesh::COMM_WORLD);

    // Build the elements we just received if we are not
    // processor 0.
    if (libMesh::processor_id() != 0)
      {
	assert (mesh.n_elem() == 0);
	
	unsigned int cnt = 0;

        // This map keeps track of elements we've previously added to the mesh 
        // to avoid O(n) lookup times for parent pointers.
        std::map<unsigned int, Elem*> parents;

	while (cnt < conn.size())
	  {
	    // Declare the element that we will add
            Elem* elem = NULL;
            
            const ElemType elem_type = static_cast<ElemType>(conn[cnt++]);
              
	    // Get the subdomain id

            // get ID info
	    const int subdomain_id = conn[cnt++];
            const int self_ID      = conn[cnt++];
            const int parent_ID    = conn[cnt++];
	    
	    
            if (parent_ID != -1) // Do a log(n) search for the parent
	      {
		Elem* my_parent;
                
                // Search for parent in the parents map (log(n))
                std::map<unsigned int, Elem*>::iterator it = parents.find(parent_ID);
                
                // If the parent was not previously added, we cannot continue.
                if (it == parents.end())
                {
                  std::cerr << "Parent element with ID " << parent_ID 
                            << " not found." << std::endl; 
                  error();
                }

                // Set the my_parent pointer
                my_parent = (*it).second;
		
		my_parent->set_refinement_flag(Elem::INACTIVE);
		
		elem = mesh.add_elem(Elem::build(elem_type,my_parent).release());
		elem->set_refinement_flag(Elem::JUST_REFINED); 
		my_parent->add_child(elem);
		assert (my_parent->type() == elem->type());
	      }

            else // level 0 element has no parent
	      {
		// should be able to just use the integer elem_type
		elem = mesh.add_elem (Elem::build(elem_type).release());
	      }

	    // Assign the IDs
            elem->subdomain_id() = subdomain_id;
            elem->set_id() = self_ID;

            // Add elem to the map of parents, since it may have
            // children to be added later
            parents[self_ID] = elem;
	    
	    // Assign the connectivity
	    for (unsigned int n=0; n<elem->n_nodes(); n++)
	      {
		assert (cnt < conn.size());
		
		elem->set_node(n) = mesh.node_ptr (conn[cnt++]);
	      }
	  } // end while cnt < conn.size

	// All the elements at each level have been added, and their node pointers
	// have been set.  Now compute the node keys to put the mesh into a state consistent
	// with the state after being constructed through normal refinements. 
	// The new nodes were added through Mesh::add_point(), which
	// allocates brand new memory for the point, and does not copy
	// over any key.
	MeshBase::element_iterator it = mesh.elements_begin();
	const MeshBase::element_iterator end = mesh.elements_end();
	for (; it!=end; ++it)
	  (*it)->compute_children_node_keys();
	
      } // end if iam != cpu 0
    
    
    assert (mesh.n_elem() == n_elem);
  } // Done distributing the elements


  STOP_LOG("broadcast_mesh()","MeshCommunication");
  
#else

  // no MPI but multiple processors? Huh??
  error();
  
#endif
}



#ifdef HAVE_MPI
void MeshCommunication::broadcast_bcs (MeshBase& mesh,
					BoundaryInfo& boundary_info) const
#else // avoid spurious gcc warnings
void MeshCommunication::broadcast_bcs (MeshBase&,
					BoundaryInfo&) const
#endif
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;
  
  
#ifdef HAVE_MPI

  START_LOG("broadcast_bcs()","MeshCommunication");

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
    MPI_Bcast (&n_bcs, 1, MPI_UNSIGNED, 0, libMesh::COMM_WORLD);

    // Only continue if we have element BCs
    if (n_bcs > 0)
      {
	// Allocate space. On CPU 0, these vectors should already have size n_bcs.
	el_id.resize   (n_bcs);
	side_id.resize (n_bcs);
	bc_id.resize   (n_bcs);
	
	// Broadcast the element identities
	MPI_Bcast (&el_id[0],   n_bcs, MPI_UNSIGNED,       0, libMesh::COMM_WORLD);

	// Broadcast the side ids for those elements
	MPI_Bcast (&side_id[0], n_bcs, MPI_UNSIGNED_SHORT, 0, libMesh::COMM_WORLD);

	// Broadcast the bc ids for each side
	MPI_Bcast (&bc_id[0],   n_bcs, MPI_SHORT,          0, libMesh::COMM_WORLD);

	// Build the boundary_info structure if we aren't processor 0
	if (libMesh::processor_id() != 0)
	  for (unsigned int e=0; e<n_bcs; e++)
	    {
	      assert (el_id[e] < mesh.n_elem());
	      
	      const Elem* elem = mesh.elem(el_id[e]);

	      assert (elem != NULL);

	      // sanity: be sure that the element returned by mesh.elem() really has id()==el_id[e]
	      assert(elem->id() == el_id[e]);

	      assert (side_id[e] < elem->n_sides());
	    
	      boundary_info.add_side (elem, side_id[e], bc_id[e]);
	    }
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
    MPI_Bcast (&n_bcs, 1, MPI_UNSIGNED, 0, libMesh::COMM_WORLD);

    // Only continue if we have nodal BCs
    if (n_bcs > 0)
      {      
	// Allocate space, again on CPU 0 this should be a no-op.
	node_id.resize (n_bcs);
	bc_id.resize   (n_bcs);
	
	// Broadcast the node ids
	MPI_Bcast (&node_id[0], n_bcs, MPI_UNSIGNED, 0, libMesh::COMM_WORLD);
	
	// Broadcast the bc ids for each side
	MPI_Bcast (&bc_id[0],   n_bcs, MPI_SHORT,    0, libMesh::COMM_WORLD);

	// Build the boundary_info structure if we aren't processor 0
	if (libMesh::processor_id() != 0)
	  for (unsigned int n=0; n<n_bcs; n++)
	    {
	      assert (node_id[n] < mesh.n_nodes());
	      
	      const Node* node = mesh.node_ptr (node_id[n]);

	      assert (node != NULL);
	      
	      // sanity: be sure that the node returned by mesh.node_ptr() really has id()==node_id[n]
	      assert(node->id() == node_id[n]);
	    
	      boundary_info.add_node (node, bc_id[n]);
	    }
      }
  }

  STOP_LOG("broadcast_bcs()","MeshCommunication");
          
#else

  // no MPI but multiple processors? Huh??
  error();

#endif  
}
