// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include "parallel_mesh.h"
#include "mesh_tools.h"
#include "boundary_info.h"
#include "mesh_communication.h"
#include "parallel.h"
#include "elem.h"



// ------------------------------------------------------------
// MeshCommunication class members
void MeshCommunication::clear ()
{
  //  _neighboring_processors.clear();
}


// #ifdef HAVE_MPI
// void MeshCommunication::find_neighboring_processors (const MeshBase& mesh)
// {
//   // Don't need to do anything if there is
//   // only one processor.
//   if (libMesh::n_processors() == 1)
//     return;
  
//   _neighboring_processors.clear();

//   // Get the bounding sphere for the local processor
//   Sphere bounding_sphere =
//     MeshTools::processor_bounding_sphere (mesh, libMesh::processor_id());

//   // Just to be sure, increase its radius by 10%.  Sure would suck to
//   // miss a neighboring processor!
//   bounding_sphere.radius() *= 1.1;

//   // Collect the bounding spheres from all processors, test for intersection
//   {
//     std::vector<float>
//       send (4,                         0),
//       recv (4*libMesh::n_processors(), 0);

//     send[0] = bounding_sphere.center()(0);
//     send[1] = bounding_sphere.center()(1);
//     send[2] = bounding_sphere.center()(2);
//     send[3] = bounding_sphere.radius();

//     MPI_Allgather (&send[0], send.size(), MPI_FLOAT,
// 		   &recv[0], send.size(), MPI_FLOAT,
// 		   libMesh::COMM_WORLD);


//     for (unsigned int proc=0; proc<libMesh::n_processors(); proc++)
//       {
// 	const Point center (recv[4*proc+0],
// 			    recv[4*proc+1],
// 			    recv[4*proc+2]);
	
// 	const Real radius = recv[4*proc+3];

// 	const Sphere proc_sphere (center, radius);

// 	if (bounding_sphere.intersects(proc_sphere))
// 	  _neighboring_processors.push_back(proc);
//       }

//     // Print out the _neighboring_processors list
//     std::cout << "Processor " << libMesh::processor_id()
// 	      << " intersects:" << std::endl;
//     for (unsigned int p=0; p<_neighboring_processors.size(); p++)
//       std::cout << " " << _neighboring_processors[p] << std::endl;
//   }
// }
// #else
// void MeshCommunication::find_neighboring_processors (const MeshBase&)
// {
// }
// #endif



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
  unsigned int total_weight = MeshTools::total_weight(mesh);

  // Broadcast the sizes
  {
    std::vector<unsigned int> buf (3);
    
    buf[0] = n_nodes;
    buf[1] = n_elem;
    buf[2] = total_weight;
    
    // Broadcast
    Parallel::broadcast (buf);

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
    Parallel::broadcast (pts);

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
    // [ level etype subdomain_id self_ID parent_ID node_0 node_1 ... node_n]
    // We cannot use unsigned int because parent_ID can be negative
    std::vector<int> conn;

    // If we are processor 0, we must populate this vector and
    // broadcast it to the other processors.
    if (libMesh::processor_id() == 0)
      {
	conn.reserve (5*n_elem + total_weight);
	
	// We start from level 0. This is a bit simpler than in xdr_io.C
	// because we do not have to worry about economizing by group elements
	// of the same type. Element type is simply specified as an
	// entry in the connectivity vector, "conn".
	// By filling conn in order of levels, parents should exist before children
	// are built when we reconstruct the elements on the other processors.
	
	for (unsigned int level=0; level<=n_levels; ++level)
	  {
	    MeshBase::element_iterator it = mesh.level_elements_begin(level);
	    const MeshBase::element_iterator it_end = mesh.level_elements_end(level);
	    
	    for (; it != it_end; ++it)
	      {
		const Elem* elem = *it;
		pack_element (conn, elem);
	      }
	  }
      }
    else
      conn.resize (5*n_elem + total_weight);
    
    // Sanity check for all processors
    assert (conn.size() == (5*n_elem + total_weight));
    
    // Broadcast the element connectivity
    Parallel::broadcast (conn);

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

	    // Unpack the element header
	    const int level          = conn[cnt++];
            const ElemType elem_type = static_cast<ElemType>(conn[cnt++]);
	    const int subdomain_ID   = conn[cnt++];
            const int self_ID        = conn[cnt++];
            const int parent_ID      = conn[cnt++];
	    
#ifdef ENABLE_AMR

            if (parent_ID != -1) // Do a log(n) search for the parent
	      {
		Elem* my_parent = parents.count(parent_ID) ? parents[parent_ID] : NULL;
                
                // If the parent was not previously added, we cannot continue.
                if (my_parent == NULL)
		  {
		    std::cerr << "Parent element with ID " << parent_ID 
			      << " not found." << std::endl; 
		    error();
		  }
		
		my_parent->set_refinement_flag(Elem::INACTIVE);
		
		elem = Elem::build(elem_type,my_parent).release();
		elem->set_refinement_flag(Elem::JUST_REFINED); 
		my_parent->add_child(elem);
		assert (my_parent->type() == elem->type());
	      }
	    
            else // level 0 element has no parent
#endif 
	      {
		assert (level == 0);
		
		// should be able to just use the integer elem_type
		elem = Elem::build(elem_type).release();
	      }

	    // Assign the IDs
	    assert (elem->level() == static_cast<unsigned int>(level));
            elem->subdomain_id() = subdomain_ID;
            elem->set_id() = self_ID;
	    
            // Add elem to the map of parents, since it may have
            // children to be added later
            parents.insert(std::make_pair(self_ID,elem));
	    
	    // Assign the connectivity
	    for (unsigned int n=0; n<elem->n_nodes(); n++)
	      {
		assert (cnt < conn.size());
		
		elem->set_node(n) = mesh.node_ptr (conn[cnt++]);
	      }
	  } // end while cnt < conn.size

        // Iterate in ascending elem ID order
        for (std::map<unsigned int, Elem *>::iterator i =
             parents.begin();
             i != parents.end(); ++i)
          {
            Elem *elem = i->second;
            if (elem)
              mesh.add_elem(elem);
            else
              // We can probably handle this, but we don't expect it
              error();
          }

#ifdef ENABLE_AMR
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

#endif // #ifdef ENABLE_AMR
	
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
void MeshCommunication::broadcast_bcs (const MeshBase& mesh,
				       BoundaryInfo& boundary_info) const
#else // avoid spurious gcc warnings
void MeshCommunication::broadcast_bcs (const MeshBase&,
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
    Parallel::broadcast (n_bcs);

    // Only continue if we have element BCs
    if (n_bcs > 0)
      {
	// Allocate space. On CPU 0, these vectors should already have size n_bcs.
	el_id.resize   (n_bcs);
	side_id.resize (n_bcs);
	bc_id.resize   (n_bcs);
	
	// Broadcast the element identities
	Parallel::broadcast (el_id);

	// Broadcast the side ids for those elements
	Parallel::broadcast (side_id);

	// Broadcast the bc ids for each side
	Parallel::broadcast (bc_id);

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
    Parallel::broadcast (n_bcs);

    // Only continue if we have nodal BCs
    if (n_bcs > 0)
      {      
	// Allocate space, again on CPU 0 this should be a no-op.
	node_id.resize (n_bcs);
	bc_id.resize   (n_bcs);
	
	// Broadcast the node ids
	Parallel::broadcast (node_id);
	
	// Broadcast the bc ids for each side
	Parallel::broadcast (bc_id);

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



void MeshCommunication::allgather (ParallelMesh& mesh) const
{
  this->allgather_mesh (mesh);
  this->allgather_bcs  (mesh, *(mesh.boundary_info));
}

#ifndef HAVE_MPI
  
void MeshCommunication::allgather_mesh (ParallelMesh&) const
{
  // NO MPI == one processor, no need for this method
  return;
}
  
void MeshCommunication::allgather_bcs (const ParallelMesh&,
				       BoundaryInfo&) const
{
  // NO MPI == one processor, no need for this method
  return;
}

#else

void MeshCommunication::allgather_mesh (ParallelMesh& mesh) const
{
  // Check for quick return
  if (libMesh::n_processors() == 1)
    return;

  START_LOG ("allgather_mesh()","MeshCommunication");
  
  // Gather the number of nodes and elements on each processor.
  std::vector<unsigned int>
    n_nodes(libMesh::n_processors()), n_elem(libMesh::n_processors());
  
  {
    std::vector<unsigned int> n_objects(2);
    n_objects[0] = mesh.n_local_nodes();
    n_objects[1] = mesh.n_local_elem();

    Parallel::allgather(n_objects);
    
    for (unsigned int p=0, idx=0; p<libMesh::n_processors(); p++)
      {
	n_nodes[p] = n_objects[idx++];
	n_elem[p]  = n_objects[idx++];	
      }

    assert (mesh.n_local_nodes() == n_nodes[libMesh::processor_id()]);
    assert (mesh.n_local_elem()  ==  n_elem[libMesh::processor_id()]);
  }

  std::vector<unsigned int>
    node_offsets(libMesh::n_processors(), 0),
    elem_offsets(libMesh::n_processors(), 0);
  
  // Compute the global sizes to cross-check the results of the
  // operations that follow.
  unsigned int
    global_n_nodes = n_nodes[0], 
    global_n_elem  = n_elem[0];

  for (unsigned int p=1; p<libMesh::n_processors(); p++)
    {
      node_offsets[p] = node_offsets[p-1] + n_nodes[p-1];
      elem_offsets[p] = elem_offsets[p-1] +  n_elem[p-1];

      global_n_nodes += n_nodes[p];
      global_n_elem  += n_elem[p];
    }

  
  
  //-------------------------------------------------
  // Gather the nodal coordinates from each processor.
  {
    std::vector<Real> xyz; xyz.reserve(3*n_nodes[libMesh::processor_id()]);
    
    ParallelMesh::node_iterator       it  = mesh.local_nodes_begin();
    const ParallelMesh::node_iterator end = mesh.local_nodes_end();

    for (; it != end; ++it)
      {
	assert (*it != NULL);
	    
	const Point &p = **it;

	xyz.push_back(p(0));
	xyz.push_back(p(1));
	xyz.push_back(p(2));
      }

    assert (xyz.size() == 3*n_nodes[libMesh::processor_id()]);

    // Get values from other processors
    Parallel::allgather (xyz);

    // And add them to our mesh.
    for (unsigned int p=0; p<libMesh::n_processors(); p++)
      if (p == libMesh::processor_id()) continue; // We've already got our
                                                  // own local nodes!
      else
	{
	  const unsigned int
	    first_global_idx = node_offsets[p],
	    last_global_idx  = first_global_idx + n_nodes[p];

	  // Extract the coordinates for each node belonging to processor p
	  // and add it to our mesh.
	  for (unsigned int global_idx = first_global_idx; global_idx<last_global_idx; global_idx++)
	    {	      
	      Node *node = Node::build(xyz[3*global_idx + 0],
				       xyz[3*global_idx + 1],
				       xyz[3*global_idx + 2],
				       global_idx).release();
	      
	      assert (node != NULL);
	      assert (node->id() == global_idx);
	      
	      node->processor_id() = p;

	      mesh.insert_node(node);
	    }
	}
    
    // Check the result
    assert (global_n_nodes == mesh.n_nodes());
  }
  
  
  //----------------------------------------------------
  // Gather the element connectivity from each processor.
  {
    // Get the sum of elem->n_nodes() for all local elements.  This
    // will allow for efficient preallocation.
    const unsigned int
      local_weight   = MeshTools::weight(mesh),
      local_n_levels = MeshTools::n_levels(mesh); // Strictly speaking, this looks at all elements,
                                                  // not just the local ones.  That is OK, though. 

    unsigned int global_n_levels = local_n_levels;
    Parallel::max (global_n_levels);
    
    // The conn array contains the information needed to construct each element.
    std::vector<int> conn; conn.reserve (5*n_elem[libMesh::processor_id()] + local_weight);
						
    for (unsigned int level=0; level<=local_n_levels; level++)
      {
	// TODO:[BSK] implement local_level_elements iterators
	ParallelMesh::element_iterator        it  = mesh.level_elements_begin(level);
	const ParallelMesh::element_iterator  end = mesh.level_elements_end(level);

	for (; it != end; ++it)
	  {
	    const Elem* elem = *it;

	    assert (elem != NULL);
	    assert (elem->level() == level);
	    
	    // Only local elements!
	    if (elem->processor_id() != libMesh::processor_id()) continue;

	    pack_element (conn, elem);	    
	  }
      } // ...that was easy.

    assert (conn.size() == 5*n_elem[libMesh::processor_id()] + local_weight);

    // Get the size of the connectivity array on each processor
    std::vector<unsigned int>
      conn_size   (libMesh::n_processors(), 0),
      conn_offset (libMesh::n_processors(), 0);
    
    Parallel::allgather (static_cast<unsigned int>(conn.size()), conn_size);

    for (unsigned int p=1; p<libMesh::n_processors(); p++)
      conn_offset[p] = conn_offset[p-1] + conn_size[p-1];    
    
    // Get the element connectivity from all the other processors
    Parallel::allgather (conn);


        
    // ...and add them to our mesh.
    // This is a little tricky.  We need to insure that parents are added before children. 
    // So we need to add elements level-wise to handle, for example, the case where a child on
    // processor [0] has a parent on processor [1].  But we also need to add the elements 
    // processor-wise so that we can set the processor_id() properly.  
    // So, loop on levels/processors
    for (unsigned int level=0; level<=global_n_levels; level++)
      for (unsigned int p=0; p<libMesh::n_processors(); p++)
	if (p == libMesh::processor_id()) continue; // We've already got our
                                                    // own local elements!
	else
	  {
	    unsigned int cnt = conn_offset[p]; // counter into the conn[] array.
	    
	    const unsigned int
	      first_global_idx = elem_offsets[p],
	      last_global_idx  = first_global_idx + n_elem[p];

	    // Process each element for processor p.
	    // Note this must work in the case when conn_size[p] == 0.
	    while (cnt < (conn_offset[p] + conn_size[p]))
	      {
		// Unpack the element header
		const unsigned int elem_level = conn[cnt++];
		const ElemType elem_type      = static_cast<ElemType>(conn[cnt++]);
		const int subdomain_ID        = conn[cnt++];
		const unsigned int self_ID    = conn[cnt++];
		const int parent_ID           = conn[cnt++];
	       
		// We require contiguous numbering on each processor
		// for elements.
		assert (self_ID >= first_global_idx);
		assert (self_ID  < last_global_idx);
		assert ((level == 0) || (parent_ID != -1));
  
		// Ignore elements not matching the current level.  We
		// have to do this in a somewhat expensive fashion since
		// there is no good way to determine the number of nodes
		// in an element type without constructing one.
		if (elem_level > level) // build an element of elem_type and skip elem->n_nodes()
		  {                     // entries in the conn array.
		    AutoPtr<Elem> elem = Elem::build (elem_type);

		    cnt += elem->n_nodes();
		  }

		else if (elem_level < level) // we should already have this element, so there
		  {                          // is no need to construct a dummy element of this type
		    const Elem* elem = mesh.elem(self_ID);

		    assert (elem->subdomain_id() == subdomain_ID);
		    assert (elem->id()           == self_ID);

		    cnt += elem->n_nodes();
		  }

		// Those are the easy cases...
		else // elem_level == level
		  {
		    // Declare the element we will add
		    Elem* elem = NULL;

		    // Maybe find its parent
		    if (level > 0)
		      {
			Elem* my_parent = mesh.elem(parent_ID);

			// If the parent was not previously added, we cannot continue.
			if (my_parent == NULL)
			  {
			    std::cerr << "Parent element with ID " << parent_ID 
				      << " not found." << std::endl; 
			    error();
			  }
		
			elem = Elem::build(elem_type,my_parent).release();
#ifdef ENABLE_AMR
			my_parent->set_refinement_flag(Elem::INACTIVE);
			elem->set_refinement_flag(Elem::JUST_REFINED);
			my_parent->add_child(elem);
#endif // ENABLE_AMR
			assert (my_parent->type() == elem->type());
		      }
		    else
		      elem = Elem::build(elem_type).release();
		      
		    // Assign the IDs
		    assert (elem->level() == static_cast<unsigned int>(level));
		    elem->subdomain_id() = subdomain_ID;
		    elem->processor_id() = p;
		    elem->set_id()       = self_ID;
	    
		    // Assign the connectivity
		    for (unsigned int n=0; n<elem->n_nodes(); n++)
		      {
			assert (cnt < conn.size());
			
			elem->set_node(n) = mesh.node_ptr (conn[cnt++]);
		      }

		    // Good to go.  Add to the mesh.
		    mesh.insert_elem(elem);
		    
		  } // end elem_level == level
	      }
	    
	  }   
#ifdef ENABLE_AMR
    // All the elements at each level have been added, and their node pointers
    // have been set.  Now compute the node keys to put the mesh into a state consistent
    // with the state after being constructed through normal refinements. 
    // The new nodes were added through Mesh::add_point(), which
    // allocates brand new memory for the point, and does not copy
    // over any key.
    ParallelMesh::element_iterator it = mesh.elements_begin();
    const ParallelMesh::element_iterator end = mesh.elements_end();

    for (; it!=end; ++it)
      (*it)->compute_children_node_keys();

#endif // #ifdef ENABLE_AMR

    // Check the result
    assert (global_n_elem == mesh.n_elem());
  }

  // All done!
  STOP_LOG ("allgather_mesh()","MeshCommunication");
}



void MeshCommunication::allgather_bcs (const ParallelMesh& mesh,
				       BoundaryInfo& boundary_info) const
{
  // Check for quick return
  if (libMesh::n_processors() == 1)
    return;

  START_LOG ("allgather_bcs()","MeshCommunication");

  std::vector<int>
    xfer_elem_bcs,
    xfer_node_bcs;
  
  
  // Get the element boundary conditions
  {    
    std::vector<unsigned int>       el_id;
    std::vector<unsigned short int> side_id;
    std::vector<short int>          bc_id;
    
    boundary_info.build_side_list (el_id, side_id, bc_id);

    assert (el_id.size() == side_id.size());
    assert (el_id.size() == bc_id.size());
    
    const unsigned int n_bcs = el_id.size();

    // reserve an upper bound for the number of BCs
    xfer_elem_bcs.reserve(3*n_bcs);

    // populate the xfer_elem_bcs list with *local* elements only.
    for (unsigned int bc=0; bc<n_bcs; bc++)
      {
	const Elem* elem = mesh.elem(el_id[bc]);
	
	// sanity: be sure that the element returned by mesh.elem() really has id()==el_id[e]
	assert(elem != NULL);
	assert(elem->id() == el_id[bc]);
	assert(elem->level() == 0);
	assert(side_id[bc] < elem->n_sides());

	if (elem->processor_id() == libMesh::processor_id())
	  {
	    xfer_elem_bcs.push_back(el_id[bc]);
	    xfer_elem_bcs.push_back(side_id[bc]);
	    xfer_elem_bcs.push_back(bc_id[bc]);
	  }
      }
  
  } // done with element boundary conditions


  // Get the nodal boundary conditions
  {
    std::vector<unsigned int> node_id;
    std::vector<short int>    bc_id;
    
    boundary_info.build_node_list (node_id, bc_id);

    assert (node_id.size() == bc_id.size());

    const unsigned int n_bcs = node_id.size();

    // reserve an upper bound for the number of BCs
    xfer_node_bcs.reserve(2*n_bcs);

    // populate the xfer_node_bcs witl *local* nodes only
    for (unsigned int bc=0; bc<n_bcs; bc++)
      {
	const Node* node = mesh.node_ptr(node_id[bc]);
	
	assert(node != NULL);
	assert(node->id() == node_id[bc]);

	if (node->processor_id() == libMesh::processor_id())
	  {
	    xfer_node_bcs.push_back(node_id[bc]);
	    xfer_node_bcs.push_back(bc_id[bc]);
	  }
      }
  } // done with nodal boundary conditions


  // The xfer arrays now contain all the information for our
  // local bcs, and we are about to get all the information for
  // remote bcs.  Go ahead and clear the current boundary_info
  // information and rebuild it after we get the remote data.
  boundary_info.clear();
  
  // Get the boundary condition information from adacent processors
  Parallel::allgather (xfer_elem_bcs);
  Parallel::allgather (xfer_node_bcs);


  // Insert the elements
  {
    const unsigned int n_bcs = xfer_elem_bcs.size()/3;

    for (unsigned int bc=0, cnt=0; bc<n_bcs; bc++)
      {
	const Elem* elem = mesh.elem(xfer_elem_bcs[cnt++]);
	const unsigned short int side_id = xfer_elem_bcs[cnt++];
	const short int bc_id = xfer_elem_bcs[cnt++];

	boundary_info.add_side (elem, side_id, bc_id);	
      }

    // no need for this any more
    xfer_elem_bcs.resize(0);
  }

  
  // Insert the nodes
  {
    const unsigned int n_bcs = xfer_node_bcs.size()/2;

    for (unsigned int bc=0, cnt=0; bc<n_bcs; bc++)
      {
	const Node* node = mesh.node_ptr (xfer_node_bcs[cnt++]);
	const short int bc_id = xfer_node_bcs[cnt++];

	boundary_info.add_node (node, bc_id);
      }
  }

 
#ifndef NDEBUG
  
  // Make sure all processors agree on the number of boundary ids.
  const unsigned int n_bc_ids = boundary_info.n_boundary_ids();
  unsigned int global_n_bc_ids = n_bc_ids;
  
  Parallel::max (global_n_bc_ids);
  assert (n_bc_ids == global_n_bc_ids);
  
#endif
  

  STOP_LOG  ("allgather_bcs()","MeshCommunication");
}
#endif // HAVE_MPI



void MeshCommunication::delete_remote_elements(ParallelMesh& mesh) const
{
  START_LOG("delete_remote_elements()", "MeshCommunication");

  std::vector<bool> local_nodes(mesh.max_node_id(), false);
  std::vector<bool> semilocal_elems(mesh.max_elem_id(), false);

  // We don't want to delete any element that shares a node
  // with a local element.
  MeshBase::const_element_iterator l_elem_it = mesh.local_elements_begin(),
                                   l_end     = mesh.local_elements_end();
  for (; l_elem_it != l_end; ++l_elem_it)
    {
      const Elem *elem = *l_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        local_nodes[elem->node(n)] = true;
    }

  // We don't want to delete any element that shares a node
  // with an unpartitioned element either.
  MeshBase::const_element_iterator u_elem_it =
    mesh.pid_elements_begin(DofObject::invalid_processor_id),
                                   u_end     =
    mesh.pid_elements_end(DofObject::invalid_processor_id);
  for (; u_elem_it != u_end; ++u_elem_it)
    {
      const Elem *elem = *u_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        local_nodes[elem->node(n)] = true;
    }

  MeshBase::element_iterator nl_elem_it = mesh.not_local_elements_begin(),
                             nl_end     = mesh.not_local_elements_end();
  for (; nl_elem_it != nl_end; ++nl_elem_it)
    {
      Elem *elem = *nl_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        if (local_nodes[elem->node(n)])
          {
            semilocal_elems[elem->id()] = true;
            break;
          }
      if (!semilocal_elems[elem->id()])
        {
          // delete_elem doesn't currently invalidate element
          // iterators... that had better not change
          mesh.delete_elem(elem);
        }
    }

  STOP_LOG("delete_remote_elements()", "MeshCommunication");
}




// Pack all this information into one communication to avoid two latency hits
// For each element it is of the form
// [ level etype subdomain_id self_ID parent_ID node_0 node_1 ... node_n]
// We cannot use unsigned int because parent_ID can be negative
void MeshCommunication::pack_element (std::vector<int> &conn, const Elem* &elem) const
{
  assert (elem != NULL);
  
  conn.push_back (static_cast<int>(elem->level()));
  conn.push_back (static_cast<int>(elem->type()));
  conn.push_back (static_cast<int>(elem->subdomain_id()));
  conn.push_back (elem->id());
		
  // use parent_ID of -1 to indicate a level 0 element
  if (elem->level() == 0)
    conn.push_back(-1);
  else
    conn.push_back(elem->parent()->id());
  
  for (unsigned int n=0; n<elem->n_nodes(); n++)
    conn.push_back (elem->node(n));		
}
