// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "boundary_info.h"
#include "elem.h"
#include "libmesh_config.h"
#include "libmesh_common.h"
#include "libmesh_logging.h"
#include "location_maps.h"
#include "mesh_base.h"
#include "mesh_communication.h"
#include "mesh_tools.h"
#include "parallel.h"
#include "parallel_mesh.h"
#include "parallel_ghost_sync.h"
#include "utility.h"
#include "remote_elem.h"



//-----------------------------------------------
// anonymous namespace for implementation details
namespace {

  using libMesh::Elem;

  /**
   * Specific weak ordering for Elem*'s to be used in a set.
   * We use the id, but first sort by level.  This guarantees 
   * when traversing the set from beginning to end the lower 
   * level (parent) elements are encountered first.
   */
  struct CompareElemIdsByLevel
  {
    bool operator()(const Elem *a,
		    const Elem *b) const
    {
      libmesh_assert (a);
      libmesh_assert (b);
      const unsigned int
	al = a->level(), bl = b->level(),
	aid = a->id(),   bid = b->id();

      return (al == bl) ? aid < bid : al < bl;
    }
  };

  /**
   * Bitmask for which faces have boundary conditions.
   */
  const unsigned int 
    FaceNeighborMask[] = { 1,    // face 0 neighbor
			   2,    // face 1 neighbor
			   4,    // face 2 neighbor
			   8,    // face 3 neighbor
			   16,   // face 4 neighbor
			   32 }; // face 5 neighbor
}


namespace libMesh
{


// ------------------------------------------------------------
// MeshCommunication class members
void MeshCommunication::clear ()
{
  //  _neighboring_processors.clear();
}


// #ifdef LIBMESH_HAVE_MPI
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
//     libMesh::out << "Processor " << libMesh::processor_id()
// 	      << " intersects:" << std::endl;
//     for (unsigned int p=0; p<_neighboring_processors.size(); p++)
//       libMesh::out << " " << _neighboring_processors[p] << std::endl;
//   }
// }
// #else
// void MeshCommunication::find_neighboring_processors (const MeshBase&)
// {
// }
// #endif



#ifndef LIBMESH_HAVE_MPI // avoid spurious gcc warnings
// ------------------------------------------------------------
void MeshCommunication::redistribute (ParallelMesh &) const
{
  // no MPI == one processor, no redistribution
  return;
}
#else
// ------------------------------------------------------------
void MeshCommunication::redistribute (ParallelMesh &mesh) const
{
  // This method will be called after a new partitioning has been 
  // assigned to the elements.  This partitioning was defined in
  // terms of the active elements, and "trickled down" to the
  // parents and nodes as to be consistent.
  //
  // The point is that the entire concept of local elements is
  // kinda shaky in this method.  Elements which were previously
  // local may now be assigned to other processors, so we need to
  // send those off.  Similarly, we need to accept elements from
  // other processors.  
  // 
  // The approach is as follows:
  // (1) send all elements we have stored to their proper homes
  // (2) receive elements from all processors, watching for duplicates
  // (3) deleting all nonlocal elements elements
  // (4) obtaining required ghost elements from neighboring processors
  parallel_only();
  libmesh_assert (!mesh.is_serial());
  libmesh_assert (MeshTools::n_elem(mesh.unpartitioned_elements_begin(),
				    mesh.unpartitioned_elements_end()) == 0);
  
  START_LOG("redistribute()","MeshCommunication");

  // register a derived datatype to use in shipping nodes  
  Parallel::DataType packed_node_datatype = Node::PackedNode::create_mpi_datatype();
  packed_node_datatype.commit();

  // Get a few unique message tags to use in communications; we'll
  // default to some numbers around pi*1000
  Parallel::MessageTag
    nodestag   = Parallel::Communicator_World.get_unique_tag(3141),
    nodebcstag = Parallel::Communicator_World.get_unique_tag(3142),
    elemstag   = Parallel::Communicator_World.get_unique_tag(3143),
    elembcstag = Parallel::Communicator_World.get_unique_tag(3144);

  // Figure out how many nodes and elements we have which are assigned to each
  // processor.  send_n_nodes_and_elem_per_proc contains the number of nodes/elements
  // we will be sending to each processor, recv_n_nodes_and_elem_per_proc contains
  // the number of nodes/elements we will be receiving from each processor.
  // Format:
  //  send_n_nodes_and_elem_per_proc[5*pid+0] = number of nodes to send to pid
  //  send_n_nodes_and_elem_per_proc[5*pid+1] = number of elements to send to pid
  //  send_n_nodes_and_elem_per_proc[5*pid+2] = connectivity buffer size
  //  send_n_nodes_and_elem_per_proc[5*pid+3] = node bc buffer size
  //  send_n_nodes_and_elem_per_proc[5*pid+4] = element bc buffer size
  std::vector<unsigned int> send_n_nodes_and_elem_per_proc(5*libMesh::n_processors(), 0);
  
  std::vector<std::vector<Node::PackedNode> >
    nodes_sent(libMesh::n_processors());

  std::vector<std::vector<int> > 
    elements_sent(libMesh::n_processors()),
    node_bcs_sent(libMesh::n_processors()),
    element_bcs_sent(libMesh::n_processors());
    
  std::vector<Parallel::Request>
    node_send_requests, element_send_requests,
    node_bc_requests,   element_bc_requests;

  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
    if (pid != libMesh::processor_id()) // don't send to ourselves!!
      {
	// Build up a list of nodes and elements to send to processor pid.
	// We will certainly send all the active elements assigned to this processor,
	// but we will also ship off the other elements in the same family tree
	// as the active ones for data structure consistency.  We also
	// ship any nodes connected to these elements.  Note some of these nodes
	// and elements may be replicated from other processors, but that is OK.
	std::set<const Elem*, CompareElemIdsByLevel> elements_to_send;
	std::set<const Node*> connected_nodes;
	{
	  std::vector<const Elem*> family_tree;

	  MeshBase::const_element_iterator       elem_it  = mesh.active_pid_elements_begin(pid);
	  const MeshBase::const_element_iterator elem_end = mesh.active_pid_elements_end(pid);
	  
	  for (; elem_it!=elem_end; ++elem_it)
	    {
	      const Elem *top_parent = (*elem_it)->top_parent();

	      // avoid a lot of duplication -- if we already have top_parent
	      // in the set its entire family tree is already in the set.
	      if (!elements_to_send.count(top_parent))
		{
#ifdef LIBMESH_ENABLE_AMR
		  top_parent->family_tree(family_tree);
#else
		  family_tree.clear();
		  family_tree.push_back(top_parent);
#endif
		  for (unsigned int e=0; e<family_tree.size(); e++)
		    {
		      const Elem *elem = family_tree[e];
		      elements_to_send.insert (elem);
		      
		      for (unsigned int n=0; n<elem->n_nodes(); n++)
			connected_nodes.insert (elem->get_node(n));		  
		    }
		}

	      // then do the same for the face neighbors
	      for (unsigned int s=0; s<(*elem_it)->n_sides(); s++)
		if ((*elem_it)->neighbor(s) != NULL)
		  if (!(*elem_it)->neighbor(s)->is_remote())
		    {
		      top_parent = (*elem_it)->neighbor(s)->top_parent();
		      
		      if (!elements_to_send.count(top_parent))
			{
#ifdef LIBMESH_ENABLE_AMR
			  top_parent->family_tree(family_tree);
#else
			  family_tree.clear();
			  family_tree.push_back(top_parent);
#endif			  
			  for (unsigned int e=0; e<family_tree.size(); e++)
			    {
			      const Elem *elem = family_tree[e];
			      elements_to_send.insert (elem);
			      
			      for (unsigned int n=0; n<elem->n_nodes(); n++)
				connected_nodes.insert (elem->get_node(n));		  
			    }
			}
		    }
	    }
	}
	// The elements_to_send set now contains all the elements stored on the local
	// processor but owned by processor pid.  Additionally, the face neighbors
	// for these elements are also transferred.  Finally, the entire refinement
	// tree is also included.  This is a very simplistic way of ensuring data
	// structure consistency at the cost of larger communication buffers.  It is
	// worth profiling this on several parallel architectures to assess its impact.
	
	
	if (!connected_nodes.empty())
	  {
	    // the number of nodes we will ship to pid
	    send_n_nodes_and_elem_per_proc[5*pid+0] = connected_nodes.size();
	    
	    nodes_sent[pid].reserve(connected_nodes.size());
	    
	    for (std::set<const Node*>::const_iterator node_it = connected_nodes.begin();
		 node_it != connected_nodes.end(); ++node_it)
	      {
		nodes_sent[pid].push_back(Node::PackedNode(**node_it));

		// add the node if it has BCs
                std::vector<short int> bcs = mesh.boundary_info->boundary_ids(*node_it);
                
		if (!bcs.empty())
		  {
		    node_bcs_sent[pid].push_back((*node_it)->id());
		    node_bcs_sent[pid].push_back(bcs.size());
                    for(unsigned int bc_it=0; bc_it < bcs.size(); bc_it++)
                      node_bcs_sent[pid].push_back(bcs[bc_it]);
		  }
	      }

	    // the size of the node bc buffer we will ship to pid
	    send_n_nodes_and_elem_per_proc[5*pid+3] = node_bcs_sent[pid].size();	    
	    
	    // send the nodes off to the destination processor
	    node_send_requests.push_back(Parallel::request());
	  
	    Parallel::send (pid,
			    nodes_sent[pid],
			    packed_node_datatype,
			    node_send_requests.back(),
			    nodestag);

	    if (!node_bcs_sent[pid].empty())
	      {
		node_bc_requests.push_back(Parallel::request());

		Parallel::send (pid,
				node_bcs_sent[pid],
				node_bc_requests.back(),
				nodebcstag);
	      }
	  }
	
	if (!elements_to_send.empty())
	  {
	    // the number of elements we will send to this processor
	    send_n_nodes_and_elem_per_proc[5*pid+1] = elements_to_send.size();

	    for (std::set<const Elem*, CompareElemIdsByLevel>::const_iterator 
		   elem_it = elements_to_send.begin(); elem_it != elements_to_send.end(); ++elem_it)
	      {
		Elem::PackedElem::pack (elements_sent[pid], *elem_it);

		// if this is a level-0 element look for boundary conditions
	        if ((*elem_it)->level() == 0)
		  for (unsigned int s=0; s<(*elem_it)->n_sides(); s++)
		    if ((*elem_it)->neighbor(s) == NULL)
                      {
                        const std::vector<short int>& bc_ids = mesh.boundary_info->boundary_ids(*elem_it, s);
                        for (std::vector<short int>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
                          {
                            const short int bc_id = *id_it;
		            if (bc_id != mesh.boundary_info->invalid_id)
			      {
			        element_bcs_sent[pid].push_back ((*elem_it)->id());
			        element_bcs_sent[pid].push_back (s);
			        element_bcs_sent[pid].push_back (bc_id);
			      }
                          }
                      }
	      }
	    
	    // the packed connectivity size to send to this processor
	    send_n_nodes_and_elem_per_proc[5*pid+2] = elements_sent[pid].size();

	    // send the elements off to the destination processor
	    element_send_requests.push_back(Parallel::request());
	  
	    Parallel::send (pid,
			    elements_sent[pid],
			    element_send_requests.back(),
			    elemstag);

	    // the size of the element bc buffer we will ship to pid
	    send_n_nodes_and_elem_per_proc[5*pid+4] = element_bcs_sent[pid].size();

	    if (!element_bcs_sent[pid].empty())
	      {
		element_bc_requests.push_back(Parallel::request());

		Parallel::send (pid,
				element_bcs_sent[pid],
				element_bc_requests.back(),
				elembcstag);
	      }
	  }
      }
  
  std::vector<unsigned int> recv_n_nodes_and_elem_per_proc(send_n_nodes_and_elem_per_proc);

  Parallel::alltoall (recv_n_nodes_and_elem_per_proc);

  // In general we will only need to communicate with a subset of the other processors.
  // I can't immediately think of a case where we will send elements but not nodes, but
  // these are only bools and we have the information anyway...
  std::vector<bool>
    send_node_pair(libMesh::n_processors(),false), send_elem_pair(libMesh::n_processors(),false),
    recv_node_pair(libMesh::n_processors(),false), recv_elem_pair(libMesh::n_processors(),false);
  
  unsigned int
    n_send_node_pairs=0,     n_send_elem_pairs=0,
    n_send_node_bc_pairs=0,  n_send_elem_bc_pairs=0, 
    n_recv_node_pairs=0,     n_recv_elem_pairs=0,
    n_recv_node_bc_pairs=0,  n_recv_elem_bc_pairs=0, 
    max_n_nodes_received=0,  max_conn_size_received=0,
    max_node_bcs_received=0, max_elem_bcs_received=0;

  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
    {
      if (send_n_nodes_and_elem_per_proc[5*pid+0]) // we have nodes to send
	{
	  send_node_pair[pid] = true;
	  n_send_node_pairs++;

	  if (send_n_nodes_and_elem_per_proc[5*pid+3]) // and node bcs
	    n_send_node_bc_pairs++;
	}
      
      if (send_n_nodes_and_elem_per_proc[5*pid+1]) // we have elements to send
	{
	  send_elem_pair[pid] = true;
	  n_send_elem_pairs++;

	  if (send_n_nodes_and_elem_per_proc[5*pid+4]) // and element bcs
	    n_send_elem_bc_pairs++;
	}

      if (recv_n_nodes_and_elem_per_proc[5*pid+0]) // we have nodes to receive
	{
	  recv_node_pair[pid] = true;
	  n_recv_node_pairs++;
	  max_n_nodes_received = std::max(max_n_nodes_received,
					  recv_n_nodes_and_elem_per_proc[5*pid+0]);

	  if (recv_n_nodes_and_elem_per_proc[5*pid+3]) // and node bcs
	    {
	      n_recv_node_bc_pairs++;
	      max_node_bcs_received = std::max(max_node_bcs_received,
					       recv_n_nodes_and_elem_per_proc[5*pid+3]);
	    }		
	}
      
      if (recv_n_nodes_and_elem_per_proc[5*pid+1]) // we have elements to receive
	{
	  recv_elem_pair[pid] = true;
	  n_recv_elem_pairs++;
	  max_conn_size_received = std::max(max_conn_size_received,
					    recv_n_nodes_and_elem_per_proc[5*pid+2]);

	  if (recv_n_nodes_and_elem_per_proc[5*pid+4]) // and element bcs
	    {
	      n_recv_elem_bc_pairs++;
	      max_elem_bcs_received = std::max(max_elem_bcs_received,
					       recv_n_nodes_and_elem_per_proc[5*pid+4]);
	    }
	}
    }
  libmesh_assert (n_send_node_pairs    == node_send_requests.size());
  libmesh_assert (n_send_elem_pairs    == element_send_requests.size());
  libmesh_assert (n_send_node_bc_pairs == node_bc_requests.size());
  libmesh_assert (n_send_elem_bc_pairs == element_bc_requests.size());

  // Receive nodes.  Size this array for the largest message.
  {
    std::vector<Node::PackedNode> received_nodes; /**/
    received_nodes.reserve(max_n_nodes_received);

    // We now know how many processors will be sending us nodes.
    for (unsigned int node_comm_step=0; node_comm_step<n_recv_node_pairs; node_comm_step++)
      {
	// but we don't necessarily want to impose an ordering, so
	// just grab whatever message is available next.  note that
	// the blocking receive should resize the vector to fit the
	// incoming message
	Parallel::Status status =
	  Parallel::receive (Parallel::any_source,
			     received_nodes,
			     packed_node_datatype,
			     nodestag);
	const unsigned int source_pid = status.source();
	const unsigned int n_nodes_received =
	  recv_n_nodes_and_elem_per_proc[5*source_pid+0];
	libmesh_assert (n_nodes_received == received_nodes.size());
	libmesh_assert (status.size() == n_nodes_received);
	libmesh_assert (recv_node_pair[source_pid]);
	
	for (unsigned int n=0; n<n_nodes_received; n++)
	  {
	    Node *node = received_nodes[n].build_node().release();
	    mesh.insert_node(node); // insert_node works even if the
	  }                         // node already exists in the mesh,
      }                             // in which case it overwrites (x,y,z)
  }

  // Receive BCs for those nodes.
  {
    std::vector<int> received_node_bcs(max_node_bcs_received);

    // again, receive all the messages in whatever order they come.
    for (unsigned int node_comm_step=0; node_comm_step<n_recv_node_bc_pairs; node_comm_step++)
      {
	Parallel::Status status =
	  Parallel::receive (Parallel::any_source,
			     received_node_bcs,
			     nodebcstag);

	libmesh_assert (status.size() == recv_n_nodes_and_elem_per_proc[5*status.source()+3]);

	const unsigned int buffer_size = status.size();
	unsigned int cnt=0;

	while (cnt < buffer_size)
	  {
	    const unsigned int node_id = received_node_bcs[cnt++];
	    const unsigned int num_bcs = received_node_bcs[cnt++];

            libmesh_assert (mesh.node_ptr(node_id));
            
            for(unsigned int bc_it=0; bc_it < num_bcs; bc_it++)
              mesh.boundary_info->add_node (mesh.node_ptr(node_id), received_node_bcs[cnt++]);
	  }	
      }
  }


  
  // Receive elements.  Size this array for the largest message.
  {
    std::vector<int> received_elements(max_conn_size_received);

    // Similarly we know how many processors are sending us elements, 
    // but we don't really care in what order we receive them.
    for (unsigned int elem_comm_step=0; elem_comm_step<n_recv_elem_pairs; elem_comm_step++)
      {
	Parallel::Status status =
	  Parallel::receive (Parallel::any_source,
			     received_elements,
			     elemstag);

#ifndef NDEBUG
	// Avoid declaring these variables unless asserts are enabled.
	const unsigned int source_pid = status.source();
	const unsigned int n_elem_received = 
	  recv_n_nodes_and_elem_per_proc[5*source_pid+1];
#endif
	libmesh_assert (recv_elem_pair[source_pid]);
	libmesh_assert (recv_n_nodes_and_elem_per_proc[5*source_pid+2] 
		<= received_elements.size());
	libmesh_assert (recv_n_nodes_and_elem_per_proc[5*source_pid+2] 
		== status.size());

	// iterate through the input buffer and add the elements
	const unsigned int xfer_buffer_size = status.size();
	unsigned int cnt=0;
	unsigned int current_elem=0;
	while (cnt < xfer_buffer_size)
	  {
	    // Unpack the element 
	    Elem::PackedElem packed_elem (received_elements.begin()+cnt);
	    
	    // The ParallelMesh::elem(i) member will return NULL if the element
	    // is not in the mesh.  We rely on that here, so it better not change!
	    Elem *elem = mesh.elem(packed_elem.id());
	    
	    // if we already have this element, make sure its properties match
	    // but then go on
	    if (elem)
	      {
		libmesh_assert (elem->level()             == packed_elem.level());
		libmesh_assert (elem->id()                == packed_elem.id());
		libmesh_assert (elem->processor_id()      == packed_elem.processor_id());
		libmesh_assert (elem->subdomain_id()      == packed_elem.subdomain_id());
		libmesh_assert (elem->type()              == packed_elem.type());
#ifdef LIBMESH_ENABLE_AMR
		libmesh_assert (elem->p_level()           == packed_elem.p_level());
		libmesh_assert (elem->refinement_flag()   == packed_elem.refinement_flag());
		libmesh_assert (elem->p_refinement_flag() == packed_elem.p_refinement_flag());
		
		if (elem->level() > 0)
		  {
		    libmesh_assert (elem->parent()->id() == static_cast<unsigned int>(packed_elem.parent_id()));
		    libmesh_assert (elem->parent()->child(packed_elem.which_child_am_i()) == elem);
		  }	      
#endif
		libmesh_assert (elem->n_nodes() == packed_elem.n_nodes());
	      }
	    else
	      {
		// We need to add the element.
#ifdef LIBMESH_ENABLE_AMR
		// maybe find the parent
		if (packed_elem.level() > 0)
		  {
		    Elem *parent = mesh.elem(packed_elem.parent_id());
		    
		    // Note that we were very careful to construct the send connectivity
		    // so that parents are encountered before children.  So if we get here
		    // and can't find the parent that is a fatal error.
		    if (parent == NULL)
		      {
			libMesh::err << "Parent element with ID " << packed_elem.parent_id()
				      << " not found." << std::endl; 
			libmesh_error();
		      }
		    
		    elem = packed_elem.unpack (mesh, parent);
		  }
		else
		  {
		    libmesh_assert (packed_elem.parent_id() == -1);
#endif // LIBMESH_ENABLE_AMR
		    elem = packed_elem.unpack (mesh);
#ifdef LIBMESH_ENABLE_AMR
		  }
#endif		
		// Good to go.  Add to the mesh.
		libmesh_assert (elem);
		libmesh_assert (elem->n_nodes() == packed_elem.n_nodes());
		mesh.insert_elem(elem);
	      }
	    
	    // properly position cnt for the next element 
	    cnt += Elem::PackedElem::header_size + packed_elem.n_nodes();
	    current_elem++;
	  }
	libmesh_assert (current_elem == n_elem_received);
      }
  }

  // Receive boundary conditions for those elements
  {
    std::vector<int> received_element_bcs(max_elem_bcs_received);
    
    // again, receive all the messages in whatever order they come.
    for (unsigned int elem_comm_step=0; elem_comm_step<n_recv_elem_bc_pairs; elem_comm_step++)
      {
	Parallel::Status status =
	  Parallel::receive (Parallel::any_source,
			     received_element_bcs,
			     elembcstag);

	libmesh_assert (status.size() == recv_n_nodes_and_elem_per_proc[5*status.source()+4]);

	const unsigned int buffer_size = status.size();
	unsigned int cnt=0;

	while (cnt < buffer_size)
	  {
	    const unsigned int elem_id = received_element_bcs[cnt++];
	    const unsigned int side    = received_element_bcs[cnt++];
	    const int bc_id            = received_element_bcs[cnt++];

	    libmesh_assert (mesh.elem(elem_id));
	    mesh.boundary_info->add_side (mesh.elem(elem_id), side, bc_id);
	  }
      }    
  }

  // Wait for all sends to complete
  Parallel::wait (node_send_requests);
  Parallel::wait (node_bc_requests);
  Parallel::wait (element_send_requests);
  Parallel::wait (element_bc_requests);
  
  // unregister MPI datatypes
  packed_node_datatype.free();

  // Check on the redistribution consistency
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_refinement_tree(mesh);
#endif
  
  STOP_LOG("redistribute()","MeshCommunication");  
}
#endif // LIBMESH_HAVE_MPI



#ifndef LIBMESH_HAVE_MPI // avoid spurious gcc warnings
// ------------------------------------------------------------
void MeshCommunication::gather_neighboring_elements (ParallelMesh &) const
{
  // no MPI == one processor, no need for this method...
  return;
}
#else
// ------------------------------------------------------------
void MeshCommunication::gather_neighboring_elements (ParallelMesh &mesh) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;  

  // This function must be run on all processors at once
  parallel_only();
  
  START_LOG("gather_neighboring_elements()","MeshCommunication");

  //------------------------------------------------------------------
  // The purpose of this function is to provide neighbor data structure
  // consistency for a parallel, distributed mesh.  In libMesh we require
  // that each local element have access to a full set of valid face
  // neighbors.  In some cases this requires us to store "ghost elements" -
  // elements that belong to other processors but we store to provide
  // data structure consistency.  Also, it is assumed that any element
  // with a NULL neighbor resides on a physical domain boundary.  So,
  // even our "ghost elements" must have non-NULL neighbors.  To handle
  // this the concept of "RemoteElem" is used - a special construct which
  // is used to denote that an element has a face neighbor, but we do
  // not actually store detailed information about that neighbor.  This
  // is required to prevent data structure explosion.
  //
  // So when this method is called we should have only local elements.
  // These local elements will then find neighbors among the local
  // element set.  After this is completed, any element with a NULL
  // neighbor has either (i) a face on the physical boundary of the mesh,
  // or (ii) a neighboring element which lives on a remote processor.
  // To handle case (ii), we communicate the global node indices connected
  // to all such faces to our neighboring processors.  They then send us
  // all their elements with a NULL neighbor that are connected to any
  // of the nodes in our list.
  //------------------------------------------------------------------

  // Let's begin with finding consistent neighbor data information
  // for all the elements we currently have.  We'll use a clean
  // slate here - clear any existing information, including RemoteElem's.
  mesh.find_neighbors (/* reset_remote_elements = */ true,
		       /* reset_current_list    = */ true);    
  
  // register a derived datatype to use in shipping nodes  
  Parallel::DataType packed_node_datatype = Node::PackedNode::create_mpi_datatype();
  packed_node_datatype.commit();

  // Get a few unique message tags to use in communications; we'll
  // default to some numbers around pi*10000
  Parallel::MessageTag
    node_list_tag         = Parallel::Communicator_World.get_unique_tag(31415),
    element_neighbors_tag = Parallel::Communicator_World.get_unique_tag(31416);

  // Now any element with a NULL neighbor either
  // (i) lives on the physical domain boundary, or
  // (ii) lives on an inter-processor boundary.
  // We will now gather all the elements from adjacent processors
  // which are of the same state, which should address all the type (ii)
  // elements.
  
  // A list of all the processors which *may* contain neighboring elements.
  // (for development simplicity, just make this the identity map)
  std::vector<unsigned int> adjacent_processors;
  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
    if (pid != libMesh::processor_id())
      adjacent_processors.push_back (pid);


  const unsigned int n_adjacent_processors = adjacent_processors.size();

  //-------------------------------------------------------------------------
  // Let's build a list of all nodes which live on NULL-neighbor sides.
  // For simplicity, we will use a set to build the list, then transfer
  // it to a vector for communication.
  std::vector<unsigned int> my_interface_node_list;
  std::vector<const Elem*>  my_interface_elements;
  {
    std::set<unsigned int> my_interface_node_set;

    // since parent nodes are a subset of children nodes, this should be sufficient
    MeshBase::const_element_iterator       it     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator it_end = mesh.active_local_elements_end();

    for (; it != it_end; ++it)
      {
	const Elem * const elem = *it;
	libmesh_assert (elem != NULL);

	if (elem->on_boundary()) // denotes *any* side has a NULL neighbor
	  {
	    my_interface_elements.push_back(elem); // add the element, but only once, even
	                                           // if there are multiple NULL neighbors
	    for (unsigned int s=0; s<elem->n_sides(); s++)
	      if (elem->neighbor(s) == NULL)
		{
		  AutoPtr<Elem> side(elem->build_side(s));
		  
		  for (unsigned int n=0; n<side->n_vertices(); n++)
		    my_interface_node_set.insert (side->node(n));   
		}
	  }
      }

    my_interface_node_list.reserve (my_interface_node_set.size());
    my_interface_node_list.insert  (my_interface_node_list.end(),
				    my_interface_node_set.begin(),
				    my_interface_node_set.end());
  }
  
  if (false)
    libMesh::out << "[" << libMesh::processor_id() << "] "
	          << "mesh.n_nodes()=" << mesh.n_nodes() << ", "
	          << "my_interface_node_list.size()=" << my_interface_node_list.size()
	          << std::endl;
  
  // we will now send my_interface_node_list to all of the adjacent processors.
  // note that for the time being we will copy the list to a unique buffer for 
  // each processor so that we can use a nonblocking send and not access the
  // buffer again until the send completes.  it is my understanding that the
  // MPI 2.1 standard seeks to remove this restriction as unnecessary, so in
  // the future we should change this to send the same buffer to each of the
  // adjacent processors. - BSK 11/17/2008
  std::vector<std::vector<unsigned int> > 
    my_interface_node_xfer_buffers (n_adjacent_processors, my_interface_node_list);
  std::vector<Parallel::Request> my_interface_node_list_requests (n_adjacent_processors);
  std::map<unsigned int, unsigned int> n_comm_steps;

  for (unsigned int comm_step=0; comm_step<n_adjacent_processors; comm_step++)
    {
      n_comm_steps[adjacent_processors[comm_step]]=1;
      Parallel::send (adjacent_processors[comm_step],
		      my_interface_node_xfer_buffers[comm_step],
		      my_interface_node_list_requests[comm_step],
		      node_list_tag);
    }

  //-------------------------------------------------------------------------
  // processor pairings are symmetric - I expect to receive an interface node
  // list from each processor in adjacent_processors as well!
  // now we will catch an incoming node list for each of our adjacent processors.
  //
  // we are done with the adjacent_processors list - note that it is in general
  // a superset of the processors we truly share elements with.  so let's
  // clear the superset list, and we will fill it with the true list.
  adjacent_processors.clear();
  
  std::vector<unsigned int> common_interface_node_list;

  // send buffers.  we will fill these with data from the elements we own 
  // which share nodes with an adjacent processor.  we will slightly abuse 
  // the node_bcs_sent and element_bcs_sent buffers - we'll fill them and then tack
  // them on to the end of the elements_sent buffer to reduce the message count.
  std::vector<std::vector<Node::PackedNode> > nodes_sent(n_adjacent_processors);
  std::vector<std::vector<int> > elements_sent(n_adjacent_processors);
  std::vector<int> node_bcs_sent, &element_bcs_sent(node_bcs_sent);
  std::vector<std::vector<unsigned int> > element_neighbors_sent(n_adjacent_processors);
    
  std::vector<Parallel::Request> 
    node_send_requests(n_adjacent_processors), 
    element_send_requests(n_adjacent_processors),
    element_neighbor_send_requests(n_adjacent_processors);

  // receive buffers
  std::vector<Node::PackedNode> nodes_received;
  std::vector<int> elements_received;
  std::vector<unsigned int> element_neighbors_received;

  // we expect two classess of messages - 
  // (1) incoming interface node lists, to which we will reply with our elements 
  //     touching nodes in the list, and
  // (2) replies from the requests we sent off previously.  
  //  (2.a) - nodes
  //  (2.b) - element connectivity & bc info
  // so we expect 3 communications from each adjacent processor.
  // by structuring the communication in this way we hopefully impose no 
  // order on the handling of the arriving messages.  in particular, we
  // should be able to handle the case where we receive a request and
  // all replies from processor A before even receiving a request from
  // processor B.
  unsigned int 
    n_node_replies_sent=0, n_node_replies_received=0, 
    n_elem_replies_sent=0, n_elem_replies_received=0;

  for (unsigned int comm_step=0; comm_step<3*n_adjacent_processors; comm_step++)
    {
      //------------------------------------------------------------------
      // catch incoming node list
      Parallel::Status
	status(Parallel::probe (Parallel::any_source,
				node_list_tag));      
      const unsigned int
	source_pid_idx = status.source(),
	dest_pid_idx   = source_pid_idx;

      //------------------------------------------------------------------
      // first time - incoming request
      if (n_comm_steps[source_pid_idx] == 1)
	{
	  n_comm_steps[source_pid_idx]++;
	  
	  Parallel::receive (source_pid_idx,
			     common_interface_node_list,
			     node_list_tag);
	  const unsigned int	
	    their_interface_node_list_size = common_interface_node_list.size();
	  
	  // we now have the interface node list from processor source_pid_idx.
	  // now we can find all of our elements which touch any of these nodes
	  // and send copies back to this processor.  however, we can make our
	  // search more efficient by first excluding all the nodes in
	  // their list which are not also contained in
	  // my_interface_node_list.  we can do this in place as a set
	  // intersection.
	  common_interface_node_list.erase
	    (std::set_intersection (my_interface_node_list.begin(),
				    my_interface_node_list.end(),
				    common_interface_node_list.begin(),
				    common_interface_node_list.end(),
				    common_interface_node_list.begin()),
	     common_interface_node_list.end());
	  
	  if (false)
	    libMesh::out << "[" << libMesh::processor_id() << "] "
		          << "my_interface_node_list.size()="       << my_interface_node_list.size()
		          << ", [" << source_pid_idx << "] "
		          << "their_interface_node_list.size()="    << their_interface_node_list_size
		          << ", common_interface_node_list.size()=" << common_interface_node_list.size()
		          << std::endl;
	  
	  // Check for quick return?
	  if (common_interface_node_list.empty())
	    {
	      // let's try to be smart here - if we have no nodes in common,
	      // we cannot share elements.  so post the messages expected
	      // from us here and go on about our business.  
	      // note that even though these are nonblocking sends
	      // they should complete essentially instantly, because
	      // in all cases the send buffers are empty
	      Parallel::send (dest_pid_idx,
			      nodes_sent[n_node_replies_sent],
			      packed_node_datatype,
			      node_send_requests[n_node_replies_sent],
			      node_list_tag);
	      n_node_replies_sent++;

	      Parallel::send (dest_pid_idx,
			      elements_sent[n_elem_replies_sent],
			      element_send_requests[n_elem_replies_sent],
			      node_list_tag);

	      Parallel::send (dest_pid_idx,
			      element_neighbors_sent[n_elem_replies_sent],
			      element_neighbor_send_requests[n_elem_replies_sent],
			      element_neighbors_tag);
	      n_elem_replies_sent++;

	      continue;
	    }
	  // otherwise, this really *is* an adjacent processor.
	  adjacent_processors.push_back(source_pid_idx);

	  // Now we need to see which of our elements touch the nodes in the list.
	  // We built a reduced element list above, and we know the pointers are
	  // not NULL.
	  // We will certainly send all the active elements which intersect source_pid_idx,
	  // but we will also ship off the other elements in the same family tree
	  // as the active ones for data structure consistency.  We also
	  // ship any nodes connected to these elements.  Note some of these nodes
	  // and elements may be replicated from other processors, but that is OK.
	  std::set<const Elem*, CompareElemIdsByLevel> elements_to_send;
	  std::set<const Node*> connected_nodes;
	  std::vector<const Elem*> family_tree;
	  
	  for (unsigned int e=0, n_shared_nodes=0; e<my_interface_elements.size(); e++, n_shared_nodes=0)
	    {
	      const Elem * elem = my_interface_elements[e];
	      
	      for (unsigned int n=0; n<elem->n_vertices(); n++)
		if (std::binary_search (common_interface_node_list.begin(),
					common_interface_node_list.end(),
					elem->node(n)))
		  { 
		    n_shared_nodes++;
		    
		    // TBD - how many nodes do we need to share 
		    // before we care?  certainly 2, but 1?  not 
		    // sure, so let's play it safe...
		    if (n_shared_nodes > 0) break;
		  }
	      
	      if (n_shared_nodes) // share at least one node?
		{
		  elem = elem->top_parent();
		  
		  // avoid a lot of duplicated effort -- if we already have elem
		  // in the set its entire family tree is already in the set.
		  if (!elements_to_send.count(elem))
		    {
#ifdef LIBMESH_ENABLE_AMR
		      elem->family_tree(family_tree);
#else
		      family_tree.clear();
		      family_tree.push_back(elem);
#endif
		      for (unsigned int leaf=0; leaf<family_tree.size(); leaf++)
			{
			  elem = family_tree[leaf];
			  elements_to_send.insert (elem);
			  
			  for (unsigned int n=0; n<elem->n_nodes(); n++)
			    connected_nodes.insert (elem->get_node(n));		  
			}		    
		    }
		}
	    }
	  
	  // The elements_to_send and connected_nodes sets now contain all
	  // the elements and nodes we need to send to this processor.
	  // All that remains is to pack up the objects (along with
	  // any boundary conditions) and send the messages off.
	  {
	    if (elements_to_send.empty()) libmesh_assert (connected_nodes.empty());
	    if (connected_nodes.empty())  libmesh_assert (elements_to_send.empty());
	    
	    const unsigned int n_nodes_sent = connected_nodes.size();
	    nodes_sent[n_node_replies_sent].reserve (n_nodes_sent);
	    node_bcs_sent.clear();
	
	    for (std::set<const Node*>::const_iterator node_it = connected_nodes.begin();
		 node_it != connected_nodes.end(); ++node_it)
	      {
		nodes_sent[n_node_replies_sent].push_back (Node::PackedNode(**node_it));
		
		// add the node if it has BCs
                std::vector<short int> bcs = mesh.boundary_info->boundary_ids(*node_it);
                
		if (!bcs.empty())
		  {
		    node_bcs_sent.push_back((*node_it)->id());
		    node_bcs_sent.push_back(bcs.size());
                    for(unsigned int bc_it=0; bc_it < bcs.size(); bc_it++)
                      node_bcs_sent.push_back(bcs[bc_it]);
		  }
	      }
	    connected_nodes.clear();
	
	    // send the nodes off to the destination processor
	    Parallel::send (dest_pid_idx,
			    nodes_sent[n_node_replies_sent],
			    packed_node_datatype,
			    node_send_requests[n_node_replies_sent],
			    node_list_tag);
	    n_node_replies_sent++;
	
	    // let's pack the node bc data in the front of the element
	    // connectivity buffer
	    const unsigned int
	      n_node_bcs_sent = node_bcs_sent.size();
	    
	    elements_sent[n_elem_replies_sent].clear();
	    elements_sent[n_elem_replies_sent].insert (elements_sent[n_elem_replies_sent].end(),
						       node_bcs_sent.begin(),
						       node_bcs_sent.end());
	
	    // this is really just a reference to node_bcs_sent,
	    // so be careful not to clear it before packing the
	    // node_bcs_sent into the elements_sent buffer!!
	    const unsigned int n_elements_sent = elements_to_send.size();
	    element_neighbors_sent[n_elem_replies_sent].reserve(2*n_elements_sent+1);
	    element_bcs_sent.clear(); 
	    
	    for (std::set<const Elem*, CompareElemIdsByLevel>::const_iterator 
		   elem_it = elements_to_send.begin(); elem_it != elements_to_send.end(); ++elem_it)
	      {
		Elem::PackedElem::pack (elements_sent[n_elem_replies_sent], *elem_it);
		
		// let's put the element index in the element_neighbors_sent
		// buffer - we will later overwrite it with a bitmask of the
		// element neighbor information to send out.  this will be used 
		// to create RemoteElems where appropriate.
		element_neighbors_sent[n_elem_replies_sent].push_back((*elem_it)->id());
		element_neighbors_sent[n_elem_replies_sent].push_back(0);
		
		// if this is a level-0 element look for boundary conditions
		if ((*elem_it)->level() == 0)
		  for (unsigned int s=0; s<(*elem_it)->n_sides(); s++)
		    if ((*elem_it)->neighbor(s) == NULL)
                      {
                        const std::vector<short int>& bc_ids = mesh.boundary_info->boundary_ids(*elem_it, s);
                        for (std::vector<short int>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
                          {
                            const short int bc_id = *id_it;
		            if (bc_id != mesh.boundary_info->invalid_id)
			      {
			        element_bcs_sent.push_back ((*elem_it)->id());
			        element_bcs_sent.push_back (s);
			        element_bcs_sent.push_back (bc_id);
			      }
                          }
                      }
	      }
	    element_neighbors_sent[n_elem_replies_sent].push_back (dest_pid_idx);
	    elements_to_send.clear();
	    
	    // let's pack the node bc data in the front of the element
	    // connectivity buffer
	    const unsigned int
	      n_elem_bcs_sent = element_bcs_sent.size() / 3;
	    
	    elements_sent[n_elem_replies_sent].insert (elements_sent[n_elem_replies_sent].end(),
						       element_bcs_sent.begin(),
						       element_bcs_sent.end());
	    
	    // only append to the message if it is not empty
	    if (!elements_sent[n_elem_replies_sent].empty())
	      {
		// let's tack three ints on to
		// the end of the elements_sent
		// buffer for use on the receiving end
		elements_sent[n_elem_replies_sent].push_back (n_elements_sent);
		elements_sent[n_elem_replies_sent].push_back (n_elem_bcs_sent);
		elements_sent[n_elem_replies_sent].push_back (n_node_bcs_sent);
	      }
	    
	    // send the elements off to the destination processor
	    Parallel::send (dest_pid_idx,
			    elements_sent[n_elem_replies_sent],
			    element_send_requests[n_elem_replies_sent],
			    node_list_tag);
	    n_elem_replies_sent++;
	  }
	}



      //------------------------------------------------------------------
      // second time - reply of nodes
      else if (n_comm_steps[source_pid_idx] == 2)
	{
	  n_comm_steps[source_pid_idx]++;
	
	  Parallel::receive (source_pid_idx,
			     nodes_received,
			     packed_node_datatype,
			     node_list_tag);
	  n_node_replies_received++;
	  
	  // add the nodes we just received
	  for (unsigned int n=0; n<nodes_received.size(); n++)
	    mesh.insert_node (nodes_received[n].build_node().release());
	}



      //------------------------------------------------------------------
      // third time - elements & bcs
      else if (n_comm_steps[source_pid_idx] == 3)
	{ 
	  n_comm_steps[source_pid_idx]++;

	  Parallel::receive (source_pid_idx,
			     elements_received,
			     node_list_tag);	  
	  n_elem_replies_received++;
	  
	  if (elements_received.empty())
	    continue;
	  
	  libmesh_assert (elements_received.size() > 3);
	  
	  ////////////////////////////////////////////////////////////////////////////////
	  // remember: elements_sent = 
	  //          { [node bd data] [element conntectivity] [element bc data] # # # }
	  ////////////////////////////////////////////////////////////////////////////////
	  const unsigned int n_node_bcs_received = elements_received.back(); elements_received.pop_back();
	  const unsigned int n_elem_bcs_received = elements_received.back(); elements_received.pop_back();
	  const unsigned int n_elements_received = elements_received.back(); elements_received.pop_back();
	  
	  // counter into the bc/element connectivty buffer
	  unsigned int cnt=0;
	  
	  // add any node bcs
	  while (cnt < n_node_bcs_received)
	    {
              const unsigned int node_id = elements_received[cnt++];
              const unsigned int num_bcs = elements_received[cnt++];

              libmesh_assert (mesh.node_ptr(node_id));
            
              for(unsigned int bc_it=0; bc_it < num_bcs; bc_it++)
                mesh.boundary_info->add_node (mesh.node_ptr(node_id), elements_received[cnt++]);
	    }
      
	  // add the elements we just received
	  for (unsigned int e=0; e<n_elements_received; e++)
	    {
	      // Unpack the element 
	      Elem::PackedElem packed_elem (elements_received.begin()+cnt);
	  
	      // The ParallelMesh::elem(i) member will return NULL if the element
	      // is not in the mesh.  We rely on that here, so it better not change!
	      Elem *elem = mesh.elem(packed_elem.id());
	  
	      // if we already have this element, make sure its properties match
	      // but then go on
	      if (elem)
		{
		  libmesh_assert (elem->level()             == packed_elem.level());
		  libmesh_assert (elem->id()                == packed_elem.id());
		  libmesh_assert (elem->processor_id()      == packed_elem.processor_id());
		  libmesh_assert (elem->subdomain_id()      == packed_elem.subdomain_id());
		  libmesh_assert (elem->type()              == packed_elem.type());
#ifdef LIBMESH_ENABLE_AMR
		  libmesh_assert (elem->p_level()           == packed_elem.p_level());
		  libmesh_assert (elem->refinement_flag()   == packed_elem.refinement_flag());
		  libmesh_assert (elem->p_refinement_flag() == packed_elem.p_refinement_flag());
		  
		  if (elem->level() > 0)
		    {
		      libmesh_assert (elem->parent()->id() == static_cast<unsigned int>(packed_elem.parent_id()));
		      libmesh_assert (elem->parent()->child(packed_elem.which_child_am_i()) == elem);
		    }	      
#endif
		  libmesh_assert (elem->n_nodes() == packed_elem.n_nodes());
		}
	      else
		{
		  // We need to add the element.
#ifdef LIBMESH_ENABLE_AMR
		  // maybe find the parent
		  if (packed_elem.level() > 0)
		    {
		      Elem *parent = mesh.elem(packed_elem.parent_id());
		      
		      // Note that we were very careful to construct the send connectivity
		      // so that parents are encountered before children.  So if we get here
		      // and can't find the parent that is a fatal error.
		      if (parent == NULL)
			{
			  libMesh::err << "Parent element with ID " << packed_elem.parent_id()
				        << " not found." << std::endl; 
			  libmesh_error();
			}
		      
		      elem = packed_elem.unpack (mesh, parent);
		    }
		  else
		    {
		      libmesh_assert (packed_elem.parent_id() == -1);
		      elem = packed_elem.unpack (mesh);
		    }
#else // !defined(LIBMESH_ENABLE_AMR)
		  elem = packed_elem.unpack (mesh);
#endif		
		  // Good to go.  Add to the mesh.
		  libmesh_assert (elem);
		  libmesh_assert (elem->n_nodes() == packed_elem.n_nodes());
		  mesh.insert_elem(elem);
		}
	  
	      // properly position cnt for the next element 
	      cnt += Elem::PackedElem::header_size + packed_elem.n_nodes();	      
	    } // done adding elements
      
	  // add any element bcs
	  unsigned int n_elem_bcs=0;
	  while (cnt < elements_received.size())
	    {
	      n_elem_bcs++;
	      const unsigned int elem_id = elements_received[cnt++];
	      const unsigned int side    = elements_received[cnt++];
	      const int bc_id            = elements_received[cnt++];
	      
	      libmesh_assert (mesh.elem(elem_id));
	      mesh.boundary_info->add_side (mesh.elem(elem_id), side, bc_id);
	    }
	  libmesh_assert (n_elem_bcs_received == n_elem_bcs);
	}

      // Huh?
      else
	{
	  libMesh::err << "ERROR:  unexpected number of replies: "
		        << n_comm_steps[source_pid_idx]
		    << std::endl;
	}
    } // done catching & processing replies associated with tag ~ 100,000pi
  
  // Update neighbor information with the new elements,
  // but don't throw away old information.
  mesh.find_neighbors (/* reset_remote_elements = */ true,
		       /* reset_current_list    = */ false);    
  
  // OK, now at least all local elements should have a full set
  // of neighbor information.  we still need to resolve remote
  // element neighbor information though. to do this, we need
  // to get the neighbor information from the processors that
  // sent elements to us.
  // 
  // first, send element face information to requisite processors
  for (unsigned int comm_step=0; comm_step<n_adjacent_processors; comm_step++)
    if (!element_neighbors_sent[comm_step].empty()) // if it is empty we already sent it!
      {
	const unsigned int dest_pid_idx = 
	  element_neighbors_sent[comm_step].back();
	element_neighbors_sent[comm_step].pop_back();
	
	for (unsigned int pos=0; pos<element_neighbors_sent[comm_step].size(); pos+=2)
	  {
	    const unsigned int elem_global_idx = element_neighbors_sent[comm_step][pos+0];
	    const Elem *elem = mesh.elem(elem_global_idx);
	    unsigned int face_neighbor_bitmask = 0;
	    
	    for (unsigned int s=0; s<elem->n_sides(); s++)
	      if (elem->neighbor(s) != NULL)
		face_neighbor_bitmask |= FaceNeighborMask[s];
	    
	    element_neighbors_sent[comm_step][pos+1] = face_neighbor_bitmask;
	  }
      
	Parallel::send (dest_pid_idx,
			element_neighbors_sent[comm_step],
			element_neighbor_send_requests[comm_step],
			element_neighbors_tag);
      } // done sending neighbor information

  // second, receive and process all face neighbor information
  // from the processors which sent us elements
  for (unsigned int comm_step=0; comm_step<n_elem_replies_received; comm_step++)
    {
      Parallel::receive (Parallel::any_source,
			 element_neighbors_received,
			 element_neighbors_tag);

      for (unsigned int pos=0; pos<element_neighbors_received.size(); pos+=2)
	{
	  const unsigned int 
	    elem_global_idx       = element_neighbors_received[pos+0],
	    face_neighbor_bitmask = element_neighbors_received[pos+1];

	  Elem *elem = mesh.elem(elem_global_idx);

	  for (unsigned int s=0; s<elem->n_sides(); s++)
	    if ((face_neighbor_bitmask & FaceNeighborMask[s]) // if there is a face neighbor...
		&& (elem->neighbor(s) == NULL))               // ...and we do not have it
	      elem->set_neighbor(s,const_cast<RemoteElem*>(remote_elem)); // then it is a remote element
	}      
    } // done catching & processing neighbor information

  // allow any pending requests to complete
  Parallel::wait (my_interface_node_list_requests);
  Parallel::wait (node_send_requests);
  Parallel::wait (element_send_requests);
  Parallel::wait (element_neighbor_send_requests);  

  // unregister MPI datatypes
  packed_node_datatype.free();
  
  STOP_LOG("gather_neighboring_elements()","MeshCommunication");  
}
#endif // LIBMESH_HAVE_MPI



// ------------------------------------------------------------
void MeshCommunication::broadcast (MeshBase& mesh) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;
  
  this->broadcast_mesh (mesh);
  this->broadcast_bcs  (mesh, *(mesh.boundary_info));
}



#ifndef LIBMESH_HAVE_MPI
// ------------------------------------------------------------
void MeshCommunication::broadcast_mesh (MeshBase&) const // avoid spurious gcc warnings
{
  // no MPI, no need for this method...
  return;
}
#else 
// ------------------------------------------------------------
void MeshCommunication::broadcast_mesh (MeshBase& mesh) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;  

  // This function must be run on all processors at once
  parallel_only();

  START_LOG("broadcast_mesh()","MeshCommunication");

  // Explicitly clear the mesh on all but processor 0.
  if (libMesh::processor_id() != 0)
    mesh.clear();
  
  // Get important sizes
  unsigned int n_nodes      = mesh.n_nodes();
  unsigned int n_elem       = mesh.n_elem();
  unsigned int n_levels     = MeshTools::n_levels(mesh);
  unsigned int total_weight = MeshTools::total_weight(mesh);
  unsigned int dimension    = mesh.mesh_dimension();

  // Broadcast the sizes
  {
    std::vector<unsigned int> buf (4);
    
    buf[0] = n_nodes;
    buf[1] = n_elem;
    buf[2] = total_weight;
    buf[3] = dimension;
    
    // Broadcast
    Parallel::broadcast (buf);

    if (libMesh::processor_id() != 0)
      {
	n_nodes      = buf[0];
	n_elem       = buf[1];
	total_weight = buf[2];
        mesh.set_mesh_dimension(buf[3]);
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
	pts.reserve (LIBMESH_DIM*n_nodes);

	MeshBase::node_iterator       it     = mesh.nodes_begin();
	const MeshBase::node_iterator it_end = mesh.nodes_end();

	for (; it != it_end; ++it)
	  {
	    libmesh_assert (*it != NULL);
            libmesh_assert (!(*it)->valid_processor_id());
            libmesh_assert ((*it)->id()*LIBMESH_DIM == pts.size());
	    
	    const Point& p = **it;
	    
	    pts.push_back ( p(0) ); // x
#if LIBMESH_DIM > 1
	    pts.push_back ( p(1) ); // y
#endif
#if LIBMESH_DIM > 2
	    pts.push_back ( p(2) ); // z	  
#endif
	  }
      }
    else
      pts.resize (LIBMESH_DIM*n_nodes);

    // Sanity check for all processors
    libmesh_assert (pts.size() == (LIBMESH_DIM*n_nodes));
    
    // Broadcast the pts vector
    Parallel::broadcast (pts);

    // Add the nodes we just received if we are not
    // processor 0.
    if (libMesh::processor_id() != 0)
      {
	libmesh_assert (mesh.n_nodes() == 0);
	
	for (unsigned int i=0; i<pts.size(); i += LIBMESH_DIM)
	  mesh.add_point (Point(pts[i+0]
#if LIBMESH_DIM > 1
				, pts[i+1]
#endif
#if LIBMESH_DIM > 2
				, pts[i+2]
#endif
                               ),
			  i/LIBMESH_DIM);
      }
    
    libmesh_assert (mesh.n_nodes() == n_nodes);
  } // Done distributing the nodes

  
  // Now build up the elements vector which
  // contains the element types and connectivity
  {
    // The conn array contains the information needed to construct each element.
    // Pack all this information into one communication to avoid two latency hits
    // For each element it is of the form
    // [ level p_level r_flag p_flag etype subdomain_id 
    //   self_ID parent_ID which_child node_0 node_1 ... node_n]
    // We cannot use unsigned int because parent_ID can be negative
    std::vector<int> conn;

    // If we are processor 0, we must populate this vector and
    // broadcast it to the other processors.
    if (libMesh::processor_id() == 0)
      {
	conn.reserve (Elem::PackedElem::header_size*n_elem + total_weight);
	
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
                libmesh_assert (*it);
                libmesh_assert (!(*it)->valid_processor_id());
		const Elem* elem = *it;
		Elem::PackedElem::pack (conn, elem);
	      }
	  }
      }
    else
      conn.resize (Elem::PackedElem::header_size*n_elem + total_weight);
    
    // Sanity check for all processors
    libmesh_assert (conn.size() == (Elem::PackedElem::header_size*n_elem + total_weight));
    
    // Broadcast the element connectivity
    Parallel::broadcast (conn);

    // Build the elements we just received if we are not
    // processor 0.
    if (libMesh::processor_id() != 0)
      {
	libmesh_assert (mesh.n_elem() == 0);
	
	unsigned int cnt = 0;

        // This map keeps track of elements we've previously added to the mesh 
        // to avoid O(n) lookup times for parent pointers.
        std::map<unsigned int, Elem*> parents;

	while (cnt < conn.size())
	  {
	    // Unpack the element
	    Elem::PackedElem packed_elem (conn.begin()+cnt);

	    // Declare the element that we will add
            Elem* elem = NULL;

// 	    // Unpack the element header
// #ifdef LIBMESH_ENABLE_AMR
// 	    const int level             = conn[cnt++];
// 	    const int p_level           = conn[cnt++];
// 	    const Elem::RefinementState refinement_flag =
//               static_cast<Elem::RefinementState>(conn[cnt++]);
// 	    const Elem::RefinementState p_refinement_flag =
//               static_cast<Elem::RefinementState>(conn[cnt++]);
// #endif
//             const ElemType elem_type    = static_cast<ElemType>(conn[cnt++]);
// 	    const unsigned int elem_PID = conn[cnt++];
// 	    const int subdomain_ID      = conn[cnt++];
//             const int self_ID           = conn[cnt++];
// #ifdef LIBMESH_ENABLE_AMR
//             const int parent_ID         = conn[cnt++];
//             const int which_child       = conn[cnt++];

#ifdef LIBMESH_ENABLE_AMR

            if (packed_elem.parent_id() != -1) // Do a log(n) search for the parent
	      {
		Elem* my_parent = 
		  parents.count(packed_elem.parent_id()) ? parents[packed_elem.parent_id()] : NULL;
                
                // If the parent was not previously added, we cannot continue.
                if (my_parent == NULL)
		  {
		    libMesh::err << "Parent element with ID " << packed_elem.parent_id()
			          << " not found." << std::endl; 
		    libmesh_error();
		  }
		
                libmesh_assert (my_parent->refinement_flag() == Elem::INACTIVE);		
		
		elem = packed_elem.unpack (mesh, my_parent);
	      }
	    
            else // level 0 element has no parent
	      {
		libmesh_assert (packed_elem.level() == 0);
#endif 		
		elem = packed_elem.unpack (mesh);
#ifdef LIBMESH_ENABLE_AMR
	      }
#endif	    
            // Add elem to the map of parents, since it may have
            // children to be added later
            parents.insert(std::make_pair(elem->id(),elem));

	    cnt += Elem::PackedElem::header_size + packed_elem.n_nodes();
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
              libmesh_error();
          }

      } // end if iam != cpu 0
    
    
    libmesh_assert (mesh.n_elem() == n_elem);
  } // Done distributing the elements


  STOP_LOG("broadcast_mesh()","MeshCommunication");
}
#endif



#ifndef LIBMESH_HAVE_MPI
// ------------------------------------------------------------
void MeshCommunication::broadcast_bcs (const MeshBase&,
				       BoundaryInfo&) const // avoid spurious gcc warnings
{
}
#else
// ------------------------------------------------------------
void MeshCommunication::broadcast_bcs (const MeshBase& mesh,
				       BoundaryInfo& boundary_info) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (libMesh::n_processors() == 1)
    return;
  
  // This function must be run on all processors at once
  parallel_only();

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

    libmesh_assert (el_id.size() == side_id.size());
    libmesh_assert (el_id.size() == bc_id.size());
    
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
	      libmesh_assert (el_id[e] < mesh.n_elem());
	      
	      const Elem* elem = mesh.elem(el_id[e]);

	      libmesh_assert (elem != NULL);

	      // sanity: be sure that the element returned by mesh.elem() really has id()==el_id[e]
	      libmesh_assert(elem->id() == el_id[e]);

	      libmesh_assert (side_id[e] < elem->n_sides());
	    
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

    libmesh_assert (node_id.size() == bc_id.size());
    
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
	      libmesh_assert (node_id[n] < mesh.n_nodes());
	      
	      const Node* node = mesh.node_ptr (node_id[n]);

	      libmesh_assert (node != NULL);
	      
	      // sanity: be sure that the node returned by mesh.node_ptr() really has id()==node_id[n]
	      libmesh_assert(node->id() == node_id[n]);
	    
	      boundary_info.add_node (node, bc_id[n]);
	    }
      }
  }

  STOP_LOG("broadcast_bcs()","MeshCommunication");
}
#endif  



// ------------------------------------------------------------
void MeshCommunication::allgather (ParallelMesh& mesh) const
{
  START_LOG("allgather()","MeshCommunication");

  // The mesh should know it's about to be serialized
  libmesh_assert (mesh.is_serial());

  this->allgather_mesh (mesh);
  this->allgather_bcs  (mesh, *(mesh.boundary_info));

  // Inform new elements of their neighbors,
  // while resetting all remote_elem links
  mesh.find_neighbors(true);

  STOP_LOG("allgather()","MeshCommunication");
}

#ifndef LIBMESH_HAVE_MPI
  
// ------------------------------------------------------------
void MeshCommunication::allgather_mesh (ParallelMesh&) const
{
  // NO MPI == one processor, no need for this method
  return;
}
  
// ------------------------------------------------------------
void MeshCommunication::allgather_bcs (const ParallelMesh&,
				       BoundaryInfo&) const
{
  // NO MPI == one processor, no need for this method
  return;
}

#else

// ------------------------------------------------------------
void MeshCommunication::allgather_mesh (ParallelMesh& mesh) const
{
  // Check for quick return
  if (libMesh::n_processors() == 1)
    return;

  // This function must be run on all processors at once
  parallel_only();

  START_LOG ("allgather_mesh()","MeshCommunication");
  
  // Gather the number of nodes and elements on each processor.
  std::vector<unsigned int>
    n_nodes(libMesh::n_processors()), n_elem(libMesh::n_processors());
  
  {
    std::vector<unsigned int> n_objects(2);
    n_objects[0] = mesh.n_local_nodes();
    n_objects[1] = mesh.n_local_elem();

    Parallel::allgather(n_objects, /* identical_buffer_sizes = */ true);
    
    for (unsigned int p=0, idx=0; p<libMesh::n_processors(); p++)
      {
	n_nodes[p] = n_objects[idx++];
	n_elem[p]  = n_objects[idx++];	
      }

    libmesh_assert (mesh.n_local_nodes() == n_nodes[libMesh::processor_id()]);
    libmesh_assert (mesh.n_local_elem()  ==  n_elem[libMesh::processor_id()]);
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
	libmesh_assert (*it != NULL);
	    
	const Point &p = **it;

	xyz.push_back(p(0));
	xyz.push_back(p(1));
	xyz.push_back(p(2));
      }

    libmesh_assert (xyz.size() == 3*n_nodes[libMesh::processor_id()]);

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
	      libmesh_assert ((3*global_idx + 2) < xyz.size());
	      
	      Node *node = Node::build(xyz[3*global_idx + 0],
				       xyz[3*global_idx + 1],
				       xyz[3*global_idx + 2],
				       global_idx).release();
	      
	      libmesh_assert (node != NULL);
	      libmesh_assert (node->id() == global_idx);
	      
	      node->processor_id() = p;

	      mesh.insert_node(node);
	    }
	}
    
    // Check the result
    libmesh_assert (global_n_nodes == mesh.n_nodes());
  }
  
  //----------------------------------------------------
  // Gather the element connectivity from each processor.
  {
    // Get the sum of elem->n_nodes() for all local elements.  This
    // will allow for efficient preallocation.
    const unsigned int
      local_weight   = MeshTools::weight(mesh),
      local_n_levels = MeshTools::n_local_levels(mesh);

    unsigned int global_n_levels = local_n_levels;
    Parallel::max (global_n_levels);
    
    // The conn array contains the information needed to construct each element.
    std::vector<int> conn; conn.reserve 
      (Elem::PackedElem::header_size*n_elem[libMesh::processor_id()] + local_weight);
						
    for (unsigned int level=0; level<=local_n_levels; level++)
      {
	ParallelMesh::element_iterator        it  = mesh.local_level_elements_begin(level);
	const ParallelMesh::element_iterator  end = mesh.local_level_elements_end(level);

	for (; it != end; ++it)
	  {
	    const Elem* elem = *it;

	    libmesh_assert (elem != NULL);
	    libmesh_assert (elem->level() == level);
	    
	    // We're not handling unpartitioned elements!
	    libmesh_assert (elem->processor_id() != DofObject::invalid_processor_id);

	    // Only local elements!
	    libmesh_assert (elem->processor_id() == libMesh::processor_id());

	    Elem::PackedElem::pack (conn, elem);	    
	  }
      } // ...that was easy.

    libmesh_assert (conn.size() == 
      Elem::PackedElem::header_size*n_elem[libMesh::processor_id()] + local_weight);

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
	if (p != libMesh::processor_id()) // We've already got our
          {                               // own local elements!
	    unsigned int cnt = conn_offset[p]; // counter into the conn[] array.
	    
#ifndef NDEBUG
	    // The first and last global indices are only used for error-checking.
	    // Avoid declaring them unless asserts are enabled.
	    const unsigned int
	      first_global_idx = elem_offsets[p],
	      last_global_idx  = first_global_idx + n_elem[p];
#endif
	    
	    // Process each element for processor p.
	    // Note this must work in the case when conn_size[p] == 0.
	    while (cnt < (conn_offset[p] + conn_size[p]))
	      {
		// Unpack the element
		Elem::PackedElem packed_elem (conn.begin()+cnt);

 		// We require contiguous numbering on each processor
 		// for elements.
 		libmesh_assert (packed_elem.id() >= first_global_idx);
 		libmesh_assert (packed_elem.id()  < last_global_idx);

#ifdef LIBMESH_ENABLE_AMR
		libmesh_assert ((packed_elem.level() == 0) || (packed_elem.parent_id() != -1));
  
		// Ignore elements not matching the current level.
		if (packed_elem.level() > level) // skip all entries in the conn array for this element.
		  cnt += Elem::PackedElem::header_size + packed_elem.n_nodes();

                else
#endif
		if (packed_elem.level() < level ||     // we should already have
		    (packed_elem.level() == level &&   // lower level and some
                     mesh.elem(packed_elem.id())))     // ghost elements
		  {
                    // No need to construct a dummy element of this type
#ifndef NDEBUG
		    // The elem information need not be extracted unless asserts are enabled.
		    const Elem* elem = mesh.elem(packed_elem.id());
#endif
                    libmesh_assert (elem);
#ifdef LIBMESH_ENABLE_AMR
		    libmesh_assert (!elem->parent() ||
				    elem->parent()->dim() != elem->dim() ||
				    elem->parent()->id() ==
				    static_cast<unsigned int>(packed_elem.parent_id()));
		    libmesh_assert (elem->p_level()           == packed_elem.p_level());
		    libmesh_assert (elem->refinement_flag()   == packed_elem.refinement_flag());
		    libmesh_assert (elem->p_refinement_flag() == packed_elem.p_refinement_flag());
#endif
		    libmesh_assert (elem->type()              == packed_elem.type());
		    libmesh_assert (elem->subdomain_id()      == packed_elem.subdomain_id());
		    libmesh_assert (elem->id()                == packed_elem.id());
		    libmesh_assert (elem->n_nodes()           == packed_elem.n_nodes());

		    cnt += Elem::PackedElem::header_size + packed_elem.n_nodes();
		  }
		// Those are the easy cases...
		// now elem_level == level and we don't have it
		else
		  {
		    // Declare the element we will add
		    Elem* elem = NULL;

#ifdef LIBMESH_ENABLE_AMR
		    // Maybe find its parent
		    if (level > 0)
		      {
			Elem* my_parent = mesh.elem(packed_elem.parent_id());

			// If the parent was not previously added, we
			// cannot continue.
			if (my_parent == NULL)
			  {
			    libMesh::err << "Parent element with ID " << packed_elem.parent_id()
				          << " not found." << std::endl; 
			    libmesh_error();
			  }
			
			elem = packed_elem.unpack (mesh, my_parent);
		      }
		    else
                      {
                        libmesh_assert (packed_elem.parent_id() == -1);
#endif // LIBMESH_ENABLE_AMR
		        elem = packed_elem.unpack (mesh);
#ifdef LIBMESH_ENABLE_AMR
                      }
#endif

		    // Good to go.  Add to the mesh.
		    libmesh_assert (elem);
		    mesh.insert_elem(elem);

		    libmesh_assert (elem->n_nodes() == packed_elem.n_nodes());

		    cnt += Elem::PackedElem::header_size + packed_elem.n_nodes();
		    
		  } // end elem_level == level
	      }
	    
	  }   

    // Check the result
    libmesh_assert (global_n_elem == mesh.n_elem());
  }

  // All done!
  STOP_LOG ("allgather_mesh()","MeshCommunication");
}



// ------------------------------------------------------------
void MeshCommunication::allgather_bcs (const ParallelMesh& mesh,
				       BoundaryInfo& boundary_info) const
{
  // Check for quick return
  if (libMesh::n_processors() == 1)
    return;

  // This function must be run on all processors at once
  parallel_only();

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

    libmesh_assert (el_id.size() == side_id.size());
    libmesh_assert (el_id.size() == bc_id.size());
    
    const unsigned int n_bcs = el_id.size();

    // reserve an upper bound for the number of BCs
    xfer_elem_bcs.reserve(3*n_bcs);

    // populate the xfer_elem_bcs list with *local* elements only.
    for (unsigned int bc=0; bc<n_bcs; bc++)
      {
	const Elem* elem = mesh.elem(el_id[bc]);
	
	// sanity: be sure that the element returned by mesh.elem() really has id()==el_id[e]
	libmesh_assert(elem != NULL);
	libmesh_assert(elem->id() == el_id[bc]);
	libmesh_assert(elem->level() == 0);
	libmesh_assert(side_id[bc] < elem->n_sides());

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

    libmesh_assert (node_id.size() == bc_id.size());

    const unsigned int n_bcs = node_id.size();

    // reserve an upper bound for the number of BCs
    xfer_node_bcs.reserve(2*n_bcs);

    // populate the xfer_node_bcs witl *local* nodes only
    for (unsigned int bc=0; bc<n_bcs; bc++)
      {
	const Node* node = mesh.node_ptr(node_id[bc]);
	
	libmesh_assert(node != NULL);
	libmesh_assert(node->id() == node_id[bc]);

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
  // Insert the elements
  {
    Parallel::allgather (xfer_elem_bcs);

    const unsigned int n_bcs = xfer_elem_bcs.size()/3;

    for (unsigned int bc=0, cnt=0; bc<n_bcs; bc++)
      {
	const Elem* elem = mesh.elem(xfer_elem_bcs[cnt++]);
	const unsigned short int side_id = xfer_elem_bcs[cnt++];
	const short int bc_id = xfer_elem_bcs[cnt++];

	boundary_info.add_side (elem, side_id, bc_id);	
      }

    // no need for this any more
    Utility::deallocate (xfer_elem_bcs);
  }

  
  // Insert the nodes
  {
    Parallel::allgather (xfer_node_bcs);

    const unsigned int n_bcs = xfer_node_bcs.size()/2;

    for (unsigned int bc=0, cnt=0; bc<n_bcs; bc++)
      {
	const Node* node = mesh.node_ptr (xfer_node_bcs[cnt++]);
	const short int bc_id = xfer_node_bcs[cnt++];

	boundary_info.add_node (node, bc_id);
      }
 
    // no need for this any more
    Utility::deallocate (xfer_node_bcs);
  }

 
#ifndef NDEBUG
  
  // Make sure all processors agree on the number of boundary ids.
  const unsigned int n_bc_ids = boundary_info.n_boundary_ids();
  unsigned int global_n_bc_ids = n_bc_ids;
  
  Parallel::max (global_n_bc_ids);
  libmesh_assert (n_bc_ids == global_n_bc_ids);
  
#endif
  

  STOP_LOG  ("allgather_bcs()","MeshCommunication");
}


#endif // LIBMESH_HAVE_MPI



// Functor for make_elems_parallel_consistent and
// make_node_ids_parallel_consistent
namespace {

struct SyncIds
{
typedef unsigned int datum;
typedef void (MeshBase::*renumber_obj)(unsigned int, unsigned int);

SyncIds(MeshBase &_mesh, renumber_obj _renumberer) :
  mesh(_mesh),
  renumber(_renumberer) {}

MeshBase &mesh;
renumber_obj renumber;
// renumber_obj &renumber;

// Find the id of each requested DofObject -
// sync_dofobject_data_by_xyz already did the work for us
void gather_data (const std::vector<unsigned int>& ids,
                  std::vector<datum>& ids_out)
{
  ids_out = ids;
}

void act_on_data (const std::vector<unsigned int>& old_ids,
                  std::vector<datum>& new_ids)
{
  for (unsigned int i=0; i != old_ids.size(); ++i)
    if (old_ids[i] != new_ids[i])
      (mesh.*renumber)(old_ids[i], new_ids[i]);
}
};
}



// ------------------------------------------------------------
void MeshCommunication::make_node_ids_parallel_consistent
  (MeshBase &mesh,
   LocationMap<Node> &loc_map)
{
  // This function must be run on all processors at once
  parallel_only();

  START_LOG ("make_node_ids_parallel_consistent()", "MeshCommunication");

  SyncIds syncids(mesh, &MeshBase::renumber_node);
  Parallel::sync_dofobject_data_by_xyz
    (mesh.nodes_begin(), mesh.nodes_end(),
     loc_map, syncids);

  STOP_LOG ("make_node_ids_parallel_consistent()", "MeshCommunication");
}



// ------------------------------------------------------------
void MeshCommunication::make_elems_parallel_consistent(MeshBase &mesh)
{
  // This function must be run on all processors at once
  parallel_only();

  START_LOG ("make_elems_parallel_consistent()", "MeshCommunication");

  SyncIds syncids(mesh, &MeshBase::renumber_elem);
  Parallel::sync_element_data_by_parent_id
    (mesh, mesh.active_elements_begin(),
     mesh.active_elements_end(), syncids);

  STOP_LOG ("make_elems_parallel_consistent()", "MeshCommunication");
}



// Functors for make_node_proc_ids_parallel_consistent
namespace {

struct SyncProcIds
{
typedef unsigned int datum;

SyncProcIds(MeshBase &_mesh) : mesh(_mesh) {}

MeshBase &mesh;

// ------------------------------------------------------------
void gather_data (const std::vector<unsigned int>& ids,
                  std::vector<datum>& data)
{
  // Find the processor id of each requested node
  data.resize(ids.size());

  for (unsigned int i=0; i != ids.size(); ++i)
    {
      // Look for this point in the mesh
      Node *node = mesh.node_ptr(ids[i]);

      // We'd better find every node we're asked for
      libmesh_assert (node);

      // Return the node's correct processor id,
      data[i] = node->processor_id();
    }
}

// ------------------------------------------------------------
void act_on_data (const std::vector<unsigned int>& ids,
                  std::vector<datum> proc_ids)
{
  // Set the ghost node processor ids we've now been informed of
  for (unsigned int i=0; i != ids.size(); ++i)
    {
      Node *node = mesh.node_ptr(ids[i]);
      node->processor_id() = proc_ids[i];
    }
}
};
}



// ------------------------------------------------------------
void MeshCommunication::make_node_proc_ids_parallel_consistent
  (MeshBase& mesh,
   LocationMap<Node>& loc_map)
{
  START_LOG ("make_node_proc_ids_parallel_consistent()", "MeshCommunication");

  // This function must be run on all processors at once
  parallel_only();

  // When this function is called, each section of a parallelized mesh
  // should be in the following state:
  //
  // All nodes should have the exact same physical location on every
  // processor where they exist.
  //
  // Local nodes should have unique authoritative ids,
  // and processor ids consistent with all processors which own
  // an element touching them.
  //
  // Ghost nodes touching local elements should have processor ids
  // consistent with all processors which own an element touching
  // them.
  
  SyncProcIds sync(mesh);
  Parallel::sync_dofobject_data_by_xyz
    (mesh.nodes_begin(), mesh.nodes_end(), loc_map, sync);

  STOP_LOG ("make_node_proc_ids_parallel_consistent()", "MeshCommunication");
}



// ------------------------------------------------------------
void MeshCommunication::make_nodes_parallel_consistent
  (MeshBase &mesh,
   LocationMap<Node> &loc_map)
{
  // This function must be run on all processors at once
  parallel_only();

  // Create the loc_map if it hasn't been done already
  bool need_map_update = (mesh.nodes_begin() != mesh.nodes_end() && 
                          loc_map.empty());
  Parallel::max(need_map_update);

  if (need_map_update)
    loc_map.init(mesh);

  // When this function is called, each section of a parallelized mesh
  // should be in the following state:
  //
  // All nodes should have the exact same physical location on every
  // processor where they exist.
  //
  // Local nodes should have unique authoritative ids,
  // and processor ids consistent with all processors which own
  // an element touching them.
  //
  // Ghost nodes touching local elements should have processor ids
  // consistent with all processors which own an element touching
  // them.
  //
  // Ghost nodes should have ids which are either already correct
  // or which are in the "unpartitioned" id space.

  // First, let's sync up processor ids.  Some of these processor ids
  // may be "wrong" from coarsening, but they're right in the sense
  // that they'll tell us who has the authoritative dofobject ids for
  // each node.
  this->make_node_proc_ids_parallel_consistent(mesh, loc_map);

  // Second, sync up dofobject ids.
  this->make_node_ids_parallel_consistent(mesh, loc_map);

  // Finally, correct the processor ids to make DofMap happy
  MeshTools::correct_node_proc_ids(mesh, loc_map);
}



// ------------------------------------------------------------
void MeshCommunication::delete_remote_elements(ParallelMesh& mesh, const std::set<Elem *> & extra_ghost_elem_ids) const
{
  // The mesh should know it's about to be parallelized
  libmesh_assert (!mesh.is_serial());

  START_LOG("delete_remote_elements()", "MeshCommunication");

#ifdef DEBUG
  // We expect maximum ids to be in sync so we can use them to size
  // vectors
  Parallel::verify(mesh.max_node_id());
  Parallel::verify(mesh.max_elem_id());
  libmesh_assert(mesh.parallel_max_node_id() == mesh.max_node_id());
  libmesh_assert(mesh.parallel_max_elem_id() == mesh.max_elem_id());
#endif

  // FIXME - should these be "unsorted_set"s?  O(N) is O(N)...
  std::vector<bool> local_nodes(mesh.max_node_id(), false);
  std::vector<bool> semilocal_elems(mesh.max_elem_id(), false);

  // We don't want to delete any element that shares a node
  // with or is an ancestor of a local element.
  MeshBase::const_element_iterator l_elem_it = mesh.local_elements_begin(),
                                   l_end     = mesh.local_elements_end();
  for (; l_elem_it != l_end; ++l_elem_it)
    {
      const Elem *elem = *l_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        {
          unsigned int nodeid = elem->node(n);
          libmesh_assert(nodeid < local_nodes.size());
          local_nodes[nodeid] = true;
        }
      while (elem)
        {
          unsigned int elemid = elem->id();
          libmesh_assert(elemid < semilocal_elems.size());
          semilocal_elems[elemid] = true;

          const Elem *parent = elem->parent();
	  // Don't proceed from a boundary mesh to an interior mesh
          if (parent && parent->dim() != elem->dim())
            break;

          elem = parent;
        }
    }

  // We don't want to delete any element that shares a node
  // with or is an ancestor of an unpartitioned element either.
  MeshBase::const_element_iterator
    u_elem_it = mesh.unpartitioned_elements_begin(),
    u_end     = mesh.unpartitioned_elements_end();
  
  for (; u_elem_it != u_end; ++u_elem_it)
    {
      const Elem *elem = *u_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        local_nodes[elem->node(n)] = true;
      while (elem)
        {
          semilocal_elems[elem->id()] = true;

          const Elem *parent = elem->parent();
	  // Don't proceed from a boundary mesh to an interior mesh
          if (parent && parent->dim() != elem->dim())
            break;

          elem = parent;
        }
    }

  // Flag all the elements that share nodes with
  // local and unpartitioned elements, along with their ancestors
  MeshBase::element_iterator nl_elem_it = mesh.not_local_elements_begin(),
                             nl_end     = mesh.not_local_elements_end();
  for (; nl_elem_it != nl_end; ++nl_elem_it)
    {
      const Elem *elem = *nl_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        if (local_nodes[elem->node(n)])
          {
            while (elem)
              {
                semilocal_elems[elem->id()] = true;

                const Elem *parent = elem->parent();
	        // Don't proceed from a boundary mesh to an interior mesh
                if (parent && parent->dim() != elem->dim())
                  break;

                elem = parent;
              }
            break;
          }
    }

  // Don't delete elements that we were explicitly told not to
  for(std::set<Elem *>::iterator it = extra_ghost_elem_ids.begin();
      it != extra_ghost_elem_ids.end();
      ++it)
    semilocal_elems[(*it)->id()] = true;

  // Delete all the elements we have no reason to save,
  // starting with the most refined so that the mesh
  // is valid at all intermediate steps
  unsigned int n_levels = MeshTools::n_levels(mesh);

  for (int l = n_levels - 1; l >= 0; --l)
    {
      MeshBase::element_iterator lev_elem_it = mesh.level_elements_begin(l),
                                 lev_end     = mesh.level_elements_end(l);
      for (; lev_elem_it != lev_end; ++lev_elem_it)
        {
          Elem *elem = *lev_elem_it;
          libmesh_assert (elem);
          // Make sure we don't leave any invalid pointers
          if (!semilocal_elems[elem->id()])
            elem->make_links_to_me_remote();

          // Subactive neighbor pointers aren't preservable here
          if (elem->subactive())
            for (unsigned int s=0; s != elem->n_sides(); ++s)
              elem->set_neighbor(s, NULL);

          // delete_elem doesn't currently invalidate element
          // iterators... that had better not change
          if (!semilocal_elems[elem->id()])
            mesh.delete_elem(elem);
        }
    }

#ifdef DEBUG
  MeshTools::libmesh_assert_valid_refinement_tree(mesh);
#endif

  STOP_LOG("delete_remote_elements()", "MeshCommunication");

  // Now make sure the containers actually shrink - strip
  // any newly-created NULL voids out of the element array
  mesh.renumber_nodes_and_elements();
}




// // Pack all this information into one communication to avoid two latency hits
// // For each element it is of the form
// // [ level p_level r_flag p_flag etype subdomain_id 
// //   self_ID parent_ID which_child node_0 node_1 ... node_n]
// // We cannot use unsigned int because parent_ID can be negative
// void MeshCommunication::pack_element (std::vector<int> &conn, const Elem* elem) const
// {
//   libmesh_assert (elem != NULL);

// #ifdef LIBMESH_ENABLE_AMR
//   conn.push_back (static_cast<int>(elem->level()));
//   conn.push_back (static_cast<int>(elem->p_level()));
//   conn.push_back (static_cast<int>(elem->refinement_flag()));
//   conn.push_back (static_cast<int>(elem->p_refinement_flag()));
// #endif
//   conn.push_back (static_cast<int>(elem->type()));
//   conn.push_back (static_cast<int>(elem->processor_id()));
//   conn.push_back (static_cast<int>(elem->subdomain_id()));
//   conn.push_back (elem->id());
		
// #ifdef LIBMESH_ENABLE_AMR
//   // use parent_ID of -1 to indicate a level 0 element
//   if (elem->level() == 0)
//     {
//       conn.push_back(-1);
//       conn.push_back(-1);
//     }
//   else
//     {
//       conn.push_back(elem->parent()->id());
//       conn.push_back(elem->parent()->which_child_am_i(elem));
//     }
// #endif
  
//   for (unsigned int n=0; n<elem->n_nodes(); n++)
//     conn.push_back (elem->node(n));		
// }

} // namespace libMesh
