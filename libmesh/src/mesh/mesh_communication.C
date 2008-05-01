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



//-----------------------------------------------
// anonymous namespace for implementation details
namespace {

#ifdef ENABLE_AMR
  const unsigned int packed_elem_header_size = 10;
#else
  const unsigned int packed_elem_header_size = 4;
#endif


#ifdef HAVE_MPI
  /**
   * Convenient way to communicate nodes.
   */
  struct PackedNode
  {
    unsigned int id;
    unsigned int pid;
    Real x;
    Real y;
    Real z;

    PackedNode () :
      id(0),
      x(0.),
      y(0.),
      z(0.)
    {}

    PackedNode (const Node &node) :
      id(node.id()),
      pid(node.processor_id()),
      x(node(0)),
      y(node(1)),
      z(node(2))
    {}

    AutoPtr<Node> build_node () const
    {
      AutoPtr<Node> node(new Node(x,y,z,id));
      node->processor_id() = pid;
      return node;
    }

    Point build_point () const
    {
      return Point(x,y,z);
    }
    
    static MPI_Datatype create_mpi_datatype ();
    
  };
  
  MPI_Datatype PackedNode::create_mpi_datatype ()
  {
    MPI_Datatype packed_node_type;
    MPI_Datatype types[] = { MPI_UNSIGNED, MPI_REAL };
    int blocklengths[] = { 2, 3 };
    MPI_Aint displs[2];
    
    // create a Packed node and get the addresses of the elements.
    // this will properly handle id/pid getting padded, for example,
    // in which case id and x may not be 2*sizeof(unsigned int) apart.
    PackedNode pn;

    MPI_Address (&pn.id, &displs[0]);
    MPI_Address (&pn.x,  &displs[1]);
    displs[1] -= displs[0];
    displs[0] = 0;

#if MPI_VERSION > 1
      MPI_Type_create_struct (2, blocklengths, displs, types, &packed_node_type);
#else
      MPI_Type_struct (2, blocklengths, displs, types, &packed_node_type);
#endif

      return packed_node_type;
    }
#endif

  /**
   * Specific weak ordering for Elem*'s to be used in a set.
   * We use the id, but first sort by level.  This guarantees 
   * when traversing the set from beginning to end the lower 
   * level (parent) elements are encountered first. Additionally,
   * order siblings so that child(0) is added before child(1).
   */
  struct CompareElemIdsByLevel
  {
    bool operator()(const Elem *a,
		    const Elem *b) const
    {
      libmesh_assert (a);
      libmesh_assert (b);
      
      if (a->level() == b->level())
	{
// 	  // order siblings sequentially
// 	  if (a->parent() && b->parent())   // a & b have parents
// 	    if (a->parent() == b->parent()) // which are the same
// 	      return (a->parent()->which_child_am_i(a) <
// 		      b->parent()->which_child_am_i(b));

	  // otherwise order by id
	  return a->id() < b->id();
	}

      // otherwise order by level
      return a->level() < b->level();
    }
  };    
}



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



#ifndef HAVE_MPI // avoid spurious gcc warnings
void MeshCommunication::redistribute (ParallelMesh &) const
{
  // no MPI, no redistribution
  return;
}
#else
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
  MPI_Datatype packed_node_datatype = PackedNode::create_mpi_datatype();
  MPI_Type_commit (&packed_node_datatype);

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
  
  std::vector<std::vector<PackedNode> >
    nodes_sent(libMesh::n_processors());

  std::vector<std::vector<int> > 
    elements_sent(libMesh::n_processors()),
    node_bcs_sent(libMesh::n_processors()),
    element_bcs_sent(libMesh::n_processors());
    
  std::vector<Parallel::request>
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
#ifdef ENABLE_AMR
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
#ifdef ENABLE_AMR
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
		nodes_sent[pid].push_back(PackedNode(**node_it));

		// add the node if it has BCs
		if (mesh.boundary_info->boundary_id(*node_it) !=
		    mesh.boundary_info->invalid_id)
		  {
		    node_bcs_sent[pid].push_back((*node_it)->id());
		    node_bcs_sent[pid].push_back(mesh.boundary_info->boundary_id(*node_it));
		  }
	      }

	    // the size of the node bc buffer we will ship to pid
	    send_n_nodes_and_elem_per_proc[5*pid+3] = node_bcs_sent[pid].size();	    
	    
	    // send the nodes off to the destination processor
	    node_send_requests.push_back(Parallel::request());
	  
	    Parallel::isend (pid,
			     nodes_sent[pid],
			     packed_node_datatype,
			     node_send_requests.back(),
			     /* tag = */ 0);

	    if (!node_bcs_sent[pid].empty())
	      {
		node_bc_requests.push_back(Parallel::request());

		Parallel::isend (pid,
				 node_bcs_sent[pid],
				 node_bc_requests.back(),
				 /* tag = */ 2);
	      }
	  }
	
	if (!elements_to_send.empty())
	  {
	    // the number of elements we will send to this processor
	    send_n_nodes_and_elem_per_proc[5*pid+1] = elements_to_send.size();

	    for (std::set<const Elem*, CompareElemIdsByLevel>::const_iterator 
		   elem_it = elements_to_send.begin(); elem_it != elements_to_send.end(); ++elem_it)
	      {
		pack_element (elements_sent[pid], *elem_it);

		// if this is a level-0 element look for boundary conditions
	        if ((*elem_it)->level() == 0)
		  for (unsigned int s=0; s<(*elem_it)->n_sides(); s++)
		    if ((*elem_it)->neighbor(s) == NULL)
		      if (mesh.boundary_info->boundary_id (*elem_it, s) !=
			  mesh.boundary_info->invalid_id)
			{
			  element_bcs_sent[pid].push_back ((*elem_it)->id());
			  element_bcs_sent[pid].push_back (s);
			  element_bcs_sent[pid].push_back (mesh.boundary_info->boundary_id (*elem_it, s));
			}
	      }
	    
	    // the packed connectivity size to send to this processor
	    send_n_nodes_and_elem_per_proc[5*pid+2] = elements_sent[pid].size();

	    // send the elements off to the destination processor
	    element_send_requests.push_back(Parallel::request());
	  
	    Parallel::isend (pid,
			     elements_sent[pid],
			     element_send_requests.back(),
			     /* tag = */ 1);

	    // the size of the element bc buffer we will ship to pid
	    send_n_nodes_and_elem_per_proc[5*pid+4] = element_bcs_sent[pid].size();

	    if (!element_bcs_sent[pid].empty())
	      {
		element_bc_requests.push_back(Parallel::request());

		Parallel::isend (pid,
				 element_bcs_sent[pid],
				 element_bc_requests.back(),
				 /* tag = */ 3);
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
    n_send_node_pairs=0,    n_send_elem_pairs=0,
    n_send_node_bc_pairs=0, n_send_elem_bc_pairs=0, 
    n_recv_node_pairs=0,    n_recv_elem_pairs=0,
    n_recv_node_bc_pairs=0, n_recv_elem_bc_pairs=0, 
    max_n_nodes_received=0, max_conn_size_received=0,
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
    std::vector<PackedNode> received_nodes(max_n_nodes_received);

    // We now know how many processors will be sending us nodes.
    for (unsigned int node_comm_step=0; node_comm_step<n_recv_node_pairs; node_comm_step++)
      {
	// but we don't necessarily want to impose an ordering, so
	// just grab whatever message is available next.
	Parallel::Status status =
	  Parallel::recv (Parallel::any_source,
			  received_nodes,
			  packed_node_datatype,
			  /* tag = */ 0);
	const unsigned int source_pid = status.source();
	const unsigned int n_nodes_received =
	  recv_n_nodes_and_elem_per_proc[5*source_pid+0];
	libmesh_assert (n_nodes_received <= received_nodes.size());
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
	  Parallel::recv (Parallel::any_source,
			  received_node_bcs,
			  /* tag = */ 2);
	const unsigned int source_pid = status.source();
	libmesh_assert (status.size() == recv_n_nodes_and_elem_per_proc[5*source_pid+3]);

	const unsigned int buffer_size = status.size();
	unsigned int cnt=0;

	while (cnt < buffer_size)
	  {
	    const unsigned int node_id = received_node_bcs[cnt++];
	    const int bc_id            = received_node_bcs[cnt++];

	    libmesh_assert (mesh.node_ptr(node_id));
	    mesh.boundary_info->add_node (mesh.node_ptr(node_id), bc_id);
	  }	
      }
  }


  
  // Receive elements.  Size this array for the largest message.
  {
    std::vector<unsigned int> received_elements(max_conn_size_received);

    // Similarly we know how many processors are sending us elements, 
    // but we don't really care in what order we receive them.
    for (unsigned int elem_comm_step=0; elem_comm_step<n_recv_elem_pairs; elem_comm_step++)
      {
	Parallel::Status status =
	  Parallel::recv (Parallel::any_source,
			  received_elements,
			  /* tag = */ 1);
	const unsigned int source_pid = status.source();
	const unsigned int n_elem_received = 
	  recv_n_nodes_and_elem_per_proc[5*source_pid+1];
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
	    // Unpack the element header
#ifdef ENABLE_AMR
	    const unsigned int elem_level  = received_elements[cnt++];
	    const unsigned int p_level     = received_elements[cnt++];
	    const Elem::RefinementState refinement_flag =
	      static_cast<Elem::RefinementState>(received_elements[cnt++]);
	    const Elem::RefinementState p_refinement_flag =
	      static_cast<Elem::RefinementState>(received_elements[cnt++]);
#else
	    const unsigned int elem_level = 0;
#endif
	    const ElemType elem_type       = static_cast<ElemType>(received_elements[cnt++]);
	    const unsigned int elem_PID    = received_elements[cnt++];
	    const int subdomain_ID         = received_elements[cnt++];
	    const unsigned int self_ID     = received_elements[cnt++];
#ifdef ENABLE_AMR
	    const int parent_ID            = received_elements[cnt++];
	    const unsigned int which_child = received_elements[cnt++];
#endif
	    
	    // The ParallelMesh::elem(i) member will return NULL if the element
	    // is not in the mesh.  We rely on that here, so it better not change!
	    Elem *elem = mesh.elem(self_ID);
	    
	    // if we already have this element, make sure its properties match
	    // but then go on
	    if (elem)
	      {
		libmesh_assert (elem->level() == elem_level);
		libmesh_assert (elem->id() == self_ID);
		libmesh_assert (elem->processor_id() == elem_PID);
		libmesh_assert (elem->subdomain_id() == subdomain_ID);
		libmesh_assert (elem->type() == elem_type);
#ifdef ENABLE_AMR
		libmesh_assert (elem->p_level() == p_level);
		libmesh_assert (elem->refinement_flag() == refinement_flag);
		libmesh_assert (elem->p_refinement_flag() == p_refinement_flag);
		
		if (elem->level() > 0)
		  {
		    libmesh_assert (elem->parent()->id() == static_cast<unsigned int>(parent_ID));
		    libmesh_assert (elem->parent()->child(which_child) == elem);
		  }	      
#endif
		libmesh_assert (elem->n_nodes() == Elem::type_to_n_nodes_map[elem_type]);
	      
		// skip the connectivity for this element
		cnt += elem->n_nodes();
	      }
	    else
	      {
		// We need to add the element.
#ifdef ENABLE_AMR
		// maybe find the parent
		if (elem_level > 0)
		  {
		    Elem *parent = mesh.elem(parent_ID);
		    
		    // Note that we were very careful to construct the send connectivity
		    // so that parents are encountered before children.  So if we get here
		    // and can't find the parent that is a fatal error.
		    if (parent == NULL)
		      {
			std::cerr << "Parent element with ID " << parent_ID 
				  << " not found." << std::endl; 
			libmesh_error();
		      }
		    
		    elem = Elem::build(elem_type, parent).release();
		    libmesh_assert (elem);
		    parent->add_child(elem, which_child);
		    libmesh_assert (parent->type() == elem->type());
		    libmesh_assert (parent->child(which_child) == elem);
		  }
		else
		  {
		    libmesh_assert (parent_ID == -1);
#endif // ENABLE_AMR
		    elem = Elem::build(elem_type).release();
		    libmesh_assert (elem);
#ifdef ENABLE_AMR
		  }
		
		// Assign the IDs
		elem->set_p_level(p_level);
		elem->set_refinement_flag(refinement_flag);
		elem->set_p_refinement_flag(p_refinement_flag);
		libmesh_assert (elem->level() == static_cast<unsigned int>(elem_level));
#endif
		elem->processor_id() = elem_PID;
		elem->subdomain_id() = subdomain_ID;
		elem->set_id()       = self_ID;
		
		// Assign the connectivity
		for (unsigned int n=0; n<elem->n_nodes(); n++)
		  {
		    libmesh_assert (cnt < received_elements.size());		  
		    elem->set_node(n) = mesh.node_ptr (received_elements[cnt++]);
		  }
		
		// Good to go.  Add to the mesh.
		mesh.insert_elem(elem);
	      }
	    
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
	  Parallel::recv (Parallel::any_source,
			  received_element_bcs,
			  /* tag = */ 3);
	const unsigned int source_pid = status.source();
	libmesh_assert (status.size() == recv_n_nodes_and_elem_per_proc[5*source_pid+4]);

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
  MPI_Type_free (&packed_node_datatype);
  
  STOP_LOG("redistribute()","MeshCommunication");  
}
#endif // HAVE_MPI




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
	    libmesh_assert (*it != NULL);
            libmesh_assert (!(*it)->valid_processor_id());
            libmesh_assert ((*it)->id()*3 == pts.size());
	    
	    const Point& p = **it;
	    
	    pts.push_back ( p(0) ); // x
	    pts.push_back ( p(1) ); // y
	    pts.push_back ( p(2) ); // z	  
	  }
      }
    else
      pts.resize (3*n_nodes);

    // Sanity check for all processors
    libmesh_assert (pts.size() == (3*n_nodes));
    
    // Broadcast the pts vector
    Parallel::broadcast (pts);

    // Add the nodes we just received if we are not
    // processor 0.
    if (libMesh::processor_id() != 0)
      {
	libmesh_assert (mesh.n_nodes() == 0);
	
	for (unsigned int i=0; i<pts.size(); i += 3)
	  mesh.add_point (Point(pts[i+0],
				pts[i+1],
				pts[i+2]),
			  i/3);
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
	conn.reserve (packed_elem_header_size*n_elem + total_weight);
	
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
		pack_element (conn, elem);
	      }
	  }
      }
    else
      conn.resize (packed_elem_header_size*n_elem + total_weight);
    
    // Sanity check for all processors
    libmesh_assert (conn.size() == (packed_elem_header_size*n_elem + total_weight));
    
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
	    // Declare the element that we will add
            Elem* elem = NULL;

	    // Unpack the element header
#ifdef ENABLE_AMR
	    const int level             = conn[cnt++];
	    const int p_level           = conn[cnt++];
	    const Elem::RefinementState refinement_flag =
              static_cast<Elem::RefinementState>(conn[cnt++]);
	    const Elem::RefinementState p_refinement_flag =
              static_cast<Elem::RefinementState>(conn[cnt++]);
#endif
            const ElemType elem_type    = static_cast<ElemType>(conn[cnt++]);
	    const unsigned int elem_PID = conn[cnt++];
	    const int subdomain_ID      = conn[cnt++];
            const int self_ID           = conn[cnt++];
#ifdef ENABLE_AMR
            const int parent_ID         = conn[cnt++];
            const int which_child       = conn[cnt++];

            if (parent_ID != -1) // Do a log(n) search for the parent
	      {
		Elem* my_parent = parents.count(parent_ID) ? parents[parent_ID] : NULL;
                
                // If the parent was not previously added, we cannot continue.
                if (my_parent == NULL)
		  {
		    std::cerr << "Parent element with ID " << parent_ID 
			      << " not found." << std::endl; 
		    libmesh_error();
		  }
		
                libmesh_assert (my_parent->refinement_flag() == Elem::INACTIVE);
		
		elem = Elem::build(elem_type,my_parent).release();
		my_parent->add_child(elem);
		libmesh_assert (my_parent->type() == elem->type());
                libmesh_assert (my_parent->child(which_child) == elem);
	      }
	    
            else // level 0 element has no parent
	      {
		libmesh_assert (level == 0);
#endif 
		
		// should be able to just use the integer elem_type
		elem = Elem::build(elem_type).release();
#ifdef ENABLE_AMR
	      }

	    // Assign the IDs
	    libmesh_assert (elem->level() == static_cast<unsigned int>(level));
	    elem->set_refinement_flag(refinement_flag); 
	    elem->set_p_refinement_flag(p_refinement_flag); 
	    elem->set_p_level(p_level); 
#endif
	    elem->processor_id() = elem_PID;
            elem->subdomain_id() = subdomain_ID;
            elem->set_id() = self_ID;
	    
            // Add elem to the map of parents, since it may have
            // children to be added later
            parents.insert(std::make_pair(self_ID,elem));
	    
	    // Assign the connectivity
	    for (unsigned int n=0; n<elem->n_nodes(); n++)
	      {
		libmesh_assert (cnt < conn.size());
		
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
              libmesh_error();
          }

      } // end if iam != cpu 0
    
    
    libmesh_assert (mesh.n_elem() == n_elem);
  } // Done distributing the elements


  STOP_LOG("broadcast_mesh()","MeshCommunication");
  
#else

  // no MPI but multiple processors? Huh??
  libmesh_error();
  
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
          
#else

  // no MPI but multiple processors? Huh??
  libmesh_error();

#endif  
}



void MeshCommunication::allgather (ParallelMesh& mesh) const
{
  // The mesh should know it's about to be serialized
  libmesh_assert (mesh.is_serial());

  this->allgather_mesh (mesh);
  this->allgather_bcs  (mesh, *(mesh.boundary_info));

  // Inform new elements of their neighbors,
  // while resetting all remote_elem links
  mesh.find_neighbors(true);
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

    Parallel::allgather(n_objects);
    
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
      (packed_elem_header_size*n_elem[libMesh::processor_id()] + local_weight);
						
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

	    pack_element (conn, elem);	    
	  }
      } // ...that was easy.

    libmesh_assert (conn.size() == 
      packed_elem_header_size*n_elem[libMesh::processor_id()] + local_weight);

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
	    
	    const unsigned int
	      first_global_idx = elem_offsets[p],
	      last_global_idx  = first_global_idx + n_elem[p];
	    
	    // Process each element for processor p.
	    // Note this must work in the case when conn_size[p] == 0.
	    while (cnt < (conn_offset[p] + conn_size[p]))
	      {
		// Unpack the element header
#ifdef ENABLE_AMR
	        const unsigned int elem_level = conn[cnt++];
	        const unsigned int p_level    = conn[cnt++];
	        const Elem::RefinementState refinement_flag   =
                  static_cast<Elem::RefinementState>(conn[cnt++]);
	        const Elem::RefinementState p_refinement_flag =
                  static_cast<Elem::RefinementState>(conn[cnt++]);
#else
	        const unsigned int elem_level = 0;
#endif
                const ElemType elem_type      = static_cast<ElemType>(conn[cnt++]);
		const unsigned int elem_PID   = conn[cnt++];
	        const int subdomain_ID        = conn[cnt++];
                const unsigned int self_ID    = conn[cnt++];
		// We require contiguous numbering on each processor
		// for elements.
		libmesh_assert (self_ID >= first_global_idx);
		libmesh_assert (self_ID  < last_global_idx);
#ifdef ENABLE_AMR
                const int parent_ID           = conn[cnt++];
                const int which_child         = conn[cnt++];

		libmesh_assert ((elem_level == 0) || (parent_ID != -1));
  
		// Ignore elements not matching the current level.
		if (elem_level > level) // skip elem->n_nodes() entries in the conn array.
		  cnt += Elem::type_to_n_nodes_map[elem_type];

                else
#endif
		if (elem_level < level ||     // we should already have
		    (elem_level == level &&   // lower level and some
                     mesh.elem(self_ID)))     // ghost elements
		  {
                    // No need to construct a dummy element of this type
		    const Elem* elem = mesh.elem(self_ID);

                    libmesh_assert (elem);
#ifdef ENABLE_AMR
		    libmesh_assert (!elem->parent() ||
                            elem->parent()->id() ==
                            static_cast<unsigned int>(parent_ID));
		    libmesh_assert (elem->p_level()           == p_level);
		    libmesh_assert (elem->refinement_flag()   == refinement_flag);
		    libmesh_assert (elem->p_refinement_flag() == p_refinement_flag);
#endif
		    libmesh_assert (elem->type()              == elem_type);
		    libmesh_assert (elem->subdomain_id()      == subdomain_ID);
		    libmesh_assert (elem->id()                == self_ID);

		    cnt += elem->n_nodes();
		  }
		// Those are the easy cases...
		// now elem_level == level and we don't have it
		else
		  {
		    // Declare the element we will add
		    Elem* elem = NULL;

#ifdef ENABLE_AMR
		    // Maybe find its parent
		    if (level > 0)
		      {
			Elem* my_parent = mesh.elem(parent_ID);

			// If the parent was not previously added, we
			// cannot continue.
			if (my_parent == NULL)
			  {
			    std::cerr << "Parent element with ID " << parent_ID 
				      << " not found." << std::endl; 
			    libmesh_error();
			  }
		
			elem = Elem::build(elem_type,my_parent).release();
		        libmesh_assert (elem);
			my_parent->add_child(elem, which_child);
			libmesh_assert (my_parent->type() == elem->type());
                        libmesh_assert (my_parent->child(which_child) == elem);
		      }
		    else
                      {
                        libmesh_assert (parent_ID == -1);
#endif // ENABLE_AMR
		        elem = Elem::build(elem_type).release();
		        libmesh_assert (elem);
#ifdef ENABLE_AMR
                      }
		      
		    // Assign the IDs
		    elem->set_p_level(p_level);
		    elem->set_refinement_flag(refinement_flag);
		    elem->set_p_refinement_flag(p_refinement_flag);
		    libmesh_assert (elem->level() == static_cast<unsigned int>(level));
#endif
		    elem->subdomain_id() = subdomain_ID;
		    libmesh_assert (elem_PID == p);
		    elem->processor_id() = elem_PID;
		    elem->set_id()       = self_ID;
	    
		    // Assign the connectivity
		    for (unsigned int n=0; n<elem->n_nodes(); n++)
		      {
			libmesh_assert (cnt < conn.size());
			
			elem->set_node(n) = mesh.node_ptr (conn[cnt++]);
		      }

		    // Good to go.  Add to the mesh.
		    mesh.insert_elem(elem);
		    
		  } // end elem_level == level
	      }
	    
	  }   

    // Check the result
    libmesh_assert (global_n_elem == mesh.n_elem());
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
  libmesh_assert (n_bc_ids == global_n_bc_ids);
  
#endif
  

  STOP_LOG  ("allgather_bcs()","MeshCommunication");
}
#endif // HAVE_MPI



void MeshCommunication::delete_remote_elements(ParallelMesh& mesh) const
{
  // The mesh should know it's about to be parallelized
  libmesh_assert (!mesh.is_serial());

  START_LOG("delete_remote_elements()", "MeshCommunication");

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
        local_nodes[elem->node(n)] = true;
      const Elem *a = elem;
      while (a)
        {
          semilocal_elems[a->id()] = true;
          a = a->parent();
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
      const Elem *a = elem;
      while (a)
        {
          semilocal_elems[a->id()] = true;
          a = a->parent();
        }
    }

  // Flag all the elements that share nodes with
  // local and unpartitioned elements, along with their ancestors
  MeshBase::element_iterator nl_elem_it = mesh.not_local_elements_begin(),
                             nl_end     = mesh.not_local_elements_end();
  for (; nl_elem_it != nl_end; ++nl_elem_it)
    {
      Elem *elem = *nl_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        if (local_nodes[elem->node(n)])
          {
            const Elem *a = elem;
            while (a)
              {
                semilocal_elems[a->id()] = true;
                a = a->parent();
              }
            break;
          }
    }

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
          if (!semilocal_elems[elem->id()])
            {
              // Make sure we don't leave any invalid pointers
              elem->make_links_to_me_remote();

              // delete_elem doesn't currently invalidate element
              // iterators... that had better not change
              mesh.delete_elem(elem);
            }
        }
    }

  STOP_LOG("delete_remote_elements()", "MeshCommunication");

  // Now make sure the containers actually shrink - strip
  // any newly-created NULL voids out of the element array
  mesh.renumber_nodes_and_elements();
}




// Pack all this information into one communication to avoid two latency hits
// For each element it is of the form
// [ level p_level r_flag p_flag etype subdomain_id 
//   self_ID parent_ID which_child node_0 node_1 ... node_n]
// We cannot use unsigned int because parent_ID can be negative
void MeshCommunication::pack_element (std::vector<int> &conn, const Elem* elem) const
{
  libmesh_assert (elem != NULL);

#ifdef ENABLE_AMR
  conn.push_back (static_cast<int>(elem->level()));
  conn.push_back (static_cast<int>(elem->p_level()));
  conn.push_back (static_cast<int>(elem->refinement_flag()));
  conn.push_back (static_cast<int>(elem->p_refinement_flag()));
#endif
  conn.push_back (static_cast<int>(elem->type()));
  conn.push_back (static_cast<int>(elem->processor_id()));
  conn.push_back (static_cast<int>(elem->subdomain_id()));
  conn.push_back (elem->id());
		
#ifdef ENABLE_AMR
  // use parent_ID of -1 to indicate a level 0 element
  if (elem->level() == 0)
    {
      conn.push_back(-1);
      conn.push_back(-1);
    }
  else
    {
      conn.push_back(elem->parent()->id());
      conn.push_back(elem->parent()->which_child_am_i(elem));
    }
#endif
  
  for (unsigned int n=0; n<elem->n_nodes(); n++)
    conn.push_back (elem->node(n));		
}
