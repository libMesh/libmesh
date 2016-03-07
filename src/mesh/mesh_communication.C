// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <numeric>

// Local Includes -----------------------------------
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_elem.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/parallel_node.h"
#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/utility.h"
#include "libmesh/remote_elem.h"



//-----------------------------------------------
// anonymous namespace for implementation details
namespace {

using namespace libMesh;

/**
 * Specific weak ordering for Elem *'s to be used in a set.
 * We use the id, but first sort by level.  This guarantees
 * when traversing the set from beginning to end the lower
 * level (parent) elements are encountered first.
 */
struct CompareElemIdsByLevel
{
  bool operator()(const Elem * a,
                  const Elem * b) const
  {
    libmesh_assert (a);
    libmesh_assert (b);
    const unsigned int
      al = a->level(), bl = b->level();
    const dof_id_type
      aid = a->id(),   bid = b->id();

    return (al == bl) ? aid < bid : al < bl;
  }
};

struct SyncNeighbors
{
  typedef std::vector<dof_id_type> datum;

  SyncNeighbors(MeshBase & _mesh) :
    mesh(_mesh) {}

  MeshBase & mesh;

  // Find the neighbor ids for each requested element
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & neighbors)
  {
    neighbors.resize(ids.size());

    for (std::size_t i=0; i != ids.size(); ++i)
      {
        // Look for this element in the mesh
        // We'd better find every element we're asked for
        const Elem * elem = mesh.elem(ids[i]);

        // Return the element's neighbors
        const unsigned int n_neigh = elem->n_neighbors();
        neighbors[i].resize(n_neigh);
        for (unsigned int n = 0; n != n_neigh; ++n)
          {
            const Elem * neigh = elem->neighbor(n);
            if (neigh)
              {
                libmesh_assert_not_equal_to(neigh, remote_elem);
                neighbors[i][n] = neigh->id();
              }
            else
              neighbors[i][n] = DofObject::invalid_id;
          }
      }
  }

  void act_on_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & neighbors)
  {
    for (std::size_t i=0; i != ids.size(); ++i)
      {
        Elem * elem = mesh.elem(ids[i]);

        datum & new_neigh = neighbors[i];

        const unsigned int n_neigh = elem->n_neighbors();
        libmesh_assert_equal_to (n_neigh, new_neigh.size());

        for (unsigned int n = 0; n != n_neigh; ++n)
          {
            const dof_id_type new_neigh_id = new_neigh[n];
            const Elem * old_neigh = elem->neighbor(n);
            if (old_neigh && old_neigh != remote_elem)
              {
                libmesh_assert_equal_to(old_neigh->id(), new_neigh_id);
              }
            else if (new_neigh_id == DofObject::invalid_id)
              {
                libmesh_assert (!old_neigh);
              }
            else
              {
                Elem * neigh = mesh.query_elem(new_neigh_id);
                if (neigh)
                  elem->set_neighbor(n, neigh);
                else
                  elem->set_neighbor(n, const_cast<RemoteElem *>(remote_elem));
              }
          }
      }
  }
};

}


namespace libMesh
{


// ------------------------------------------------------------
// MeshCommunication class members
void MeshCommunication::clear ()
{
  //  _neighboring_processors.clear();
}



#ifndef LIBMESH_HAVE_MPI // avoid spurious gcc warnings
// ------------------------------------------------------------
void MeshCommunication::redistribute (ParallelMesh &) const
{
  // no MPI == one processor, no redistribution
  return;
}
#else
// ------------------------------------------------------------
void MeshCommunication::redistribute (ParallelMesh & mesh) const
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
  libmesh_parallel_only(mesh.comm());
  libmesh_assert (!mesh.is_serial());
  libmesh_assert (MeshTools::n_elem(mesh.unpartitioned_elements_begin(),
                                    mesh.unpartitioned_elements_end()) == 0);

  START_LOG("redistribute()","MeshCommunication");

  // Get a few unique message tags to use in communications; we'll
  // default to some numbers around pi*1000
  Parallel::MessageTag
    nodestag   = mesh.comm().get_unique_tag(3141),
    elemstag   = mesh.comm().get_unique_tag(3142);

  // Figure out how many nodes and elements we have which are assigned to each
  // processor.  send_n_nodes_and_elem_per_proc contains the number of nodes/elements
  // we will be sending to each processor, recv_n_nodes_and_elem_per_proc contains
  // the number of nodes/elements we will be receiving from each processor.
  // Format:
  //  send_n_nodes_and_elem_per_proc[2*pid+0] = number of nodes to send to pid
  //  send_n_nodes_and_elem_per_proc[2*pid+1] = number of elements to send to pid
  std::vector<dof_id_type> send_n_nodes_and_elem_per_proc(2*mesh.n_processors(), 0);

  std::vector<Parallel::Request>
    node_send_requests, element_send_requests;

  for (processor_id_type pid=0; pid<mesh.n_processors(); pid++)
    if (pid != mesh.processor_id()) // don't send to ourselves!!
      {
        // Build up a list of nodes and elements to send to processor pid.
        // We will certainly send all the elements assigned to this processor,
        // but we will also ship off any other elements which touch
        // their nodes.
        std::set<const Node *> connected_nodes;
        {
          MeshBase::const_element_iterator       elem_it  = mesh.pid_elements_begin(pid);
          const MeshBase::const_element_iterator elem_end = mesh.pid_elements_end(pid);

          for (; elem_it!=elem_end; ++elem_it)
            {
              const Elem * elem = *elem_it;

              for (unsigned int n=0; n<elem->n_nodes(); n++)
                connected_nodes.insert (elem->get_node(n));
            }
        }

        std::set<const Elem *, CompareElemIdsByLevel> elements_to_send;
        {
          MeshBase::const_element_iterator       elem_it  = mesh.elements_begin();
          const MeshBase::const_element_iterator elem_end = mesh.elements_end();

          for (; elem_it!=elem_end; ++elem_it)
            {
              const Elem * elem = *elem_it;

              for (unsigned int n=0; n<elem->n_nodes(); n++)
                if (connected_nodes.count(elem->get_node(n)))
                  {
                    elements_to_send.insert (elem);

                    // The remote processor needs all its semilocal
                    // elements' ancestors, and it's possible that
                    // they might not all be connected to the remote
                    // processor's nodes.
                    for (const Elem * parent = elem->parent(); parent;
                         parent = parent->parent())
                      elements_to_send.insert(parent);
                  }
            }
        }

        connected_nodes.clear();
        {
          std::set<const Elem *, CompareElemIdsByLevel>::iterator
            elem_it  = elements_to_send.begin(),
            elem_end = elements_to_send.end();

          for (; elem_it!=elem_end; ++elem_it)
            {
              const Elem * elem = *elem_it;

              for (unsigned int n=0; n<elem->n_nodes(); n++)
                connected_nodes.insert(elem->get_node(n));
            }
        }

        // the number of nodes we will ship to pid
        send_n_nodes_and_elem_per_proc[2*pid+0] =
          cast_int<dof_id_type>(connected_nodes.size());

        // send any nodes off to the destination processor
        if (!connected_nodes.empty())
          {
            node_send_requests.push_back(Parallel::request());

            mesh.comm().send_packed_range (pid,
                                           &mesh,
                                           connected_nodes.begin(),
                                           connected_nodes.end(),
                                           node_send_requests.back(),
                                           nodestag);
          }

        // the number of elements we will send to this processor
        send_n_nodes_and_elem_per_proc[2*pid+1] =
          cast_int<dof_id_type>(elements_to_send.size());

        if (!elements_to_send.empty())
          {
            // send the elements off to the destination processor
            element_send_requests.push_back(Parallel::request());

            mesh.comm().send_packed_range (pid,
                                           &mesh,
                                           elements_to_send.begin(),
                                           elements_to_send.end(),
                                           element_send_requests.back(),
                                           elemstag);
          }
      }

  std::vector<dof_id_type> recv_n_nodes_and_elem_per_proc(send_n_nodes_and_elem_per_proc);

  mesh.comm().alltoall (recv_n_nodes_and_elem_per_proc);

  // In general we will only need to communicate with a subset of the other processors.
  // I can't immediately think of a case where we will send elements but not nodes, but
  // these are only bools and we have the information anyway...
  std::vector<bool>
    send_node_pair(mesh.n_processors(),false), send_elem_pair(mesh.n_processors(),false),
    recv_node_pair(mesh.n_processors(),false), recv_elem_pair(mesh.n_processors(),false);

  unsigned int
    n_send_node_pairs=0,      n_send_elem_pairs=0,
    n_recv_node_pairs=0,      n_recv_elem_pairs=0;

  for (processor_id_type pid=0; pid<mesh.n_processors(); pid++)
    {
      if (send_n_nodes_and_elem_per_proc[2*pid+0]) // we have nodes to send
        {
          send_node_pair[pid] = true;
          n_send_node_pairs++;
        }

      if (send_n_nodes_and_elem_per_proc[2*pid+1]) // we have elements to send
        {
          send_elem_pair[pid] = true;
          n_send_elem_pairs++;
        }

      if (recv_n_nodes_and_elem_per_proc[2*pid+0]) // we have nodes to receive
        {
          recv_node_pair[pid] = true;
          n_recv_node_pairs++;
        }

      if (recv_n_nodes_and_elem_per_proc[2*pid+1]) // we have elements to receive
        {
          recv_elem_pair[pid] = true;
          n_recv_elem_pairs++;
        }
    }
  libmesh_assert_equal_to (n_send_node_pairs, node_send_requests.size());
  libmesh_assert_equal_to (n_send_elem_pairs, element_send_requests.size());

  // Receive nodes.
  // We now know how many processors will be sending us nodes.
  for (unsigned int node_comm_step=0; node_comm_step<n_recv_node_pairs; node_comm_step++)
    // but we don't necessarily want to impose an ordering, so
    // just process whatever message is available next.
    mesh.comm().receive_packed_range (Parallel::any_source,
                                      &mesh,
                                      mesh_inserter_iterator<Node>(mesh),
                                      (Node**)libmesh_nullptr,
                                      nodestag);

  // Receive elements.
  // Similarly we know how many processors are sending us elements,
  // but we don't really care in what order we receive them.
  for (unsigned int elem_comm_step=0; elem_comm_step<n_recv_elem_pairs; elem_comm_step++)
    mesh.comm().receive_packed_range (Parallel::any_source,
                                      &mesh,
                                      mesh_inserter_iterator<Elem>(mesh),
                                      (Elem**)libmesh_nullptr,
                                      elemstag);

  // Wait for all sends to complete
  Parallel::wait (node_send_requests);
  Parallel::wait (element_send_requests);

  // Check on the redistribution consistency
#ifdef DEBUG
  MeshTools::libmesh_assert_equal_n_systems(mesh);

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
void MeshCommunication::gather_neighboring_elements (ParallelMesh & mesh) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (mesh.n_processors() == 1)
    return;

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

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

  // Get a unique message tag to use in communications; we'll default
  // to some numbers around pi*10000
  Parallel::MessageTag
    element_neighbors_tag = mesh.comm().get_unique_tag(31416);

  // Now any element with a NULL neighbor either
  // (i) lives on the physical domain boundary, or
  // (ii) lives on an inter-processor boundary.
  // We will now gather all the elements from adjacent processors
  // which are of the same state, which should address all the type (ii)
  // elements.

  // A list of all the processors which *may* contain neighboring elements.
  // (for development simplicity, just make this the identity map)
  std::vector<processor_id_type> adjacent_processors;
  for (processor_id_type pid=0; pid<mesh.n_processors(); pid++)
    if (pid != mesh.processor_id())
      adjacent_processors.push_back (pid);


  const processor_id_type n_adjacent_processors =
    cast_int<processor_id_type>(adjacent_processors.size());

  //-------------------------------------------------------------------------
  // Let's build a list of all nodes which live on NULL-neighbor sides.
  // For simplicity, we will use a set to build the list, then transfer
  // it to a vector for communication.
  std::vector<dof_id_type> my_interface_node_list;
  std::vector<const Elem *>  my_interface_elements;
  {
    std::set<dof_id_type> my_interface_node_set;

    // since parent nodes are a subset of children nodes, this should be sufficient
    MeshBase::const_element_iterator       it     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator it_end = mesh.active_local_elements_end();

    for (; it != it_end; ++it)
      {
        const Elem * const elem = *it;
        libmesh_assert(elem);

        if (elem->on_boundary()) // denotes *any* side has a NULL neighbor
          {
            my_interface_elements.push_back(elem); // add the element, but only once, even
            // if there are multiple NULL neighbors
            for (unsigned int s=0; s<elem->n_sides(); s++)
              if (elem->neighbor(s) == libmesh_nullptr)
                {
                  UniquePtr<Elem> side(elem->build_side(s));

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

  // we will now send my_interface_node_list to all of the adjacent processors.
  // note that for the time being we will copy the list to a unique buffer for
  // each processor so that we can use a nonblocking send and not access the
  // buffer again until the send completes.  it is my understanding that the
  // MPI 2.1 standard seeks to remove this restriction as unnecessary, so in
  // the future we should change this to send the same buffer to each of the
  // adjacent processors. - BSK 11/17/2008
  std::vector<std::vector<dof_id_type> >
    my_interface_node_xfer_buffers (n_adjacent_processors, my_interface_node_list);
  std::map<processor_id_type, unsigned char> n_comm_steps;

  std::vector<Parallel::Request> send_requests (3*n_adjacent_processors);
  unsigned int current_request = 0;

  for (unsigned int comm_step=0; comm_step<n_adjacent_processors; comm_step++)
    {
      n_comm_steps[adjacent_processors[comm_step]]=1;
      mesh.comm().send (adjacent_processors[comm_step],
                        my_interface_node_xfer_buffers[comm_step],
                        send_requests[current_request++],
                        element_neighbors_tag);
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

  std::vector<dof_id_type> common_interface_node_list;

  // we expect two classess of messages -
  // (1) incoming interface node lists, to which we will reply with our elements
  //     touching nodes in the list, and
  // (2) replies from the requests we sent off previously.
  //  (2.a) - nodes
  //  (2.b) - elements
  // so we expect 3 communications from each adjacent processor.
  // by structuring the communication in this way we hopefully impose no
  // order on the handling of the arriving messages.  in particular, we
  // should be able to handle the case where we receive a request and
  // all replies from processor A before even receiving a request from
  // processor B.

  for (unsigned int comm_step=0; comm_step<3*n_adjacent_processors; comm_step++)
    {
      //------------------------------------------------------------------
      // catch incoming node list
      Parallel::Status
        status(mesh.comm().probe (Parallel::any_source,
                                  element_neighbors_tag));
      const processor_id_type
        source_pid_idx = cast_int<processor_id_type>(status.source()),
        dest_pid_idx   = source_pid_idx;

      //------------------------------------------------------------------
      // first time - incoming request
      if (n_comm_steps[source_pid_idx] == 1)
        {
          n_comm_steps[source_pid_idx]++;

          mesh.comm().receive (source_pid_idx,
                               common_interface_node_list,
                               element_neighbors_tag);
          const std::size_t
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
            libMesh::out << "[" << mesh.processor_id() << "] "
                         << "my_interface_node_list.size()="       << my_interface_node_list.size()
                         << ", [" << source_pid_idx << "] "
                         << "their_interface_node_list.size()="    << their_interface_node_list_size
                         << ", common_interface_node_list.size()=" << common_interface_node_list.size()
                         << std::endl;

          // Now we need to see which of our elements touch the nodes in the list.
          // We will certainly send all the active elements which intersect source_pid_idx,
          // but we will also ship off the other elements in the same family tree
          // as the active ones for data structure consistency.
          //
          // FIXME - shipping full family trees is unnecessary and inefficient.
          //
          // We also ship any nodes connected to these elements.  Note
          // some of these nodes and elements may be replicated from
          // other processors, but that is OK.
          std::set<const Elem *, CompareElemIdsByLevel> elements_to_send;
          std::set<const Node *> connected_nodes;

          // Check for quick return?
          if (common_interface_node_list.empty())
            {
              // let's try to be smart here - if we have no nodes in common,
              // we cannot share elements.  so post the messages expected
              // from us here and go on about our business.
              // note that even though these are nonblocking sends
              // they should complete essentially instantly, because
              // in all cases the send buffers are empty
              mesh.comm().send_packed_range (dest_pid_idx,
                                             &mesh,
                                             connected_nodes.begin(),
                                             connected_nodes.end(),
                                             send_requests[current_request++],
                                             element_neighbors_tag);

              mesh.comm().send_packed_range (dest_pid_idx,
                                             &mesh,
                                             elements_to_send.begin(),
                                             elements_to_send.end(),
                                             send_requests[current_request++],
                                             element_neighbors_tag);

              continue;
            }
          // otherwise, this really *is* an adjacent processor.
          adjacent_processors.push_back(source_pid_idx);

          std::vector<const Elem *> family_tree;

          for (dof_id_type e=0, n_shared_nodes=0; e<my_interface_elements.size(); e++, n_shared_nodes=0)
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
            libmesh_assert (connected_nodes.empty() || !elements_to_send.empty());
            libmesh_assert (!connected_nodes.empty() || elements_to_send.empty());

            // send the nodes off to the destination processor
            mesh.comm().send_packed_range (dest_pid_idx,
                                           &mesh,
                                           connected_nodes.begin(),
                                           connected_nodes.end(),
                                           send_requests[current_request++],
                                           element_neighbors_tag);

            // send the elements off to the destination processor
            mesh.comm().send_packed_range (dest_pid_idx,
                                           &mesh,
                                           elements_to_send.begin(),
                                           elements_to_send.end(),
                                           send_requests[current_request++],
                                           element_neighbors_tag);
          }
        }
      //------------------------------------------------------------------
      // second time - reply of nodes
      else if (n_comm_steps[source_pid_idx] == 2)
        {
          n_comm_steps[source_pid_idx]++;

          mesh.comm().receive_packed_range (source_pid_idx,
                                            &mesh,
                                            mesh_inserter_iterator<Node>(mesh),
                                            (Node**)libmesh_nullptr,
                                            element_neighbors_tag);
        }
      //------------------------------------------------------------------
      // third time - reply of elements
      else if (n_comm_steps[source_pid_idx] == 3)
        {
          n_comm_steps[source_pid_idx]++;

          mesh.comm().receive_packed_range (source_pid_idx,
                                            &mesh,
                                            mesh_inserter_iterator<Elem>(mesh),
                                            (Elem**)libmesh_nullptr,
                                            element_neighbors_tag);
        }
      //------------------------------------------------------------------
      // fourth time - shouldn't happen
      else
        {
          libMesh::err << "ERROR:  unexpected number of replies: "
                       << n_comm_steps[source_pid_idx]
                       << std::endl;
        }
    } // done catching & processing replies associated with tag ~ 100,000pi

  // allow any pending requests to complete
  Parallel::wait (send_requests);

  // We can now find neighbor information for the interfaces between
  // local elements and ghost elements.
  mesh.find_neighbors (/* reset_remote_elements = */ true,
                       /* reset_current_list    = */ false);

  // Ghost elements may not have correct remote_elem neighbor links,
  // and we may not be able to locally infer correct neighbor links to
  // remote elements.  So we synchronize ghost element neighbor links.
  SyncNeighbors nsync(mesh);

  Parallel::sync_dofobject_data_by_id
    (mesh.comm(), mesh.elements_begin(), mesh.elements_end(), nsync);

  STOP_LOG("gather_neighboring_elements()","MeshCommunication");
}
#endif // LIBMESH_HAVE_MPI


#ifndef LIBMESH_HAVE_MPI // avoid spurious gcc warnings
// ------------------------------------------------------------
void MeshCommunication::broadcast (MeshBase &) const
{
  // no MPI == one processor, no need for this method...
  return;
}
#else
// ------------------------------------------------------------
void MeshCommunication::broadcast (MeshBase & mesh) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (mesh.n_processors() == 1)
    return;

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  START_LOG("broadcast()","MeshCommunication");

  // Explicitly clear the mesh on all but processor 0.
  if (mesh.processor_id() != 0)
    mesh.clear();

  // Broadcast nodes
  mesh.comm().broadcast_packed_range(&mesh,
                                     mesh.nodes_begin(),
                                     mesh.nodes_end(),
                                     &mesh,
                                     mesh_inserter_iterator<Node>(mesh));

  // Broadcast elements from coarsest to finest, so that child
  // elements will see their parents already in place.
  unsigned int n_levels = MeshTools::n_levels(mesh);
  mesh.comm().broadcast(n_levels);

  for (unsigned int l=0; l != n_levels; ++l)
    mesh.comm().broadcast_packed_range(&mesh,
                                       mesh.level_elements_begin(l),
                                       mesh.level_elements_end(l),
                                       &mesh,
                                       mesh_inserter_iterator<Elem>(mesh));

  // Make sure mesh_dimension and elem_dimensions are consistent.
  mesh.cache_elem_dims();

  // Broadcast all of the named entity information
  mesh.comm().broadcast(mesh.set_subdomain_name_map());
  mesh.comm().broadcast(mesh.get_boundary_info().set_sideset_name_map());
  mesh.comm().broadcast(mesh.get_boundary_info().set_nodeset_name_map());

  libmesh_assert (mesh.comm().verify(mesh.n_elem()));
  libmesh_assert (mesh.comm().verify(mesh.n_nodes()));

#ifdef DEBUG
  MeshTools::libmesh_assert_valid_procids<Elem>(mesh);
  MeshTools::libmesh_assert_valid_procids<Node>(mesh);
#endif

  STOP_LOG("broadcast()","MeshCommunication");
}
#endif // LIBMESH_HAVE_MPI



#ifndef LIBMESH_HAVE_MPI // avoid spurious gcc warnings
// ------------------------------------------------------------
void MeshCommunication::gather (const processor_id_type, ParallelMesh &) const
{
  // no MPI == one processor, no need for this method...
  return;
}
#else
// ------------------------------------------------------------
void MeshCommunication::gather (const processor_id_type root_id, ParallelMesh & mesh) const
{
  // The mesh should know it's about to be serialized
  libmesh_assert (mesh.is_serial());

  // Check for quick return
  if (mesh.n_processors() == 1)
    return;

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  START_LOG("(all)gather()","MeshCommunication");

  (root_id == DofObject::invalid_processor_id) ?

    mesh.comm().allgather_packed_range (&mesh,
                                        mesh.nodes_begin(),
                                        mesh.nodes_end(),
                                        mesh_inserter_iterator<Node>(mesh)) :

    mesh.comm().gather_packed_range (root_id,
                                     &mesh,
                                     mesh.nodes_begin(),
                                     mesh.nodes_end(),
                                     mesh_inserter_iterator<Node>(mesh));

  // Gather elements from coarsest to finest, so that child
  // elements will see their parents already in place.
  // rank 0 should know n_levels regardless, so this is
  // safe independent of root_id
  unsigned int n_levels = MeshTools::n_levels(mesh);
  mesh.comm().broadcast(n_levels);

  for (unsigned int l=0; l != n_levels; ++l)
    (root_id == DofObject::invalid_processor_id) ?

      mesh.comm().allgather_packed_range (&mesh,
                                          mesh.level_elements_begin(l),
                                          mesh.level_elements_end(l),
                                          mesh_inserter_iterator<Elem>(mesh)) :

      mesh.comm().gather_packed_range (root_id,
                                       &mesh,
                                       mesh.level_elements_begin(l),
                                       mesh.level_elements_end(l),
                                       mesh_inserter_iterator<Elem>(mesh));


  // If we are doing an allgather(), perform sanity check on the result.
  if (root_id == DofObject::invalid_processor_id)
    {
      libmesh_assert (mesh.comm().verify(mesh.n_elem()));
      libmesh_assert (mesh.comm().verify(mesh.n_nodes()));
    }

  // Inform new elements of their neighbors,
  // while resetting all remote_elem links on
  // the appropriate ranks.
  if ((root_id == DofObject::invalid_processor_id) ||
      (mesh.comm().rank() == root_id))
    mesh.find_neighbors(true);

  // All done!
  STOP_LOG("(all)gather()","MeshCommunication");
}
#endif // LIBMESH_HAVE_MPI



// Functor for make_elems_parallel_consistent and
// make_node_ids_parallel_consistent
namespace {

struct SyncIds
{
  typedef dof_id_type datum;
  typedef void (MeshBase::*renumber_obj)(dof_id_type, dof_id_type);

  SyncIds(MeshBase & _mesh, renumber_obj _renumberer) :
    mesh(_mesh),
    renumber(_renumberer) {}

  MeshBase & mesh;
  renumber_obj renumber;
  // renumber_obj & renumber;

  // Find the id of each requested DofObject -
  // Parallel::sync_* already did the work for us
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & ids_out)
  {
    ids_out = ids;
  }

  void act_on_data (const std::vector<dof_id_type> & old_ids,
                    std::vector<datum> & new_ids)
  {
    for (unsigned int i=0; i != old_ids.size(); ++i)
      if (old_ids[i] != new_ids[i])
        (mesh.*renumber)(old_ids[i], new_ids[i]);
  }
};


struct SyncPLevels
{
  typedef unsigned char datum;

  SyncPLevels(MeshBase & _mesh) :
    mesh(_mesh) {}

  MeshBase & mesh;

  // Find the p_level of each requested Elem
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & ids_out)
  {
    ids_out.reserve(ids.size());

    for (unsigned int i=0; i != ids.size(); ++i)
      {
        Elem *elem = mesh.elem(ids[i]);

        ids_out.push_back(elem->p_level());
      }
  }

  void act_on_data (const std::vector<dof_id_type> & old_ids,
                    std::vector<datum> & new_p_levels)
  {
    for (unsigned int i=0; i != old_ids.size(); ++i)
      {
        Elem *elem = mesh.elem(old_ids[i]);

        elem->set_p_level(new_p_levels[i]);
      }
  }
};


#ifdef LIBMESH_ENABLE_UNIQUE_ID
template <typename DofObjSubclass>
struct SyncUniqueIds
{
  typedef unique_id_type datum;
  typedef DofObjSubclass* (MeshBase::*query_obj)(const dof_id_type);

  SyncUniqueIds(MeshBase &_mesh, query_obj _querier) :
    mesh(_mesh),
    query(_querier) {}

  MeshBase &mesh;
  query_obj query;

  // Find the id of each requested DofObject -
  // Parallel::sync_* already did the work for us
  void gather_data (const std::vector<dof_id_type>& ids,
                    std::vector<datum>& ids_out)
  {
    ids_out.reserve(ids.size());

    for (unsigned int i=0; i != ids.size(); ++i)
      {
        DofObjSubclass *d = (mesh.*query)(ids[i]);
        libmesh_assert(d);

        ids_out.push_back(d->unique_id());
      }
  }

  void act_on_data (const std::vector<dof_id_type>& ids,
                    std::vector<datum>& unique_ids)
  {
    for (unsigned int i=0; i != ids.size(); ++i)
      {
        DofObjSubclass *d = (mesh.*query)(ids[i]);
        libmesh_assert(d);

        d->set_unique_id() = unique_ids[i];
      }
  }
};
#endif // LIBMESH_ENABLE_UNIQUE_ID
}



// ------------------------------------------------------------
void MeshCommunication::make_node_ids_parallel_consistent (MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  START_LOG ("make_node_ids_parallel_consistent()", "MeshCommunication");

  SyncIds syncids(mesh, &MeshBase::renumber_node);
  Parallel::sync_node_data_by_element_id
    (mesh, mesh.elements_begin(), mesh.elements_end(), syncids);

  STOP_LOG ("make_node_ids_parallel_consistent()", "MeshCommunication");
}



void MeshCommunication::make_node_unique_ids_parallel_consistent
(MeshBase &mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  START_LOG ("make_node_unique_ids_parallel_consistent()", "MeshCommunication");

  SyncUniqueIds<Node> syncuniqueids(mesh, &MeshBase::query_node_ptr);
  Parallel::sync_dofobject_data_by_id
    (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(),
     syncuniqueids);

  STOP_LOG ("make_node_unique_ids_parallel_consistent()", "MeshCommunication");
#endif
}




// ------------------------------------------------------------
void MeshCommunication::make_elems_parallel_consistent(MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  START_LOG ("make_elems_parallel_consistent()", "MeshCommunication");

  SyncIds syncids(mesh, &MeshBase::renumber_elem);
  Parallel::sync_element_data_by_parent_id
    (mesh, mesh.active_elements_begin(),
     mesh.active_elements_end(), syncids);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  SyncUniqueIds<Elem> syncuniqueids(mesh, &MeshBase::query_elem);
  Parallel::sync_dofobject_data_by_id
    (mesh.comm(), mesh.active_elements_begin(),
     mesh.active_elements_end(), syncuniqueids);
#endif

  STOP_LOG ("make_elems_parallel_consistent()", "MeshCommunication");
}



// ------------------------------------------------------------
void MeshCommunication::make_p_levels_parallel_consistent(MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  START_LOG ("make_p_levels_parallel_consistent()", "MeshCommunication");

  SyncPLevels syncplevels(mesh);
  Parallel::sync_dofobject_data_by_id
    (mesh.comm(), mesh.elements_begin(), mesh.elements_end(),
     syncplevels);

  STOP_LOG ("make_p_levels_parallel_consistent()", "MeshCommunication");
}



// Functors for make_node_proc_ids_parallel_consistent
namespace {

struct SyncProcIds
{
  typedef processor_id_type datum;

  SyncProcIds(MeshBase & _mesh) : mesh(_mesh) {}

  MeshBase & mesh;

  // ------------------------------------------------------------
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & data)
  {
    // Find the processor id of each requested node
    data.resize(ids.size());

    for (std::size_t i=0; i != ids.size(); ++i)
      {
        // Look for this point in the mesh
        // We'd better find every node we're asked for
        Node * node = mesh.node_ptr(ids[i]);

        // Return the node's correct processor id,
        data[i] = node->processor_id();
      }
  }

  // ------------------------------------------------------------
  void act_on_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> proc_ids)
  {
    // Set the ghost node processor ids we've now been informed of
    for (std::size_t i=0; i != ids.size(); ++i)
      {
        Node * node = mesh.node_ptr(ids[i]);
        node->processor_id() = proc_ids[i];
      }
  }
};
}



// ------------------------------------------------------------
void MeshCommunication::make_node_proc_ids_parallel_consistent(MeshBase & mesh)
{
  START_LOG ("make_node_proc_ids_parallel_consistent()", "MeshCommunication");

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

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
  Parallel::sync_node_data_by_element_id
    (mesh, mesh.elements_begin(), mesh.elements_end(), sync);

  STOP_LOG ("make_node_proc_ids_parallel_consistent()", "MeshCommunication");
}



// ------------------------------------------------------------
void MeshCommunication::make_nodes_parallel_consistent (MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

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
  this->make_node_proc_ids_parallel_consistent(mesh);

  // Second, sync up dofobject ids.
  this->make_node_ids_parallel_consistent(mesh);

  // Third, sync up dofobject unique_ids if applicable.
  this->make_node_unique_ids_parallel_consistent(mesh);

  // Finally, correct the processor ids to make DofMap happy
  MeshTools::correct_node_proc_ids(mesh);
}



// ------------------------------------------------------------
void MeshCommunication::delete_remote_elements(ParallelMesh & mesh, const std::set<Elem *> & extra_ghost_elem_ids) const
{
  // The mesh should know it's about to be parallelized
  libmesh_assert (!mesh.is_serial());

  START_LOG("delete_remote_elements()", "MeshCommunication");

#ifdef DEBUG
  // We expect maximum ids to be in sync so we can use them to size
  // vectors
  mesh.comm().verify(mesh.max_node_id());
  mesh.comm().verify(mesh.max_elem_id());
  const dof_id_type par_max_node_id = mesh.parallel_max_node_id();
  const dof_id_type par_max_elem_id = mesh.parallel_max_elem_id();
  libmesh_assert_equal_to (par_max_node_id, mesh.max_node_id());
  libmesh_assert_equal_to (par_max_elem_id, mesh.max_elem_id());
#endif

  // FIXME - should these be "unsorted_set"s?  O(N) is O(N)...
  std::vector<bool> local_nodes(mesh.max_node_id(), false);
  std::vector<bool> semilocal_nodes(mesh.max_node_id(), false);
  std::vector<bool> semilocal_elems(mesh.max_elem_id(), false);

  // We don't want to delete any element that shares a node
  // with or is an ancestor of a local element.

  // Using the local_nodes vector rather than e.g. point_neighbors()
  // gives us correct results even in corner cases, such as where two
  // elements meet only at a corner.  ;-)  Links between boundary and
  // interior elements on mixed dimensional meshes also give us
  // correct ghosting in this way.

  // We also preserve neighbors of local elements - in most cases this
  // is redundant with the node check, but for non-conforming Tet4
  // meshes and non-level-one-conforming 2D+3D meshes it is possible
  // for an element and its coarse neighbor to not share any vertices.

  MeshBase::const_element_iterator l_elem_it = mesh.local_elements_begin(),
    l_end     = mesh.local_elements_end();
  for (; l_elem_it != l_end; ++l_elem_it)
    {
      const Elem * elem = *l_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        {
          dof_id_type nodeid = elem->node(n);
          libmesh_assert_less (nodeid, local_nodes.size());
          local_nodes[nodeid] = true;
        }

      for (unsigned int s=0; s != elem->n_sides(); ++s)
        {
          const Elem * neighbor = elem->neighbor(s);
          if (neighbor)
            semilocal_elems[neighbor->id()] = true;
        }

      while (elem)
        {
          dof_id_type elemid = elem->id();
          libmesh_assert_less (elemid, semilocal_elems.size());
          semilocal_elems[elemid] = true;

          for (unsigned int n=0; n != elem->n_nodes(); ++n)
            semilocal_nodes[elem->node(n)] = true;

          const Elem * parent = elem->parent();
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
      const Elem * elem = *u_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        local_nodes[elem->node(n)] = true;
      while (elem)
        {
          semilocal_elems[elem->id()] = true;

          for (unsigned int n=0; n != elem->n_nodes(); ++n)
            semilocal_nodes[elem->node(n)] = true;

          const Elem * parent = elem->parent();
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
      const Elem * elem = *nl_elem_it;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        if (local_nodes[elem->node(n)])
          {
            while (elem)
              {
                semilocal_elems[elem->id()] = true;

                for (unsigned int nn=0; nn != elem->n_nodes(); ++nn)
                  semilocal_nodes[elem->node(nn)] = true;

                const Elem * parent = elem->parent();
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
    {
      const Elem * elem = *it;
      semilocal_elems[elem->id()] = true;
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        semilocal_nodes[elem->node(n)] = true;
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
          Elem * elem = *lev_elem_it;
          libmesh_assert (elem);
          // Make sure we don't leave any invalid pointers
          if (!semilocal_elems[elem->id()])
            elem->make_links_to_me_remote();

          // Subactive neighbor pointers aren't preservable here
          if (elem->subactive())
            for (unsigned int s=0; s != elem->n_sides(); ++s)
              elem->set_neighbor(s, libmesh_nullptr);

          // delete_elem doesn't currently invalidate element
          // iterators... that had better not change
          if (!semilocal_elems[elem->id()])
            mesh.delete_elem(elem);
        }
    }

  // Delete all the nodes we have no reason to save
  MeshBase::node_iterator node_it  = mesh.nodes_begin(),
    node_end = mesh.nodes_end();
  for (node_it = mesh.nodes_begin(); node_it != node_end; ++node_it)
    {
      Node * node = *node_it;
      libmesh_assert(node);
      if (!semilocal_nodes[node->id()])
        mesh.delete_node(node);
    }

#ifdef DEBUG
  MeshTools::libmesh_assert_valid_refinement_tree(mesh);
#endif

  STOP_LOG("delete_remote_elements()", "MeshCommunication");
}

} // namespace libMesh
