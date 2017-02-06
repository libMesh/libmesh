// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <numeric> // std::accumulate

// LibMesh includes
#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/nemesis_io_helper.h"
#include "libmesh/node.h"
#include "libmesh/parallel.h"
#include "libmesh/utility.h" // is_sorted, deallocate
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_communication.h"

namespace libMesh
{


//-----------------------------------------------
// anonymous namespace for implementation details
namespace {
struct CompareGlobalIdxMappings
{
  // strict weak ordering for a.first -> a.second mapping.  since we can only map to one
  // value only order the first entry
  bool operator()(const std::pair<unsigned int, unsigned int> & a,
                  const std::pair<unsigned int, unsigned int> & b) const
  { return a.first < b.first; }

  // strict weak ordering for a.first -> a.second mapping.  lookups will
  // be in terms of a single integer, which is why we need this method.
  bool operator()(const std::pair<unsigned int, unsigned int> & a,
                  const unsigned int b) const
  { return a.first < b; }
};

// Nemesis & ExodusII use int for all integer values, even the ones which
// should never be negative.  we like to use unsigned as a force of habit,
// this trivial little method saves some typing & also makes sure something
// is not horribly wrong.
template <typename T>
inline unsigned int to_uint ( const T & t )
{
  libmesh_assert_equal_to (t, static_cast<T>(static_cast<unsigned int>(t)));

  return static_cast<unsigned int>(t);
}

// test equality for a.first -> a.second mapping.  since we can only map to one
// value only test the first entry
#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API) && !defined(NDEBUG)
inline bool global_idx_mapping_equality (const std::pair<unsigned int, unsigned int> & a,
                                         const std::pair<unsigned int, unsigned int> & b)
{
  return a.first == b.first;
}
#endif

}



// ------------------------------------------------------------
// Nemesis_IO class members
Nemesis_IO::Nemesis_IO (MeshBase & mesh,
#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
                        bool single_precision
#else
                        bool
#endif
                        ) :
  MeshInput<MeshBase> (mesh, /*is_parallel_format=*/true),
  MeshOutput<MeshBase> (mesh, /*is_parallel_format=*/true),
  ParallelObject (mesh),
#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
  nemhelper(new Nemesis_IO_Helper(*this, false, single_precision)),
#endif
  _timestep(1),
  _verbose (false),
  _append(false)
{
}



// Destructor.  Defined in the C file so we can be sure to get away
// with a forward declaration of Nemesis_IO_Helper in the header file.
Nemesis_IO::~Nemesis_IO ()
{
}



void Nemesis_IO::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;

#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
  // Set the verbose flag in the helper object
  // as well.
  nemhelper->verbose = _verbose;
#endif
}



void Nemesis_IO::append(bool val)
{
  _append = val;
}



#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
void Nemesis_IO::read (const std::string & base_filename)
{
  // On one processor, Nemesis and ExodusII should be equivalent, so
  // let's cowardly defer to that implementation...
  if (this->n_processors() == 1)
    {
      // We can do this in one line but if the verbose flag was set in this
      // object, it will no longer be set... thus no extra print-outs for serial runs.
      // ExodusII_IO(this->mesh()).read (base_filename); // ambiguous when Nemesis_IO is multiply-inherited

      MeshBase & mesh = MeshInput<MeshBase>::mesh();
      ExodusII_IO(mesh).read (base_filename);
      return;
    }

  LOG_SCOPE ("read()","Nemesis_IO");

  // This function must be run on all processors at once
  parallel_object_only();

  if (_verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] ";
      libMesh::out << "Reading Nemesis file on processor: " << this->processor_id() << std::endl;
    }

  // Construct the Nemesis filename based on the number of processors and the
  // current processor ID.
  std::string nemesis_filename = nemhelper->construct_nemesis_filename(base_filename);

  if (_verbose)
    libMesh::out << "Opening file: " << nemesis_filename << std::endl;

  // Open the Exodus file in EX_READ mode
  nemhelper->open(nemesis_filename.c_str(), /*read_only=*/true);

  // Get a reference to the mesh.  We need to be specific
  // since Nemesis_IO is multiply-inherited
  // MeshBase & mesh = this->mesh();
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // We're reading a file on each processor, so our mesh is
  // partitioned into that many parts as it's created
  this->set_n_partitions(this->n_processors());

  // Local information: Read the following information from the standard Exodus header
  //  title[0]
  //  num_dim
  //  num_nodes
  //  num_elem
  //  num_elem_blk
  //  num_node_sets
  //  num_side_sets
  nemhelper->read_header();
  nemhelper->print_header();

  // Get global information: number of nodes, elems, blocks, nodesets and sidesets
  nemhelper->get_init_global();

  // Get "load balance" information.  This includes the number of internal & border
  // nodes and elements as well as the number of communication maps.
  nemhelper->get_loadbal_param();

  // Do some error checking
  if (nemhelper->num_external_nodes)
    libmesh_error_msg("ERROR: there should be no external nodes in an element-based partitioning!");

  libmesh_assert_equal_to (nemhelper->num_nodes,
                           (nemhelper->num_internal_nodes +
                            nemhelper->num_border_nodes));

  libmesh_assert_equal_to (nemhelper->num_elem,
                           (nemhelper->num_internal_elems +
                            nemhelper->num_border_elems));

  libmesh_assert_less_equal (nemhelper->num_nodes, nemhelper->num_nodes_global);
  libmesh_assert_less_equal (nemhelper->num_elem, nemhelper->num_elems_global);

  // Read nodes from the exodus file: this fills the nemhelper->x,y,z arrays.
  nemhelper->read_nodes();

  // Reads the nemhelper->node_num_map array, node_num_map[i] is the global node number for
  // local node number i.
  nemhelper->read_node_num_map();

  // The get_cmap_params() function reads in the:
  //  node_cmap_ids[],
  //  node_cmap_node_cnts[],
  //  elem_cmap_ids[],
  //  elem_cmap_elem_cnts[],
  nemhelper->get_cmap_params();

  // Read the IDs of the interior, boundary, and external nodes.  This function
  // fills the vectors:
  //  node_mapi[],
  //  node_mapb[],
  //  node_mape[]
  nemhelper->get_node_map();

  // Read each node communication map for this processor.  This function
  // fills the vectors of vectors named:
  //  node_cmap_node_ids[][]
  //  node_cmap_proc_ids[][]
  nemhelper->get_node_cmap();

  libmesh_assert_equal_to (to_uint(nemhelper->num_node_cmaps), nemhelper->node_cmap_node_cnts.size());
  libmesh_assert_equal_to (to_uint(nemhelper->num_node_cmaps), nemhelper->node_cmap_node_ids.size());
  libmesh_assert_equal_to (to_uint(nemhelper->num_node_cmaps), nemhelper->node_cmap_proc_ids.size());

#ifndef NDEBUG
  // We expect the communication maps to be symmetric - e.g. if processor i thinks it
  // communicates with processor j, then processor j should also be expecting to
  // communicate with i.  We can assert that here easily enough with an alltoall,
  // but let's only do it when not in optimized mode to limit unnecessary communication.
  {
    std::vector<unsigned char> pid_send_partner (this->n_processors(), 0);

    // strictly speaking, we should expect to communicate with ourself...
    pid_send_partner[this->processor_id()] = 1;

    // mark each processor id we reference with a node cmap
    for (unsigned int cmap=0; cmap<to_uint(nemhelper->num_node_cmaps); cmap++)
      {
        libmesh_assert_less (to_uint(nemhelper->node_cmap_ids[cmap]), this->n_processors());

        pid_send_partner[nemhelper->node_cmap_ids[cmap]] = 1;
      }

    // Copy the send pairing so we can catch the receive paring and
    // test for equality
    const std::vector<unsigned char> pid_recv_partner (pid_send_partner);

    this->comm().alltoall (pid_send_partner);

    libmesh_assert (pid_send_partner == pid_recv_partner);
  }
#endif

  // We now have enough information to infer node ownership.  We start by assuming
  // we own all the nodes on this processor.  We will then interrogate the
  // node cmaps and see if a lower-rank processor is associated with any of
  // our nodes.  If so, then that processor owns the node, not us...
  std::vector<processor_id_type> node_ownership (nemhelper->num_internal_nodes +
                                                 nemhelper->num_border_nodes,
                                                 this->processor_id());

  // a map from processor id to cmap number, to be used later
  std::map<unsigned int, unsigned int> pid_to_cmap_map;

  // For each node_cmap...
  for (unsigned int cmap=0; cmap<to_uint(nemhelper->num_node_cmaps); cmap++)
    {
      // Good time for error checking...
      libmesh_assert_equal_to (to_uint(nemhelper->node_cmap_node_cnts[cmap]),
                               nemhelper->node_cmap_node_ids[cmap].size());

      libmesh_assert_equal_to (to_uint(nemhelper->node_cmap_node_cnts[cmap]),
                               nemhelper->node_cmap_proc_ids[cmap].size());

      // In all the samples I have seen, node_cmap_ids[cmap] is the processor
      // rank of the remote processor...
      const int adjcnt_pid_idx = nemhelper->node_cmap_ids[cmap];

      libmesh_assert_less (adjcnt_pid_idx, this->n_processors());
      libmesh_assert_not_equal_to (adjcnt_pid_idx, this->processor_id());

      // We only expect one cmap per adjacent processor
      libmesh_assert (!pid_to_cmap_map.count(adjcnt_pid_idx));

      pid_to_cmap_map[adjcnt_pid_idx] = cmap;

      // ...and each node in that cmap...
      for (unsigned int idx=0; idx<to_uint(nemhelper->node_cmap_node_cnts[cmap]); idx++)
        {
          //  Are the node_cmap_ids and node_cmap_proc_ids really redundant?
          libmesh_assert_equal_to (adjcnt_pid_idx, nemhelper->node_cmap_proc_ids[cmap][idx]);

          // we are expecting the exodus node numbering to be 1-based...
          const unsigned int local_node_idx = nemhelper->node_cmap_node_ids[cmap][idx]-1;

          libmesh_assert_less (local_node_idx, node_ownership.size());

          // if the adjacent processor is lower rank than the current
          // owner for this node, then it will get the node...
          node_ownership[local_node_idx] =
            std::min(node_ownership[local_node_idx],
                     cast_int<processor_id_type>(adjcnt_pid_idx));
        }
    } // We now should have established proper node ownership.

  // now that ownership is established, we can figure out how many nodes we
  // will be responsible for numbering.
  unsigned int num_nodes_i_must_number = 0;

  for (std::size_t idx=0; idx<node_ownership.size(); idx++)
    if (node_ownership[idx] == this->processor_id())
      num_nodes_i_must_number++;

  // more error checking...
  libmesh_assert_greater_equal (num_nodes_i_must_number, to_uint(nemhelper->num_internal_nodes));
  libmesh_assert (num_nodes_i_must_number <= to_uint(nemhelper->num_internal_nodes +
                                                     nemhelper->num_border_nodes));
  if (_verbose)
    libMesh::out << "[" << this->processor_id() << "] "
                 << "num_nodes_i_must_number="
                 << num_nodes_i_must_number
                 << std::endl;

  // The call to get_loadbal_param() gets 7 pieces of information.  We allgather
  // these now across all processors to determine some global numberings. We should
  // also gather the number of nodes each processor thinks it will number so that
  // we can (i) determine our offset, and (ii) do some error checking.
  std::vector<int> all_loadbal_data ( 8 );
  all_loadbal_data[0] = nemhelper->num_internal_nodes;
  all_loadbal_data[1] = nemhelper->num_border_nodes;
  all_loadbal_data[2] = nemhelper->num_external_nodes;
  all_loadbal_data[3] = nemhelper->num_internal_elems;
  all_loadbal_data[4] = nemhelper->num_border_elems;
  all_loadbal_data[5] = nemhelper->num_node_cmaps;
  all_loadbal_data[6] = nemhelper->num_elem_cmaps;
  all_loadbal_data[7] = num_nodes_i_must_number;

  this->comm().allgather (all_loadbal_data, /* identical_buffer_sizes = */ true);

  // OK, we are now in a position to request new global indices for all the nodes
  // we do not own

  // Let's get a unique message tag to use for send()/receive()
  Parallel::MessageTag nodes_tag = mesh.comm().get_unique_tag(12345);

  std::vector<std::vector<int> >
    needed_node_idxs (nemhelper->num_node_cmaps); // the indices we will ask for

  std::vector<Parallel::Request>
    needed_nodes_requests (nemhelper->num_node_cmaps);

  for (unsigned int cmap=0; cmap<to_uint(nemhelper->num_node_cmaps); cmap++)
    {
      // We know we will need no more indices than there are nodes
      // in this cmap, but that number is an upper bound in general
      // since the neighboring processor associated with the cmap
      //  may not actually own it
      needed_node_idxs[cmap].reserve   (nemhelper->node_cmap_node_cnts[cmap]);

      const unsigned int adjcnt_pid_idx = nemhelper->node_cmap_ids[cmap];

      // ...and each node in that cmap...
      for (unsigned int idx=0; idx<to_uint(nemhelper->node_cmap_node_cnts[cmap]); idx++)
        {
          const unsigned int
            local_node_idx  = nemhelper->node_cmap_node_ids[cmap][idx]-1,
            owning_pid_idx  = node_ownership[local_node_idx];

          // add it to the request list for its owning processor.
          if (owning_pid_idx == adjcnt_pid_idx)
            {
              const unsigned int
                global_node_idx = nemhelper->node_num_map[local_node_idx]-1;
              needed_node_idxs[cmap].push_back(global_node_idx);
            }
        }
      // now post the send for this cmap
      this->comm().send (adjcnt_pid_idx,              // destination
                         needed_node_idxs[cmap],      // send buffer
                         needed_nodes_requests[cmap], // request
                         nodes_tag);
    } // all communication requests for getting updated global indices for border
      // nodes have been initiated

  // Figure out how many nodes each processor thinks it will number and make sure
  // that it adds up to the global number of nodes. Also, set up global node
  // index offsets for each processor.
  std::vector<unsigned int>
    all_num_nodes_i_must_number (this->n_processors());

  for (unsigned int pid=0; pid<this->n_processors(); pid++)
    all_num_nodes_i_must_number[pid] = all_loadbal_data[8*pid + 7];

  // The sum of all the entries in this vector should sum to the number of global nodes
  libmesh_assert (std::accumulate(all_num_nodes_i_must_number.begin(),
                                  all_num_nodes_i_must_number.end(),
                                  0) == nemhelper->num_nodes_global);

  unsigned int my_next_node = 0;
  for (unsigned int pid=0; pid<this->processor_id(); pid++)
    my_next_node += all_num_nodes_i_must_number[pid];

  const unsigned int my_node_offset = my_next_node;

  if (_verbose)
    libMesh::out << "[" << this->processor_id() << "] "
                 << "my_node_offset="
                 << my_node_offset
                 << std::endl;

  // Add internal nodes to the DistributedMesh, using the node ID offset we
  // computed and the current processor's ID.
  for (unsigned int i=0; i<to_uint(nemhelper->num_internal_nodes); ++i)
    {
      const unsigned int local_node_idx  = nemhelper->node_mapi[i]-1;
#ifndef NDEBUG
      const unsigned int owning_pid_idx  = node_ownership[local_node_idx];
#endif

      // an internal node we do not own? huh??
      libmesh_assert_equal_to (owning_pid_idx, this->processor_id());
      libmesh_assert_less (my_next_node, to_uint(nemhelper->num_nodes_global));

      // "Catch" the node pointer after addition, make sure the
      // ID matches the requested value.
      Node * added_node =
        mesh.add_point (Point(nemhelper->x[local_node_idx],
                              nemhelper->y[local_node_idx],
                              nemhelper->z[local_node_idx]),
                        my_next_node,
                        this->processor_id());

      // Make sure the node we added has the ID we thought it would
      if (added_node->id() != my_next_node)
        {
          libMesh::err << "Error, node added with ID " << added_node->id()
                       << ", but we wanted ID " << my_next_node << std::endl;
        }

      // update the local->global index map, when we are done
      // it will be 0-based.
      nemhelper->node_num_map[local_node_idx] = my_next_node++;
    }

  // Now, for the boundary nodes...  We may very well own some of them,
  // but there may be others for which we have requested the new global
  // id.  We expect to be asked for the ids of the ones we own, so
  // we need to create a map from the old global id to the new one
  // we are about to create.
  typedef std::vector<std::pair<unsigned int, unsigned int> > global_idx_mapping_type;
  global_idx_mapping_type old_global_to_new_global_map;
  old_global_to_new_global_map.reserve (num_nodes_i_must_number // total # i will have
                                        - (my_next_node         // amount i have thus far
                                           - my_node_offset));  // this should be exact!
  CompareGlobalIdxMappings global_idx_mapping_comp;

  for (unsigned int i=0; i<to_uint(nemhelper->num_border_nodes); ++i)
    {
      const unsigned int
        local_node_idx  = nemhelper->node_mapb[i]-1,
        owning_pid_idx  = node_ownership[local_node_idx];

      // if we own it...
      if (owning_pid_idx == this->processor_id())
        {
          const unsigned int
            global_node_idx = nemhelper->node_num_map[local_node_idx]-1;

          // we will number it, and create a mapping from its old global index to
          // the new global index, for lookup purposes when neighbors come calling
          old_global_to_new_global_map.push_back(std::make_pair(global_node_idx,
                                                                my_next_node));

          // "Catch" the node pointer after addition, make sure the
          // ID matches the requested value.
          Node * added_node =
            mesh.add_point (Point(nemhelper->x[local_node_idx],
                                  nemhelper->y[local_node_idx],
                                  nemhelper->z[local_node_idx]),
                            my_next_node,
                            this->processor_id());

          // Make sure the node we added has the ID we thought it would
          if (added_node->id() != my_next_node)
            {
              libMesh::err << "Error, node added with ID " << added_node->id()
                           << ", but we wanted ID " << my_next_node << std::endl;
            }

          // update the local->global index map, when we are done
          // it will be 0-based.
          nemhelper->node_num_map[local_node_idx] = my_next_node++;
        }
    }
  // That should cover numbering all the nodes which belong to us...
  libmesh_assert_equal_to (num_nodes_i_must_number, (my_next_node - my_node_offset));

  // Let's sort the mapping so we can efficiently answer requests
  std::sort (old_global_to_new_global_map.begin(),
             old_global_to_new_global_map.end(),
             global_idx_mapping_comp);

  // and it had better be unique...
  libmesh_assert (std::unique (old_global_to_new_global_map.begin(),
                               old_global_to_new_global_map.end(),
                               global_idx_mapping_equality) ==
                  old_global_to_new_global_map.end());

  // We can now catch incoming requests and process them. for efficiency
  // let's do whatever is available next
  std::map<unsigned int, std::vector<int> > requested_node_idxs; // the indices asked of us

  std::vector<Parallel::Request> requested_nodes_requests(nemhelper->num_node_cmaps);

  // We know we will receive the request from a given processor before
  // we receive its reply to our request. However, we may receive
  // a request and a response from one processor before getting
  // a request from another processor.  So what we are doing here
  // is processing whatever message comes next, while recognizing
  // we will receive a request from a processor before receiving
  // its reply
  std::vector<bool> processed_cmap (nemhelper->num_node_cmaps, false);

  for (unsigned int comm_step=0; comm_step<2*to_uint(nemhelper->num_node_cmaps); comm_step++)
    {
      // query the first message which is available
      const Parallel::Status
        status (this->comm().probe (Parallel::any_source,
                                    nodes_tag));
      const unsigned int
        requesting_pid_idx = status.source(),
        source_pid_idx     = status.source();

      // this had better be from a processor we are expecting...
      libmesh_assert (pid_to_cmap_map.count(requesting_pid_idx));

      // the local cmap which corresponds to the source processor
      const unsigned int cmap = pid_to_cmap_map[source_pid_idx];

      if (!processed_cmap[cmap])
        {
          processed_cmap[cmap] = true;

          // we should only get one request per paired processor
          libmesh_assert (!requested_node_idxs.count(requesting_pid_idx));

          // get a reference to the request buffer for this processor to
          // avoid repeated map lookups
          std::vector<int> & xfer_buf (requested_node_idxs[requesting_pid_idx]);

          // actually receive the message.
          this->comm().receive (requesting_pid_idx, xfer_buf, nodes_tag);

          // Fill the request
          for (std::size_t i=0; i<xfer_buf.size(); i++)
            {
              // the requested old global node index, *now 0-based*
              const unsigned int old_global_node_idx = xfer_buf[i];

              // find the new global node index for the requested node -
              // note that requesting_pid_idx thinks we own this node,
              // so we better!
              const global_idx_mapping_type::const_iterator it =
                std::lower_bound (old_global_to_new_global_map.begin(),
                                  old_global_to_new_global_map.end(),
                                  old_global_node_idx,
                                  global_idx_mapping_comp);

              libmesh_assert (it != old_global_to_new_global_map.end());
              libmesh_assert_equal_to (it->first, old_global_node_idx);
              libmesh_assert_greater_equal (it->second, my_node_offset);
              libmesh_assert_less (it->second, my_next_node);

              // overwrite the requested old global node index with the new global index
              xfer_buf[i] = it->second;
            }

          // and send the new global indices back to the processor which asked for them
          this->comm().send (requesting_pid_idx,
                             xfer_buf,
                             requested_nodes_requests[cmap],
                             nodes_tag);
        } // done processing the request

      // this is the second time we have heard from this processor,
      // so it must be its reply to our request
      else
        {
          // a long time ago, we sent off our own requests.  now it is time to catch the
          // replies and get the new global node numbering.  note that for any reply
          // we receive, the corresponding nonblocking send from above *must* have been
          // completed, since the reply is in response to that request!!

          // if we have received a reply, our send *must* have completed
          // (note we never actually need to wait on the request)
          libmesh_assert (needed_nodes_requests[cmap].test());
          libmesh_assert_equal_to (to_uint(nemhelper->node_cmap_ids[cmap]), source_pid_idx);

          // now post the receive for this cmap
          this->comm().receive (source_pid_idx,
                                needed_node_idxs[cmap],
                                nodes_tag);

          libmesh_assert_less_equal (needed_node_idxs[cmap].size(),
                                     nemhelper->node_cmap_node_ids[cmap].size());

          for (std::size_t i=0, j=0; i<nemhelper->node_cmap_node_ids[cmap].size(); i++)
            {
              const unsigned int
                local_node_idx  = nemhelper->node_cmap_node_ids[cmap][i]-1,
                owning_pid_idx  = node_ownership[local_node_idx];

              // if this node is owned by source_pid_idx, its new global id
              // is in the buffer we just received
              if (owning_pid_idx == source_pid_idx)
                {
                  libmesh_assert_less (j, needed_node_idxs[cmap].size());

                  const unsigned int // now 0-based!
                    global_node_idx = needed_node_idxs[cmap][j++];

                  // "Catch" the node pointer after addition, make sure the
                  // ID matches the requested value.
                  Node * added_node =
                    mesh.add_point (Point(nemhelper->x[local_node_idx],
                                          nemhelper->y[local_node_idx],
                                          nemhelper->z[local_node_idx]),
                                    cast_int<dof_id_type>(global_node_idx),
                                    cast_int<processor_id_type>(source_pid_idx));

                  // Make sure the node we added has the ID we thought it would
                  if (added_node->id() != global_node_idx)
                    {
                      libMesh::err << "Error, node added with ID " << added_node->id()
                                   << ", but we wanted ID " << global_node_idx << std::endl;
                    }

                  // update the local->global index map, when we are done
                  // it will be 0-based.
                  nemhelper->node_num_map[local_node_idx] = global_node_idx;

                  // we are not really going to use my_next_node again, but we can
                  // keep incrimenting it to track how many nodes we have added
                  // to the mesh
                  my_next_node++;
                }
            }
        }
    } // end of node index communication loop

  // we had better have added all the nodes we need to!
  libmesh_assert_equal_to ((my_next_node - my_node_offset), to_uint(nemhelper->num_nodes));

  // After all that, we should be done with all node-related arrays *except* the
  // node_num_map, which we have transformed to use our new numbering...
  // So let's clean up the arrays we are done with.
  {
    Utility::deallocate (nemhelper->node_mapi);
    Utility::deallocate (nemhelper->node_mapb);
    Utility::deallocate (nemhelper->node_mape);
    Utility::deallocate (nemhelper->node_cmap_ids);
    Utility::deallocate (nemhelper->node_cmap_node_cnts);
    Utility::deallocate (nemhelper->node_cmap_node_ids);
    Utility::deallocate (nemhelper->node_cmap_proc_ids);
    Utility::deallocate (nemhelper->x);
    Utility::deallocate (nemhelper->y);
    Utility::deallocate (nemhelper->z);
    Utility::deallocate (needed_node_idxs);
    Utility::deallocate (node_ownership);
  }

  Parallel::wait (requested_nodes_requests);
  requested_node_idxs.clear();

  // See what the node count is up to now.
  if (_verbose)
    {
      // Report the number of nodes which have been added locally
      libMesh::out << "[" << this->processor_id() << "] ";
      libMesh::out << "mesh.n_nodes()=" << mesh.n_nodes() << std::endl;

      // Reports the number of nodes that have been added in total.
      libMesh::out << "[" << this->processor_id() << "] ";
      libMesh::out << "mesh.parallel_n_nodes()=" << mesh.parallel_n_nodes() << std::endl;
    }



  // --------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------


  // We can now read in the elements...Exodus stores them in blocks in which all
  // elements have the same geometric type.  This code is adapted directly from exodusII_io.C

  // Assertion: The sum of the border and internal elements on all processors
  // should equal nemhelper->num_elems_global
#ifndef NDEBUG
  {
    int sum_internal_elems=0, sum_border_elems=0;
    for (unsigned int j=3,c=0; c<this->n_processors(); j+=8,++c)
      sum_internal_elems += all_loadbal_data[j];

    for (unsigned int j=4,c=0; c<this->n_processors(); j+=8,++c)
      sum_border_elems += all_loadbal_data[j];

    if (_verbose)
      {
        libMesh::out << "[" << this->processor_id() << "] ";
        libMesh::out << "sum_internal_elems=" << sum_internal_elems << std::endl;

        libMesh::out << "[" << this->processor_id() << "] ";
        libMesh::out << "sum_border_elems=" << sum_border_elems << std::endl;
      }

    libmesh_assert_equal_to (sum_internal_elems+sum_border_elems, nemhelper->num_elems_global);
  }
#endif

  // We need to set the mesh dimension, but the following...
  // mesh.set_mesh_dimension(static_cast<unsigned int>(nemhelper->num_dim));

  // ... is not sufficient since some codes report num_dim==3 for two dimensional
  // meshes living in 3D, even though all the elements are of 2D type.  Therefore,
  // we instead use the dimension of the highest element found for the Mesh dimension,
  // similar to what is done by the Exodus reader, except here it requires a
  // parallel communication.
  elems_of_dimension.resize(4, false); // will use 1-based

  // Compute my_elem_offset, the amount by which to offset the local elem numbering
  // on my processor.
  unsigned int my_next_elem = 0;
  for (unsigned int pid=0; pid<this->processor_id(); ++pid)
    my_next_elem += (all_loadbal_data[8*pid + 3]+  // num_internal_elems, proc pid
                     all_loadbal_data[8*pid + 4]); // num_border_elems, proc pid
  const unsigned int my_elem_offset = my_next_elem;

  if (_verbose)
    libMesh::out << "[" << this->processor_id() << "] "
                 << "my_elem_offset=" << my_elem_offset << std::endl;


  // Fills in the:
  // global_elem_blk_ids[] and
  // global_elem_blk_cnts[] arrays.
  nemhelper->get_eb_info_global();

  //   // Fills in the vectors
  //   // elem_mapi[num_internal_elems]
  //   // elem_mapb[num_border_elems  ]
  //   // These tell which of the (locally-numbered) elements are internal and which are border elements.
  //   // In our test example these arrays are sorted (but non-contiguous), which makes it possible to
  //   // binary search for each element ID... however I don't think we need to distinguish between the
  //   // two types, since either can have nodes the boundary!
  //   nemhelper->get_elem_map();

  // Fills in the vectors of vectors:
  // elem_cmap_elem_ids[][]
  // elem_cmap_side_ids[][]
  // elem_cmap_proc_ids[][]
  // These arrays are of size num_elem_cmaps * elem_cmap_elem_cnts[i], i = 0..num_elem_cmaps
  nemhelper->get_elem_cmap();

  // Get information about the element blocks:
  // (read in the array nemhelper->block_ids[])
  nemhelper->read_block_info();

  // Reads the nemhelper->elem_num_map array, elem_num_map[i] is the global element number for
  // local element number i.
  nemhelper->read_elem_num_map();

  // Instantiate the ElementMaps interface.  This is what translates LibMesh's
  // element numbering scheme to Exodus's.
  ExodusII_IO_Helper::ElementMaps em;

  // Read in the element connectivity for each block by
  // looping over all the blocks.
  for (unsigned int i=0; i<to_uint(nemhelper->num_elem_blk); i++)
    {
      // Read the information for block i:  For nemhelper->block_ids[i], reads
      // elem_type
      // num_elem_this_blk
      // num_nodes_per_elem
      // num_attr
      // connect <-- the nodal connectivity array for each element in the block.
      nemhelper->read_elem_in_block(i);

      // Note that with parallel files it is possible we have no elements in
      // this block!
      if (!nemhelper->num_elem_this_blk) continue;

      // Set subdomain ID based on the block ID.
      subdomain_id_type subdomain_id =
        cast_int<subdomain_id_type>(nemhelper->block_ids[i]);

      // Create a type string (this uses the null-terminated string ctor).
      const std::string type_str ( &(nemhelper->elem_type[0]) );

      // Set any relevant node/edge maps for this element
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(type_str);

      if (_verbose)
        libMesh::out << "Reading a block of " << type_str << " elements." << std::endl;

      // Loop over all the elements in this block
      for (unsigned int j=0; j<to_uint(nemhelper->num_elem_this_blk); j++)
        {
          Elem * elem = Elem::build (conv.get_canonical_type()).release();
          libmesh_assert (elem);

          // Assign subdomain and processor ID to the newly-created Elem.
          // Assigning the processor ID beforehand ensures that the Elem is
          // not added as an "unpartitioned" element.  Note that the element
          // numbering in Exodus is also 1-based.
          elem->subdomain_id() = subdomain_id;
          elem->processor_id() = this->processor_id();
          elem->set_id()       = my_next_elem++;
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          elem->set_unique_id() = elem->id();
#endif

          // Mark that we have seen an element of the current element's
          // dimension.
          elems_of_dimension[elem->dim()] = true;

          // Add the created Elem to the Mesh, catch the Elem
          // pointer that the Mesh throws back.
          elem = mesh.add_elem (elem);

          // We are expecting the element "thrown back" by libmesh to have the ID we specified for it...
          // Check to see that really is the case.  Note that my_next_elem was post-incremented, so
          // subtract 1 when performing the check.
          if (elem->id() != my_next_elem-1)
            libmesh_error_msg("Unexpected ID "  \
                              << elem->id()                             \
                              << " set by parallel mesh. (expecting "   \
                              << my_next_elem-1                         \
                              << ").");

          // Set all the nodes for this element
          if (_verbose)
            libMesh::out << "[" << this->processor_id() << "] "
                         << "Setting nodes for Elem " << elem->id() << std::endl;

          for (unsigned int k=0; k<to_uint(nemhelper->num_nodes_per_elem); k++)
            {
              const unsigned int
                gi              = (j*nemhelper->num_nodes_per_elem +       // index into connectivity array
                                   conv.get_node_map(k)),
                local_node_idx  = nemhelper->connect[gi]-1,                // local node index
                global_node_idx = nemhelper->node_num_map[local_node_idx]; // new global node index

              // Set node number
              elem->set_node(k) = mesh.node_ptr(global_node_idx);
            }
        } // for (unsigned int j=0; j<nemhelper->num_elem_this_blk; j++)
    } // end for (unsigned int i=0; i<nemhelper->num_elem_blk; i++)

  libmesh_assert_equal_to ((my_next_elem - my_elem_offset), to_uint(nemhelper->num_elem));

  if (_verbose)
    {
      // Print local elems_of_dimension information
      for (std::size_t i=1; i<elems_of_dimension.size(); ++i)
        libMesh::out << "[" << this->processor_id() << "] "
                     << "elems_of_dimension[" << i << "]=" << elems_of_dimension[i] << std::endl;
    }

  // Get the max dimension seen on the current processor
  unsigned char max_dim_seen = 0;
  for (std::size_t i=1; i<elems_of_dimension.size(); ++i)
    if (elems_of_dimension[i])
      max_dim_seen = static_cast<unsigned char>(i);

  // Do a global max to determine the max dimension seen by all processors.
  // It should match -- I don't think we even support calculations on meshes
  // with elements of different dimension...
  this->comm().max(max_dim_seen);

  if (_verbose)
    {
      // Print the max element dimension from all processors
      libMesh::out << "[" << this->processor_id() << "] "
                   << "max_dim_seen=" << +max_dim_seen << std::endl;
    }

  // Set the mesh dimension to the largest encountered for an element
  mesh.set_mesh_dimension(max_dim_seen);

#if LIBMESH_DIM < 3
  if (mesh.mesh_dimension() > LIBMESH_DIM)
    libmesh_error_msg("Cannot open dimension "   \
                      << mesh.mesh_dimension()                          \
                      << " mesh file when configured without "          \
                      << mesh.mesh_dimension()                          \
                      << "D support." );
#endif


  // Global sideset information, they are distributed as well, not sure if they will require communication...?
  nemhelper->get_ss_param_global();

  if (_verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] "
                   << "Read global sideset parameter information." << std::endl;

      // These global values should be the same on all processors...
      libMesh::out << "[" << this->processor_id() << "] "
                   << "Number of global sideset IDs: " << nemhelper->global_sideset_ids.size() << std::endl;
    }

  // Read *local* sideset info the same way it is done in
  // exodusII_io_helper.  May be called any time after
  // nem_helper->read_header(); This sets num_side_sets and resizes
  // elem_list, side_list, and id_list to num_elem_all_sidesets.  Note
  // that there appears to be the same number of sidesets in each file
  // but they all have different numbers of entries (some are empty).
  // Note that the sum of "nemhelper->num_elem_all_sidesets" over all
  // processors should equal the sum of the entries in the "num_global_side_counts" array
  // filled up by nemhelper->get_ss_param_global()
  nemhelper->read_sideset_info();

  if (_verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] "
                   << "nemhelper->num_side_sets = " << nemhelper->num_side_sets << std::endl;

      libMesh::out << "[" << this->processor_id() << "] "
                   << "nemhelper->num_elem_all_sidesets = " << nemhelper->num_elem_all_sidesets << std::endl;

      if (nemhelper->num_side_sets > 0)
        {
          libMesh::out << "Sideset names are: ";
          std::map<int, std::string>::iterator
            it = nemhelper->id_to_ss_names.begin(),
            end = nemhelper->id_to_ss_names.end();
          for (; it != end; ++it)
            libMesh::out << "(" << it->first << "," << it->second << ") ";
          libMesh::out << std::endl;
        }
    }

#ifdef DEBUG
  {
    // In DEBUG mode, check that the global number of sidesets reported
    // in each nemesis file matches the sum of all local sideset counts
    // from each processor.  This requires a small communication, so only
    // do it in DEBUG mode.
    int sum_num_global_side_counts = std::accumulate(nemhelper->num_global_side_counts.begin(),
                                                     nemhelper->num_global_side_counts.end(),
                                                     0);

    // MPI sum up the local files contributions
    int sum_num_elem_all_sidesets = nemhelper->num_elem_all_sidesets;
    this->comm().sum(sum_num_elem_all_sidesets);

    if (sum_num_global_side_counts != sum_num_elem_all_sidesets)
      libmesh_error_msg("Error! global side count reported by Nemesis does not " \
                        << "match the side count reported by the individual files!");
  }
#endif

  // Note that exodus stores sidesets in separate vectors but we want to pack
  // them all into a single vector.  So when we call read_sideset(), we pass an offset
  // into the single vector of all IDs
  for (int offset=0, i=0; i<nemhelper->num_side_sets; i++)
    {
      offset += (i > 0 ? nemhelper->num_sides_per_set[i-1] : 0); // Compute new offset
      nemhelper->read_sideset (i, offset);
    }

  // Now that we have the lists of elements, sides, and IDs, we are ready to set them
  // in the BoundaryInfo object of our Mesh object.  This is slightly different in parallel...
  // For example, I think the IDs in each of the split Exodus files are numbered locally,
  // and we have to know the appropriate ID for this processor to be able to set the
  // entry in BoundaryInfo.  This offset should be given by my_elem_offset determined in
  // this function...

  // Debugging:
  // Print entries of elem_list
  // libMesh::out << "[" << this->processor_id() << "] "
  //        << "elem_list = ";
  // for (std::size_t e=0; e<nemhelper->elem_list.size(); e++)
  //   {
  //     libMesh::out << nemhelper->elem_list[e] << ", ";
  //   }
  // libMesh::out << std::endl;

  // Print entries of side_list
  // libMesh::out << "[" << this->processor_id() << "] "
  //        << "side_list = ";
  // for (std::size_t e=0; e<nemhelper->side_list.size(); e++)
  //   {
  //     libMesh::out << nemhelper->side_list[e] << ", ";
  //   }
  // libMesh::out << std::endl;


  // Loop over the entries of the elem_list, get their pointers from the
  // Mesh data structure, and assign the appropriate side to the BoundaryInfo object.
  for (std::size_t e=0; e<nemhelper->elem_list.size(); e++)
    {
      // Calling mesh.elem_ptr() is an error if no element with that
      // id exists on this processor...
      //
      // Perhaps we should iterate over elements and look them up in
      // the elem list instead?  Note that the IDs in elem_list are
      // not necessarily in order, so if we did instead loop over the
      // mesh, we would have to search the (unsorted) elem_list vector
      // for each entry!  We'll settle for doing some error checking instead.
      Elem * elem = mesh.elem_ptr
        (my_elem_offset + 
         (nemhelper->elem_list[e]-1)/*Exodus numbering is 1-based!*/);

      // The side numberings in libmesh and exodus are not 1:1, so we need to map
      // whatever side number is stored in Exodus into a libmesh side number using
      // a conv object...
      const ExodusII_IO_Helper::Conversion conv =
        em.assign_conversion(elem->type());

      // Finally, we are ready to add the element and its side to the BoundaryInfo object.
      // Call the version of add_side which takes a pointer, since we have already gone to
      // the trouble of getting said pointer...
      mesh.get_boundary_info().add_side(elem,
                                        cast_int<unsigned short>(conv.get_side_map(nemhelper->side_list[e]-1)), // Exodus numbering is 1-based
                                        cast_int<boundary_id_type>(nemhelper->id_list[e]));
    }

  // Debugging: make sure there are as many boundary conditions in the
  // boundary ID object as expected.  Note that, at this point, the
  // mesh still thinks it's serial, so n_boundary_conds() returns the
  // local number of boundary conditions (and is therefore cheap)
  // which should match nemhelper->elem_list.size().
  {
    std::size_t nbcs = mesh.get_boundary_info().n_boundary_conds();
    if (nbcs != nemhelper->elem_list.size())
      libmesh_error_msg("[" << this->processor_id() << "] "   \
                        << "BoundaryInfo contains "                     \
                        << nbcs                                         \
                        << " boundary conditions, while the Exodus file had " \
                        << nemhelper->elem_list.size());
  }

  // Read global nodeset parameters?  We might be able to use this to verify
  // something about the local files, but I haven't figured out what yet...
  nemhelper->get_ns_param_global();

  // Read local nodeset info
  nemhelper->read_nodeset_info();

  if (_verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] ";
      libMesh::out << "nemhelper->num_node_sets=" << nemhelper->num_node_sets << std::endl;
      if (nemhelper->num_node_sets > 0)
        {
          libMesh::out << "Nodeset names are: ";
          std::map<int, std::string>::iterator
            it = nemhelper->id_to_ns_names.begin(),
            end = nemhelper->id_to_ns_names.end();
          for (; it != end; ++it)
            libMesh::out << "(" << it->first << "," << it->second << ") ";
          libMesh::out << std::endl;
        }
    }

  //  // Debugging, what is currently in nemhelper->node_num_map anyway?
  //  libMesh::out << "[" << this->processor_id() << "] "
  //       << "nemhelper->node_num_map = ";
  //
  //  for (std::size_t i=0; i<nemhelper->node_num_map.size(); ++i)
  //    libMesh::out << nemhelper->node_num_map[i] << ", ";
  //  libMesh::out << std::endl;

  // For each nodeset,
  for (int nodeset=0; nodeset<nemhelper->num_node_sets; nodeset++)
    {
      // Get the user-defined ID associcated with the nodeset
      int nodeset_id = nemhelper->nodeset_ids[nodeset];

      if (_verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] ";
          libMesh::out << "nemhelper->nodeset_ids[" << nodeset << "]=" << nodeset_id << std::endl;
        }

      // Read the nodeset from file, store them in a vector
      nemhelper->read_nodeset(nodeset);

      // Add nodes from the node_list to the BoundaryInfo object
      for (std::size_t node=0; node<nemhelper->node_list.size(); node++)
        {
          // Don't run past the end of our node map!
          if (to_uint(nemhelper->node_list[node]-1) >= nemhelper->node_num_map.size())
            libmesh_error_msg("Error, index is past the end of node_num_map array!");

          // We should be able to use the node_num_map data structure set up previously to determine
          // the proper global node index.
          unsigned global_node_id = nemhelper->node_num_map[ nemhelper->node_list[node]-1 /*Exodus is 1-based!*/ ];

          if (_verbose)
            {
              libMesh::out << "[" << this->processor_id() << "] "
                           << "nodeset " << nodeset
                           << ", local node number: " << nemhelper->node_list[node]-1
                           << ", global node id: " << global_node_id
                           << std::endl;
            }

          // Add the node to the BoundaryInfo object with the proper nodeset_id
          mesh.get_boundary_info().add_node
            (cast_int<dof_id_type>(global_node_id),
             cast_int<boundary_id_type>(nodeset_id));
        }
    }

  // See what the elem count is up to now.
  if (_verbose)
    {
      // Report the number of elements which have been added locally
      libMesh::out << "[" << this->processor_id() << "] ";
      libMesh::out << "mesh.n_elem()=" << mesh.n_elem() << std::endl;

      // Reports the number of elements that have been added in total.
      libMesh::out << "[" << this->processor_id() << "] ";
      libMesh::out << "mesh.parallel_n_elem()=" << mesh.parallel_n_elem() << std::endl;
    }

  // For DistributedMesh, it seems that _is_serial is true by default.  A hack to
  // make the Mesh think it's parallel might be to call:
  mesh.update_post_partitioning();
  MeshCommunication().make_node_unique_ids_parallel_consistent(mesh);
  mesh.delete_remote_elements();

  // And if that didn't work, then we're actually reading into a
  // ReplicatedMesh, so forget about gathering neighboring elements
  if (mesh.is_serial())
    return;

  // Gather neighboring elements so that the mesh has the proper "ghost" neighbor information.
  MeshCommunication().gather_neighboring_elements(cast_ref<DistributedMesh &>(mesh));
}

#else

void Nemesis_IO::read (const std::string &)
{
  libmesh_error_msg("ERROR, Nemesis API is not defined!");
}

#endif // #if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)





#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)

void Nemesis_IO::write (const std::string & base_filename)
{
  // Get a constant reference to the mesh for writing
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // Create the filename for this processor given the base_filename passed in.
  std::string nemesis_filename = nemhelper->construct_nemesis_filename(base_filename);

  // If the user has set the append flag here, it doesn't really make
  // sense: the intent of this function is to write a Mesh with no
  // data, while "appending" is really intended to add data to an
  // existing file.  If we're verbose, print a message to this effect.
  if (_append && _verbose)
    libMesh::out << "Warning: Appending in Nemesis_IO::write() does not make sense.\n"
                 << "Creating a new file instead!"
                 << std::endl;

  nemhelper->create(nemesis_filename);

  // Initialize data structures and write some global Nemesis-specifc data, such as
  // communication maps, to file.
  nemhelper->initialize(nemesis_filename,mesh);

  // Call the Nemesis-specialized version of write_nodal_coordinates() to write
  // the nodal coordinates.
  nemhelper->write_nodal_coordinates(mesh);

  // Call the Nemesis-specialized version of write_elements() to write
  // the elements.  Note: Must write a zero if a given global block ID has no
  // elements...
  nemhelper->write_elements(mesh);

  // Call our specialized function to write the nodesets
  nemhelper->write_nodesets(mesh);

  // Call our specialized write_sidesets() function to write the sidesets to file
  nemhelper->write_sidesets(mesh);

  // Not sure if this is really necessary, but go ahead and flush the file
  // once we have written all this stuff.
  nemhelper->ex_err = exII::ex_update(nemhelper->ex_id);

  if( (mesh.get_boundary_info().n_edge_conds() > 0) &&
      _verbose )
    {
      libMesh::out << "Warning: Mesh contains edge boundary IDs, but these "
                   << "are not supported by the Nemesis format."
                   << std::endl;
    }
}

#else

void Nemesis_IO::write (const std::string & )
{
  libmesh_error_msg("ERROR, Nemesis API is not defined!");
}

#endif // #if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)


#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)

void Nemesis_IO::write_timestep (const std::string & fname,
                                 const EquationSystems & es,
                                 const int timestep,
                                 const Real time)
{
  _timestep=timestep;
  write_equation_systems(fname,es);

  nemhelper->write_timestep(timestep, time);
}

#else

void Nemesis_IO::write_timestep (const std::string &,
                                 const EquationSystems &,
                                 const int,
                                 const Real)
{
  libmesh_error_msg("ERROR, Nemesis API is not defined!");
}

#endif // #if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)



#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)

void Nemesis_IO::prepare_to_write_nodal_data (const std::string & fname,
                                              const std::vector<std::string> & names)
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  std::string nemesis_filename = nemhelper->construct_nemesis_filename(fname);

  if (!nemhelper->opened_for_writing)
    {
      // If we're appending, open() the file with read_only=false,
      // otherwise create() it and write the contents of the mesh to
      // it.
      if (_append)
        {
          nemhelper->open(nemesis_filename.c_str(), /*read_only=*/false);
          // After opening the file, read the header so that certain
          // fields, such as the number of nodes and the number of
          // elements, are correctly initialized for the subsequent
          // call to write the nodal solution.
          nemhelper->read_header();
        }
      else
        {
          nemhelper->create(nemesis_filename);
          nemhelper->initialize(nemesis_filename,mesh);
          nemhelper->write_nodal_coordinates(mesh);
          nemhelper->write_elements(mesh);
          nemhelper->write_nodesets(mesh);
          nemhelper->write_sidesets(mesh);

          if( (mesh.get_boundary_info().n_edge_conds() > 0) &&
              _verbose )
            {
              libMesh::out << "Warning: Mesh contains edge boundary IDs, but these "
                           << "are not supported by the ExodusII format."
                           << std::endl;
            }

          // If we don't have any nodes written out on this processor,
          // Exodus seems to like us better if we don't try to write out any
          // variable names too...
#ifdef LIBMESH_USE_COMPLEX_NUMBERS

          std::vector<std::string> complex_names = nemhelper->get_complex_names(names);

          nemhelper->initialize_nodal_variables(complex_names);
#else
          nemhelper->initialize_nodal_variables(names);
#endif
        }
    }
}

#else

void Nemesis_IO::prepare_to_write_nodal_data (const std::string & fname,
                                              const std::vector<std::string> & names)
{
  libmesh_error_msg("ERROR, Nemesis API is not defined.");
}

#endif



#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)

void Nemesis_IO::write_nodal_data (const std::string & base_filename,
                                   const NumericVector<Number> & parallel_soln,
                                   const std::vector<std::string> & names)
{
  LOG_SCOPE("write_nodal_data(parallel)", "Nemesis_IO");

  this->prepare_to_write_nodal_data(base_filename, names);

  // Call the new version of write_nodal_solution() that takes a
  // NumericVector directly without localizing.
  nemhelper->write_nodal_solution(parallel_soln, names, _timestep);
}

#else

void Nemesis_IO::write_nodal_data (const std::string &,
                                   const NumericVector<Number> &,
                                   const std::vector<std::string> &)
{
  libmesh_error_msg("ERROR, Nemesis API is not defined.");
}

#endif



#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)

void Nemesis_IO::write_nodal_data (const std::string & base_filename,
                                   const std::vector<Number> & soln,
                                   const std::vector<std::string> & names)
{
  LOG_SCOPE("write_nodal_data(serialized)", "Nemesis_IO");

  this->prepare_to_write_nodal_data(base_filename, names);

  nemhelper->write_nodal_solution(soln, names, _timestep);
}

#else

void Nemesis_IO::write_nodal_data (const std::string &,
                                   const std::vector<Number> &,
                                   const std::vector<std::string> &)
{
  libmesh_error_msg("ERROR, Nemesis API is not defined.");
}

#endif




#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)

void Nemesis_IO::write_global_data (const std::vector<Number> & soln,
                                    const std::vector<std::string> & names)
{
  if (!nemhelper->opened_for_writing)
    libmesh_error_msg("ERROR, Nemesis file must be initialized before outputting global variables.");

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

  std::vector<std::string> complex_names = nemhelper->get_complex_names(names);

  nemhelper->initialize_global_variables(complex_names);

  unsigned int num_values = soln.size();
  unsigned int num_vars = names.size();
  unsigned int num_elems = num_values / num_vars;

  // This will contain the real and imaginary parts and the magnitude
  // of the values in soln
  std::vector<Real> complex_soln(3*num_values);

  for (unsigned i=0; i<num_vars; ++i)
    {
      for (unsigned int j=0; j<num_elems; ++j)
        {
          Number value = soln[i*num_vars + j];
          complex_soln[3*i*num_elems + j] = value.real();
        }
      for (unsigned int j=0; j<num_elems; ++j)
        {
          Number value = soln[i*num_vars + j];
          complex_soln[3*i*num_elems + num_elems +j] = value.imag();
        }
      for (unsigned int j=0; j<num_elems; ++j)
        {
          Number value = soln[i*num_vars + j];
          complex_soln[3*i*num_elems + 2*num_elems + j] = std::abs(value);
        }
    }

  nemhelper->write_global_values(complex_soln, _timestep);

#else

  // Call the Exodus writer implementation
  nemhelper->initialize_global_variables( names );
  nemhelper->write_global_values( soln, _timestep);

#endif

}

#else

void Nemesis_IO::write_global_data (const std::vector<Number> &,
                                    const std::vector<std::string> &)
{
  libmesh_error_msg("ERROR, Nemesis API is not defined.");
}

#endif // #if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)



#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)

void Nemesis_IO::write_information_records (const std::vector<std::string> & records)
{
  if (!nemhelper->opened_for_writing)
    libmesh_error_msg("ERROR, Nemesis file must be initialized before outputting information records.");

  // Call the Exodus writer implementation
  nemhelper->write_information_records( records );
}


#else

void Nemesis_IO::write_information_records ( const std::vector<std::string> & )
{
  libmesh_error_msg("ERROR, Nemesis API is not defined.");
}

#endif // #if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)

} // namespace libMesh
