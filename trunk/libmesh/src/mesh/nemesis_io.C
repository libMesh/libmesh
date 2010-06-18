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


// C++ includes
#include <numeric> // std::accumulate

// LibMesh includes
#include "exodusII_io.h"
#include "nemesis_io.h"
#include "nemesis_io_helper.h"
#include "parallel_mesh.h"
#include "parallel.h"
#include "utility.h" // is_sorted, deallocate


//-----------------------------------------------
// anonymous namespace for implementation details
namespace {
  struct CompareGlobalIdxMappings
  {
    // strict weak ordering for a.first -> a.second mapping.  since we can only map to one
    // value only order the first entry
    bool operator()(const std::pair<unsigned int, unsigned int> &a,
		    const std::pair<unsigned int, unsigned int> &b) const 
    { return a.first < b.first; }
    
    // strict weak ordering for a.first -> a.second mapping.  lookups will
    // be in terms of a single integer, which is why we need this method.
    bool operator()(const std::pair<unsigned int, unsigned int> &a,
		    const unsigned int b) const 
    { return a.first < b; }
  };

  // Nemesis & ExodusII use int for all integer values, even the ones which 
  // should never be negative.  we like to use unsigned as a force of habit,
  // this trivial little method saves some typing & also makes sure something
  // is not horribly wrong.
  template <typename T>
  inline unsigned int to_uint ( const T &t ) 
  {
    libmesh_assert (t == static_cast<T>(static_cast<unsigned int>(t)));

    return static_cast<unsigned int>(t); 
  }

  // test equality for a.first -> a.second mapping.  since we can only map to one
  // value only test the first entry
  inline bool global_idx_mapping_equality (const std::pair<unsigned int, unsigned int> &a,
					   const std::pair<unsigned int, unsigned int> &b)
  { return a.first == b.first; }
}



// ------------------------------------------------------------
// Nemesis_IO class members
Nemesis_IO::Nemesis_IO (ParallelMesh& mesh) :
  MeshInput<ParallelMesh> (mesh, /*is_parallel_format=*/true),
  //MeshOutput<ParallelMesh> (mesh, /*is_parallel_format=*/true)
  _verbose (false)
{
}


Nemesis_IO::~Nemesis_IO ()
{
}



void Nemesis_IO::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;

#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
  // Set the verbose flag in the helper object
  // as well.
  nemhelper.verbose(_verbose);
#endif
}



#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
void Nemesis_IO::read (const std::string& base_filename)
{
  // On one processor, Nemesis and ExodusII should be equivalent, so
  // let's cowardly defer to that implementation...
  if (libMesh::n_processors() == 1)
    {						
      ExodusII_IO(this->mesh()).read (base_filename);
      return;
    }

  START_LOG ("read()","Nemesis_IO");

  // This function must be run on all processors at once
  parallel_only();
  
  if (_verbose)
    {
      libMesh::out << "[" << libMesh::processor_id() << "] ";
      libMesh::out << "Reading Nemesis file on processor: " << libMesh::processor_id() << std::endl;
    }
  
  // Construct a filename string for this processor.
  //
  // FIXME: This assumes you are reading in a mesh on exactly the
  // same number of processors it was written out on!!
  // This should be generalized at some point...
  std::ostringstream file_oss;

  file_oss << base_filename  
	   << '.' << libMesh::n_processors() 
	   << '.' << libMesh::processor_id();

  libMesh::out << "Opening file: " << file_oss.str() << std::endl;

  // Open the Exodus file
  nemhelper.open(file_oss.str().c_str());

  // Get a reference to the ParallelMesh.  
  ParallelMesh& mesh = this->mesh();

  // Local information: Read the following information from the standard Exodus header
  //  title[0]
  //  num_dim
  //  num_nodes
  //  num_elem
  //  num_elem_blk
  //  num_node_sets
  //  num_side_sets
  nemhelper.read_header();
  nemhelper.print_header();
  
  // Be sure number of dimensions is equal to the number of dimensions in the mesh supplied.
  mesh.set_mesh_dimension(static_cast<unsigned int>(nemhelper.num_dim));
  
  // Get global information: number of nodes, elems, blocks, nodesets and sidesets
  nemhelper.get_init_global();
  
  // Get "load balance" information.  This includes the number of internal & border
  // nodes and elements as well as the number of communication maps.
  nemhelper.get_loadbal_param();
  
  // Do some error checking
  if (nemhelper.num_external_nodes)
    {
      libMesh::err << "ERROR: there should be no external nodes in an element-based partitioning!"
		    << std::endl;
      libmesh_error();
    }

  libmesh_assert (nemhelper.num_nodes ==
		  (nemhelper.num_internal_nodes + 
		   nemhelper.num_border_nodes));
  
  libmesh_assert (nemhelper.num_elem ==
		  (nemhelper.num_internal_elems +
		   nemhelper.num_border_elems));
  
  libmesh_assert (nemhelper.num_nodes <= nemhelper.num_nodes_global);
  libmesh_assert (nemhelper.num_elem  <= nemhelper.num_elems_global);
      
  // Read nodes from the exodus file: this fills the nemhelper.x,y,z arrays.
  nemhelper.read_nodes();

  // Reads the nemhelper.node_num_map array, node_num_map[i] is the global node number for
  // local node number i.
  nemhelper.read_node_num_map();
  
  // The get_cmap_params() function reads in the:
  //  node_cmap_ids[],
  //  node_cmap_node_cnts[],
  //  elem_cmap_ids[],
  //  elem_cmap_elem_cnts[],
  nemhelper.get_cmap_params();
  
  // Read the IDs of the interior, boundary, and external nodes.  This function
  // fills the vectors:
  //  node_mapi[],
  //  node_mapb[],
  //  node_mape[]
  nemhelper.get_node_map();
  
  // Read each node communication map for this processor.  This function
  // fills the vectors of vectors named:
  //  node_cmap_node_ids[][]
  //  node_cmap_proc_ids[][]
  nemhelper.get_node_cmap();

  libmesh_assert (to_uint(nemhelper.num_node_cmaps) == nemhelper.node_cmap_node_cnts.size());
  libmesh_assert (to_uint(nemhelper.num_node_cmaps) == nemhelper.node_cmap_node_ids.size());
  libmesh_assert (to_uint(nemhelper.num_node_cmaps) == nemhelper.node_cmap_proc_ids.size());

#ifndef NDEBUG
  // We expect the communication maps to be symmetric - e.g. if processor i thinks it
  // communicates with processor j, then processor j should also be expecting to
  // communicate with i.  We can assert that here easily enough with an alltoall,
  // but let's only do it when not in optimized mode to limit unnecessary communication.
  {
    std::vector<unsigned char> pid_send_partener (libMesh::n_processors(), 0);

    // strictly speaking, we should expect to communicate with ourself...
    pid_send_partener[libMesh::processor_id()] = 1;
    
    // mark each processor id we reference with a node cmap 
    for (unsigned int cmap=0; cmap<to_uint(nemhelper.num_node_cmaps); cmap++)
      {
	libmesh_assert (to_uint(nemhelper.node_cmap_ids[cmap]) < libMesh::n_processors());

	pid_send_partener[nemhelper.node_cmap_ids[cmap]] = 1;
      }

    // Copy the send pairing so we can catch the receive paring and
    // test for equality
    const std::vector<unsigned char> pid_recv_partener (pid_send_partener);
    
    Parallel::alltoall (pid_send_partener);

    libmesh_assert (pid_send_partener == pid_recv_partener);
  }
#endif
  
  // We now have enough information to infer node ownership.  We start by assuming
  // we own all the nodes on this processor.  We will then interrogate the 
  // node cmaps and see if a lower-rank processor is associated with any of
  // our nodes.  If so, then that processor owns the node, not us...
  std::vector<unsigned short int> node_ownership (nemhelper.num_internal_nodes +
						  nemhelper.num_border_nodes,
						  libMesh::processor_id());

  // a map from processor id to cmap number, to be used later
  std::map<unsigned int, unsigned int> pid_to_cmap_map;

  // For each node_cmap...
  for (unsigned int cmap=0; cmap<to_uint(nemhelper.num_node_cmaps); cmap++)
    {
      // Good time for error checking...
      libmesh_assert (to_uint(nemhelper.node_cmap_node_cnts[cmap]) ==
		      nemhelper.node_cmap_node_ids[cmap].size());
      
      libmesh_assert (to_uint(nemhelper.node_cmap_node_cnts[cmap]) ==
		      nemhelper.node_cmap_proc_ids[cmap].size());
      
      // In all the samples I have seen, node_cmap_ids[cmap] is the processor
      // rank of the remote processor...
      const unsigned short int adjcnt_pid_idx = nemhelper.node_cmap_ids[cmap];

      libmesh_assert (adjcnt_pid_idx <  libMesh::n_processors());
      libmesh_assert (adjcnt_pid_idx != libMesh::processor_id());
      
      // We only expect one cmap per adjacent processor
      libmesh_assert (!pid_to_cmap_map.count(adjcnt_pid_idx));
						
      pid_to_cmap_map[adjcnt_pid_idx] = cmap;

      // ...and each node in that cmap...
      for (unsigned int idx=0; idx<to_uint(nemhelper.node_cmap_node_cnts[cmap]); idx++)
	{
	  //  Are the node_cmap_ids and node_cmap_proc_ids really redundant?
	  libmesh_assert (adjcnt_pid_idx == nemhelper.node_cmap_proc_ids[cmap][idx]);
	  
	  // we are expecting the exodus node numbering to be 1-based...
	  const unsigned int local_node_idx = nemhelper.node_cmap_node_ids[cmap][idx]-1;

	  libmesh_assert (local_node_idx < node_ownership.size());
	  
	  // if the adjacent processor is lower rank than the current
	  // owner for this node, then it will get the node...
	  node_ownership[local_node_idx] = 
	    std::min(node_ownership[local_node_idx], adjcnt_pid_idx);
	}
    } // We now should have established proper node ownership.

  // now that ownership is established, we can figure out how many nodes we 
  // will be responsible for numbering.
  unsigned int num_nodes_i_must_number = 0;

  for (unsigned int idx=0; idx<node_ownership.size(); idx++)
    if (node_ownership[idx] == libMesh::processor_id())
      num_nodes_i_must_number++;

  // more error checking...
  libmesh_assert (num_nodes_i_must_number >= to_uint(nemhelper.num_internal_nodes));
  libmesh_assert (num_nodes_i_must_number <= to_uint(nemhelper.num_internal_nodes +
						  nemhelper.num_border_nodes));
  if (_verbose)
    libMesh::out << "[" << libMesh::processor_id() << "] "
	          << "num_nodes_i_must_number="
	          << num_nodes_i_must_number
	          << std::endl;
  
  // The call to get_loadbal_param() gets 7 pieces of information.  We allgather
  // these now across all processors to determine some global numberings. We should
  // also gather the number of nodes each processor thinks it will number so that
  // we can (i) determine our offset, and (ii) do some error checking.
  std::vector<int> all_loadbal_data ( 8 );
  all_loadbal_data[0] = nemhelper.num_internal_nodes;
  all_loadbal_data[1] = nemhelper.num_border_nodes;
  all_loadbal_data[2] = nemhelper.num_external_nodes;
  all_loadbal_data[3] = nemhelper.num_internal_elems;
  all_loadbal_data[4] = nemhelper.num_border_elems;
  all_loadbal_data[5] = nemhelper.num_node_cmaps;
  all_loadbal_data[6] = nemhelper.num_elem_cmaps;
  all_loadbal_data[7] = num_nodes_i_must_number;
  
  Parallel::allgather (all_loadbal_data, /* identical_buffer_sizes = */ true);

  // OK, we are now in a position to request new global indices for all the nodes
  // we do not own

  // Let's get a unique message tag to use for send()/receive()
  Parallel::MessageTag nodes_tag = Parallel::Communicator_World.get_unique_tag(12345);

  std::vector<std::vector<int> > 
    needed_node_idxs (nemhelper.num_node_cmaps); // the indices we will ask for
    
  std::vector<Parallel::Request> 
    needed_nodes_requests (nemhelper.num_node_cmaps);
  
  for (unsigned int cmap=0; cmap<to_uint(nemhelper.num_node_cmaps); cmap++)
    {
      // We know we will need no more indices than there are nodes
      // in this cmap, but that number is an upper bound in general
      // since the neighboring processor associated with the cmap
      //  may not actually own it
      needed_node_idxs[cmap].reserve   (nemhelper.node_cmap_node_cnts[cmap]);
            
      const unsigned int adjcnt_pid_idx = nemhelper.node_cmap_ids[cmap];

      // ...and each node in that cmap...
      for (unsigned int idx=0; idx<to_uint(nemhelper.node_cmap_node_cnts[cmap]); idx++)
	{
	  const unsigned int 
	    local_node_idx  = nemhelper.node_cmap_node_ids[cmap][idx]-1,
	    owning_pid_idx  = node_ownership[local_node_idx];

	  // add it to the request list for its owning processor.
	  if (owning_pid_idx == adjcnt_pid_idx)
	    {
	      const unsigned int
		global_node_idx = nemhelper.node_num_map[local_node_idx]-1;
	      needed_node_idxs[cmap].push_back(global_node_idx);
	    }
	}
      // now post the send for this cmap
      Parallel::send (adjcnt_pid_idx,              // destination
		      needed_node_idxs[cmap],      // send buffer
		      needed_nodes_requests[cmap], // request
		      nodes_tag);      
    } // all communication requests for getting updated global indices for border
      // nodes have been initiated

  // Figure out how many nodes each processor thinks it will number and make sure
  // that it adds up to the global number of nodes. Also, set up global node
  // index offsets for each processor.
  std::vector<unsigned int>
    all_num_nodes_i_must_number (libMesh::n_processors());
   
  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
    all_num_nodes_i_must_number[pid] = all_loadbal_data[8*pid + 7];
  
  // The sum of all the entries in this vector should sum to the number of global nodes
  libmesh_assert (std::accumulate(all_num_nodes_i_must_number.begin(),
				  all_num_nodes_i_must_number.end(),
				  0) == nemhelper.num_nodes_global);
  
  unsigned int my_next_node = 0;
  for (unsigned int pid=0; pid<libMesh::processor_id(); pid++)
    my_next_node += all_num_nodes_i_must_number[pid];
  
  const unsigned int my_node_offset = my_next_node;

  if (_verbose)
    libMesh::out << "[" << libMesh::processor_id() << "] "
	          << "my_node_offset="
	          << my_node_offset
	          << std::endl;
  
  // Add internal nodes to the ParallelMesh, using the node ID offset we
  // computed and the current processor's ID.
  for (unsigned int i=0; i<to_uint(nemhelper.num_internal_nodes); ++i)
    {
      const unsigned int 
	local_node_idx  = nemhelper.node_mapi[i]-1,
	global_node_idx = nemhelper.node_num_map[local_node_idx]-1,
	owning_pid_idx  = node_ownership[local_node_idx];

      // an internal node we do not own? huh??
      libmesh_assert (owning_pid_idx == libMesh::processor_id());
      libmesh_assert (global_node_idx < to_uint(nemhelper.num_nodes_global));

      mesh.add_point (Point(nemhelper.x[local_node_idx],
			    nemhelper.y[local_node_idx],
			    nemhelper.z[local_node_idx]), 
		      my_next_node, 
		      libMesh::processor_id());

      // update the local->global index map, when we are done
      // it will be 0-based.
      nemhelper.node_num_map[local_node_idx] = my_next_node++;
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

  for (unsigned int i=0; i<to_uint(nemhelper.num_border_nodes); ++i)
    {
      const unsigned int 
	local_node_idx  = nemhelper.node_mapb[i]-1,
	owning_pid_idx  = node_ownership[local_node_idx];	
      
      // if we own it...
      if (owning_pid_idx == libMesh::processor_id())
	{
	  const unsigned int
	    global_node_idx = nemhelper.node_num_map[local_node_idx]-1;

	  // we will number it, and create a mapping from its old global index to
	  // the new global index, for lookup purposes when neighbors come calling
	  old_global_to_new_global_map.push_back(std::make_pair(global_node_idx, 
								my_next_node));
	  mesh.add_point (Point(nemhelper.x[local_node_idx],
				nemhelper.y[local_node_idx],
				nemhelper.z[local_node_idx]), 
			  my_next_node, 
			  libMesh::processor_id());

	  // update the local->global index map, when we are done
	  // it will be 0-based.
	  nemhelper.node_num_map[local_node_idx] = my_next_node++;
	}
    }
  // That should cover numbering all the nodes which belong to us...
  libmesh_assert (num_nodes_i_must_number == (my_next_node - my_node_offset));

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

  std::vector<Parallel::Request> requested_nodes_requests(nemhelper.num_node_cmaps);
  
  // We know we will receive the request from a given processor before
  // we receive its reply to our request. However, we may receive
  // a request and a response from one processor before getting
  // a request from another processor.  So what we are doing here
  // is processing whatever message comes next, while recognizing 
  // we will receive a request from a processor before receiving
  // its reply
  std::vector<bool> processed_cmap (nemhelper.num_node_cmaps, false);

  for (unsigned int comm_step=0; comm_step<2*to_uint(nemhelper.num_node_cmaps); comm_step++)
    {
      // query the first message which is available
      const Parallel::Status 
	status (Parallel::probe (Parallel::any_source, 
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
	  std::vector<int> &xfer_buf (requested_node_idxs[requesting_pid_idx]);
	    
	  // actually receive the message.  
	  Parallel::receive (requesting_pid_idx, xfer_buf, nodes_tag);
	  
	  // Fill the request
	  for (unsigned int i=0; i<xfer_buf.size(); i++)
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
	      libmesh_assert (it->first == old_global_node_idx);
	      libmesh_assert (it->second >= my_node_offset);
	      libmesh_assert (it->second <  my_next_node);
	      
	      // overwrite the requested old global node index with the new global index
	      xfer_buf[i] = it->second;
	    }
	    
	  // and send the new global indices back to the processor which asked for them
	  Parallel::send (requesting_pid_idx,
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
	  libmesh_assert (to_uint(nemhelper.node_cmap_ids[cmap]) == source_pid_idx);
	    
	  // now post the receive for this cmap
	  Parallel::receive (source_pid_idx,
			     needed_node_idxs[cmap],
			     nodes_tag);
	  
	  libmesh_assert (needed_node_idxs[cmap].size() <= 
			  nemhelper.node_cmap_node_ids[cmap].size());
	  
	  for (unsigned int i=0,j=0; i<nemhelper.node_cmap_node_ids[cmap].size(); i++)
	    {
	      const unsigned int 
		local_node_idx  = nemhelper.node_cmap_node_ids[cmap][i]-1,
		owning_pid_idx  = node_ownership[local_node_idx];	
	      
	      // if this node is owned by source_pid_idx, its new global id 
	      // is in the buffer we just received
	      if (owning_pid_idx == source_pid_idx)
		{
		  libmesh_assert (j < needed_node_idxs[cmap].size());
		    
		  const unsigned int // now 0-based!
		    global_node_idx = needed_node_idxs[cmap][j++];
		  
		  mesh.add_point (Point(nemhelper.x[local_node_idx],
					nemhelper.y[local_node_idx],
					nemhelper.z[local_node_idx]),
				  global_node_idx,
				  source_pid_idx);
		  
		  // update the local->global index map, when we are done
		  // it will be 0-based.
		  nemhelper.node_num_map[local_node_idx] = global_node_idx;
		  
		  // we are not really going to use my_next_node again, but we can
		  // keep incrimenting it to track how many nodes we have added 
		  // to the mesh
		  my_next_node++;
		}
	    }
	}
    } // end of node index communication loop
  
  // we had better have added all the nodes we need to!
  libmesh_assert ((my_next_node - my_node_offset) == to_uint(nemhelper.num_nodes));

  // After all that, we should be done with all node-related arrays *except* the
  // node_num_map, which we have transformed to use our new numbering...
  // So let's clean up the arrays we are done with.
  {
    Utility::deallocate (nemhelper.node_mapi);
    Utility::deallocate (nemhelper.node_mapb);
    Utility::deallocate (nemhelper.node_mape);
    Utility::deallocate (nemhelper.node_cmap_ids);
    Utility::deallocate (nemhelper.node_cmap_node_cnts);
    Utility::deallocate (nemhelper.node_cmap_node_ids);
    Utility::deallocate (nemhelper.node_cmap_proc_ids);
    Utility::deallocate (nemhelper.x);
    Utility::deallocate (nemhelper.y);
    Utility::deallocate (nemhelper.z);
    Utility::deallocate (needed_node_idxs);
    Utility::deallocate (node_ownership);
  }
  
  Parallel::wait (requested_nodes_requests);
  requested_node_idxs.clear();

  // See what the node count is up to now.
  if (_verbose)
    {
      // Report the number of nodes which have been added locally
      libMesh::out << "[" << libMesh::processor_id() << "] ";
      libMesh::out << "mesh.n_nodes()=" << mesh.n_nodes() << std::endl;
      
      // Reports the number of nodes that have been added in total.
      libMesh::out << "[" << libMesh::processor_id() << "] ";
      libMesh::out << "mesh.parallel_n_nodes()=" << mesh.parallel_n_nodes() << std::endl;
    }
  


  // --------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------

      
  // We can now read in the elements...Exodus stores them in blocks in which all
  // elements have the same geometric type.  This code is adapted directly from exodusII_io.C

  // Assertion: The sum of the border and internal elements on all processors
  // should equal nemhelper.num_elems_global
#ifndef NDEBUG
  {
    int sum_internal_elems=0, sum_border_elems=0;
    for (unsigned int j=3,c=0; c<libMesh::n_processors(); j+=8,++c)
      sum_internal_elems += all_loadbal_data[j];
    
    for (unsigned int j=4,c=0; c<libMesh::n_processors(); j+=8,++c)
      sum_border_elems += all_loadbal_data[j];
    
    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() << "] ";
	libMesh::out << "sum_internal_elems=" << sum_internal_elems << std::endl;
	
	libMesh::out << "[" << libMesh::processor_id() << "] ";
	libMesh::out << "sum_border_elems=" << sum_border_elems << std::endl;
      }
    
    libmesh_assert(sum_internal_elems+sum_border_elems == nemhelper.num_elems_global);
  }
#endif

  // Compute my_elem_offset, the amount by which to offset the local elem numbering
  // on my processor.
  unsigned int my_next_elem = 0;
  for (unsigned int pid=0; pid<libMesh::processor_id(); ++pid)
    my_next_elem += (all_loadbal_data[8*pid + 3]+  // num_internal_elems, proc pid
		     all_loadbal_data[8*pid + 4]); // num_border_elems, proc pid
  const unsigned int my_elem_offset = my_next_elem;

  if (_verbose)
    libMesh::out << "[" << libMesh::processor_id() << "] "
	      << "my_elem_offset=" << my_elem_offset << std::endl;


  // Fills in the: 
  // global_elem_blk_ids[] and
  // global_elem_blk_cnts[] arrays.
  nemhelper.get_eb_info_global();

//   // Fills in the vectors
//   // elem_mapi[num_internal_elems]
//   // elem_mapb[num_border_elems  ]
//   // These tell which of the (locally-numbered) elements are internal and which are border elements.
//   // In our test example these arrays are sorted (but non-contiguous), which makes it possible to
//   // binary search for each element ID... however I don't think we need to distinguish between the
//   // two types, since either can have nodes the boundary!
//   nemhelper.get_elem_map();
      
  // Fills in the vectors of vectors:
  // elem_cmap_elem_ids[][]
  // elem_cmap_side_ids[][]
  // elem_cmap_proc_ids[][]
  // These arrays are of size num_elem_cmaps * elem_cmap_elem_cnts[i], i = 0..num_elem_cmaps
  nemhelper.get_elem_cmap();
      
  // Get information about the element blocks:
  // (read in the array nemhelper.block_ids[])
  nemhelper.read_block_info();

  // Reads the nemhelper.elem_num_map array, elem_num_map[i] is the global element number for
  // local element number i.
  nemhelper.read_elem_num_map();
      
  // Instantiate the ElementMaps interface.  This is what translates LibMesh's
  // element numbering scheme to Exodus's.
  ExodusII_IO_Helper::ElementMaps em;
      
  // Read in the element connectivity for each block by
  // looping over all the blocks.
  for (unsigned int i=0; i<to_uint(nemhelper.num_elem_blk); i++)
    {
      // Read the information for block i:  For nemhelper.block_ids[i], reads
      // elem_type
      // num_elem_this_blk
      // num_nodes_per_elem
      // num_attr
      // connect <-- the nodal connectivity array for each element in the block.
      nemhelper.read_elem_in_block(i);

      // Note that with parallel files it is possible we have no elements in
      // this block!
      if (!nemhelper.num_elem_this_blk) continue;
      
      // Set subdomain ID based on the block ID.
      int subdomain_id = nemhelper.block_ids[i];

      // Create a type string (this uses the null-terminated string ctor).
      const std::string type_str ( &(nemhelper.elem_type[0]) ); 

      // Set any relevant node/edge maps for this element
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(type_str); 

      if (_verbose)
	libMesh::out << "Reading a block of " << type_str << " elements." << std::endl;
      
      // Loop over all the elements in this block
      for (unsigned int j=0; j<to_uint(nemhelper.num_elem_this_blk); j++)
	{
	  Elem* elem = Elem::build (conv.get_canonical_type()).release();
	  libmesh_assert (elem);

	  // Assign subdomain and processor ID to the newly-created Elem.
	  // Assigning the processor ID beforehand ensures that the Elem is
	  // not added as an "unpartitioned" element.  Note that the element
	  // numbering in Exodus is also 1-based.
	  elem->subdomain_id() = subdomain_id;
	  elem->processor_id() = libMesh::processor_id();
	  elem->set_id()       = my_next_elem++;
		
	  // Add the created Elem to the Mesh, catch the Elem
	  // pointer that the Mesh throws back.
	  elem = mesh.add_elem (elem); 

	  // Set all the nodes for this element
	  if (_verbose)	    
	    libMesh::out << "[" << libMesh::processor_id() << "] "
		          << "Setting nodes for Elem " << elem->id() << std::endl;
	      
	  for (unsigned int k=0; k<to_uint(nemhelper.num_nodes_per_elem); k++)
	    {
	      const unsigned int
		gi              = (j*nemhelper.num_nodes_per_elem +       // index into connectivity array
				   conv.get_node_map(k)), 
		local_node_idx  = nemhelper.connect[gi]-1,                // local node index 
		global_node_idx = nemhelper.node_num_map[local_node_idx]; // new global node index
		  
	      // Set node number
	      elem->set_node(k) = mesh.node_ptr(global_node_idx);
	    }
	} // for (unsigned int j=0; j<nemhelper.num_elem_this_blk; j++)     
    } // end for (unsigned int i=0; i<nemhelper.num_elem_blk; i++)

  libmesh_assert ((my_next_elem - my_elem_offset) == to_uint(nemhelper.num_elem));

  // See what the elem count is up to now.
  if (_verbose)
    {
      // Report the number of elements which have been added locally
      libMesh::out << "[" << libMesh::processor_id() << "] ";
      libMesh::out << "mesh.n_elem()=" << mesh.n_elem() << std::endl;

      // Reports the number of elements that have been added in total.
      libMesh::out << "[" << libMesh::processor_id() << "] ";
      libMesh::out << "mesh.parallel_n_elem()=" << mesh.parallel_n_elem() << std::endl;
    }

  STOP_LOG ("read()","Nemesis_IO");

  // For ParallelMesh, it seems that _is_serial is true by default.  A hack to
  // make the Mesh think it's parallel might be to call:
  mesh.delete_remote_elements();
}

#else

void Nemesis_IO::read (const std::string& )
{
  libMesh::err <<  "ERROR, Nemesis API is not defined!" << std::endl;
  libmesh_error();
}

#endif
