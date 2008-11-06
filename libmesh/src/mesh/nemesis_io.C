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
#include "nemesis_io.h"
#include "nemesis_io_helper.h"
#include "parallel_mesh.h"
#include "parallel.h"
#include "utility.h" // is_sorted

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
  // This function must be run on all processors at once
  parallel_only();
  
  if (_verbose)
    {
      std::cout << "[" << libMesh::processor_id() << "] ";
      std::cout << "Reading Nemesis file on processor: " << libMesh::processor_id() << std::endl;
    }
  
  // Construct a filename string for this processor.
  //
  // FIXME: This assumes you are reading in a mesh on exactly the
  // same number of processors it was written out on!!
  // This should be generalized at some point...
  std::ostringstream file_oss;

  file_oss << base_filename;
    
  // In the 1 processor case, assume the base_filename is the filename...
  if ( libMesh::n_processors() > 1 )
    file_oss << '.' << libMesh::n_processors() << '.' << libMesh::processor_id();

  std::cout << "Opening file: " << file_oss.str() << std::endl;

  // Make it a little easier to access the ExodusII helper object...
  ExodusII_IO_Helper& ex2helper = nemhelper.ex2helper;
  
  // Open the Exodus file
  ex2helper.open(file_oss.str().c_str());

  // Only split files will have additional Nemesis information
  if ( libMesh::n_processors() > 1 )
    {
      // Get global information: number of nodes, elems, blocks, nodesets and sidesets
      nemhelper.get_init_global();

      // Get "load balance" information.  This includes the number of internal & border
      // nodes and elements as well as the number of communication maps.
      nemhelper.get_loadbal_param();

      // Communicate n_internal_elems and n_border_elems to all processors.  This
      // determines the global element numbering.  
      // std::vector<int> num_internal_elems_all (libMesh::n_processors());
      // Parallel::allgather(nemhelper.num_internal_elems, num_internal_elems_all);
      //       if (_verbose)
      // 	{
      // 	  std::cout << "[" << libMesh::processor_id() << "] ";
      // 	  for (unsigned int j=0; j<libMesh::n_processors(); ++j)
      // 	    std::cout << "num_internal_elems_all["
      // 		      <<	j << "]="
      // 		      << num_internal_elems_all[j]
      // 		      << " ";
      // 	  std::cout << std::endl;
      // 	}

      // The call to get_loadbal_param() gets 7 pieces of information.  We allgather
      // these now across all processors to determine some global numberings...
      std::vector<int> all_loadbal_data ( 7 );
      all_loadbal_data[0] = nemhelper.num_internal_nodes;
      all_loadbal_data[1] = nemhelper.num_border_nodes;
      all_loadbal_data[2] = nemhelper.num_external_nodes;
      all_loadbal_data[3] = nemhelper.num_internal_elems;
      all_loadbal_data[4] = nemhelper.num_border_elems;
      all_loadbal_data[5] = nemhelper.num_node_cmaps;
      all_loadbal_data[6] = nemhelper.num_elem_cmaps;
      Parallel::allgather(all_loadbal_data);
      
      if (_verbose)
	{
	  // Ex: print out the number of internal nodes on all processors
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  for (unsigned int j=0,c=0; c<libMesh::n_processors(); j+=7)
	    std::cout << "num_internal_nodes["
		      << c++ << "]="
		      << all_loadbal_data[j]
		      << " ";
	  std::cout << std::endl;
	}

      
      
      // Get a reference to the ParallelMesh.  
      ParallelMesh& mesh = this->mesh();


      // The get_cmap_params() function reads in the:
      //
      // node_cmap_ids[],
      // node_cmap_node_cnts[],
      // elem_cmap_ids[],
      // elem_cmap_elem_cnts[],
      nemhelper.get_cmap_params();

      // Each processor is responsible for numbering all its internal nodes and the border nodes
      // for all adjoining processors with higher IDs.  
      unsigned int nodes_i_must_number = nemhelper.num_internal_nodes;
      for (unsigned int i=0; i<nemhelper.node_cmap_ids.size(); ++i)
	{
	  // If I have a cmap entry for a processor with a higher ID than mine...
	  if (static_cast<unsigned int>(nemhelper.node_cmap_ids[i]) > libMesh::processor_id())
	    {
	      // ... then I am responsible for numbering those nodes.
	      nodes_i_must_number += nemhelper.node_cmap_node_cnts[i];
	    }
	}

      if (_verbose)
	{
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "nodes_i_must_number=" << nodes_i_must_number << std::endl;
	}
      
      // Communicate the nodes_i_must_number information to all processors.  This
      // determines a global node numbering offset.  FIXME: This should be allgathered
      // with other data for better communication efficiency.
      std::vector<unsigned int> all_nodes_i_must_number(libMesh::n_processors());
      Parallel::allgather(nodes_i_must_number, all_nodes_i_must_number);      

      // The sum of all the entries in this vector should sum to the number of global nodes
      libmesh_assert(std::accumulate(all_nodes_i_must_number.begin(),
				     all_nodes_i_must_number.end(),
				     0) == nemhelper.num_nodes_global);

      // Compute the global_node_offsets, the amount by which to offset the local node numbering
      // on my processor.  The offset is determined by summing all counts of smaller-ID'd processors'
      // nodes they must number.
      std::vector<unsigned int> global_node_offsets (libMesh::n_processors());
      for (unsigned int i=1; i<global_node_offsets.size(); ++i)
	global_node_offsets[i] = global_node_offsets[i-1]+all_nodes_i_must_number[i-1]; 

      const unsigned int my_node_offset = global_node_offsets[libMesh::processor_id()];
      if (_verbose)
	{
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "my_node_offset="
		    << my_node_offset
		    << std::endl;
	}

      
      // Read the IDs of the interior, boundary, and external nodes.  This function
      // fills the vectors:
      // node_mapi[],
      // node_mapb[],
      // node_mape[]
      nemhelper.get_node_map();

      // Read each node communication map for this processor.  This function
      // fills the vectors of vectors named:
      // node_cmap_node_ids[][]
      // node_cmap_proc_ids[][]
      nemhelper.get_node_cmap();

      // Local information: Read the following information from the standard Exodus header
      // title[0]
      // num_dim
      // num_nodes
      // num_elem
      // num_elem_blk
      // num_node_sets
      // num_side_sets
      ex2helper.read_header();
      ex2helper.print_header();

      // Be sure number of dimensions is equal to the number of dimensions in the mesh supplied.
      libmesh_assert(static_cast<unsigned int>(ex2helper.num_dim) == mesh.mesh_dimension());

      // Read nodes from the exodus file: this fills the ex2helper.x,y,z arrays.
      ex2helper.read_nodes();

      // Add internal nodes to the ParallelMesh, using the node ID offset we computed and the current
      // processor's ID.
      for (int i=0; i<nemhelper.num_internal_nodes; ++i)
	{
	  const Real x = ex2helper.x[ nemhelper.node_mapi[i] ];
	  const Real y = ex2helper.y[ nemhelper.node_mapi[i] ];
	  const Real z = ex2helper.z[ nemhelper.node_mapi[i] ];
	  mesh.add_point (Point(x,y,z), my_node_offset+i, libMesh::processor_id());
	}

      // How many nodes are now in the mesh?
      // Forces mesh.n_nodes() to return the correct value.
      // mesh.update_parallel_id_counts(); 

      if (_verbose)
	{
	  // Report the number of nodes which have been added locally
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "mesh.n_nodes()=" << mesh.n_nodes() << std::endl;

	  // Reports the number of nodes that have been added in total.
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "mesh.parallel_n_nodes()=" << mesh.parallel_n_nodes() << std::endl;
	}
      
      /*
      // Let's check if our hypothesis about the node_cmap_node_ids is valid.  Print the
      // x,y,z location of the entries in node_cmap_node_ids.  This is hard-coded to print
      // entries in the zeroth node cmap on processors 0 and 1, since I happen to
      // know they share a boundary....
      if ((libMesh::processor_id()==0) || (libMesh::processor_id()==1))
	{
	  
// 	  // Print index and x,y,z location of [0][0] element 
// 	  unsigned int test_index = nemhelper.node_cmap_node_ids[0][0];
// 	  std::cout << "[" << libMesh::processor_id() << "] ";
// 	  std::cout << "nemhelper.node_cmap_node_ids[0][0]="
// 		    << test_index
// 		    << std::endl;

	  // Print about 10 entries to see if they match
	  for (unsigned int i=0; i<std::min(static_cast<unsigned int>(10),
					    nemhelper.node_cmap_node_ids[0].size()); ++i)
	    {
	      unsigned int test_index =
		nemhelper.node_cmap_node_ids[0][i];
	      
	      std::cout << "[" << libMesh::processor_id() << "] ";
	      std::cout << "Node " << test_index
			<< " = ("
			<< ex2helper.x[ test_index ]
			<< ", "
			<< ex2helper.y[ test_index ]
			<< ", "
			<< ex2helper.z[ test_index ]
			<< ")" << std::endl;
	    }
	}
      */


      // Create a symmetric container for storing node_cmap_node_ids.  This is necessary
      // in order to add identical nodes with the correct global ID, which is determined by
      // a lower-numbered processor.
      std::vector<std::vector<int> > symm_node_cmap_node_ids( nemhelper.node_cmap_node_ids.size() );

      // Loop over the existing node_cmap_ids, allocate space where we will need to receive data.
      for (unsigned int i=0; i<nemhelper.node_cmap_node_ids.size(); ++i)
	{
	  unsigned int other_proc = nemhelper.node_cmap_ids[i];
	  if (other_proc < libMesh::processor_id())
	    {
	      // We will need to pull data from the lower-ID processor.  Let's allocate
	      // space.
	      symm_node_cmap_node_ids[i].resize( nemhelper.node_cmap_node_ids[i].size() );

	      if (_verbose)
		{
		  std::cout << "[" << libMesh::processor_id() << "] ";
		  std::cout << "Allocating space for " 
			    << nemhelper.node_cmap_node_ids[i].size()
			    << " entries in symm_node_cmap_node_ids["<<i<<"]."
			    << " Values will come from processor "
			    << other_proc
			    << "."
			    << std::endl;
		}
	    }
	}

      // Loop again, this time initiate the communication.  We are gonna try posting
      // all the non-blocking receives followed by all the non-blocking sends later,
      // on the recommendation of Bill and Robert...
      std::vector<Parallel::request> isend_requests, irecv_requests;
      
      for (unsigned int i=0; i<nemhelper.node_cmap_node_ids.size(); ++i)
	{
	  unsigned int other_proc = nemhelper.node_cmap_ids[i];
	  if (other_proc < libMesh::processor_id())
	    {

	      // A symmetric hash of other_proc and libMesh::processor_id.  See also, elem.h.
	      // Unless there is a very large number of processors (~2045) this hash
	      // basically boils down to n0 + 32*n1, where n0 (resp. n1) is the min (resp. max)
	      // of "other_proc" and libMesh::processor_id().
	      unsigned int tag =
		(std::min(other_proc, libMesh::processor_id())     ) % 65449 +
		(std::max(other_proc, libMesh::processor_id()) << 5) % 65449;
	      
	      if (_verbose)
		std::cout << "[" << libMesh::processor_id() << "] Receiving, tag=" << tag << std::endl;
	      
	      // We need to receive data from other_proc
	      irecv_requests.push_back(Parallel::request());
	      Parallel::nonblocking_receive
                (other_proc,                 // sending proc id
		 symm_node_cmap_node_ids[i], // recv buffer
		 irecv_requests.back(),      // request
		 tag);                       // tag
	    }
	}


      for (unsigned int i=0; i<nemhelper.node_cmap_node_ids.size(); ++i)
	{
	  unsigned int other_proc = nemhelper.node_cmap_ids[i];
	  if (other_proc > libMesh::processor_id())
	    {
	      // A symmetric hash of other_proc and libMesh::processor_id.  See also, elem.h.
	      unsigned int tag =
		(std::min(other_proc, libMesh::processor_id())     ) % 65449 +
		(std::max(other_proc, libMesh::processor_id()) << 5) % 65449;

	      if (_verbose)
		std::cout << "[" << libMesh::processor_id() << "] Sending, tag=" << tag << std::endl;
		
	      // We need to send our data to other_proc
	      isend_requests.push_back(Parallel::request());
	      Parallel::nonblocking_send(other_proc,                      // destination proc id
			                 nemhelper.node_cmap_node_ids[i], // send buffer
			                 isend_requests.back(),           // request
			                 tag                              // tag
			                 );
	    }
	}

      // Wait for sends and receives to complete
      Parallel::wait (isend_requests);
      Parallel::wait (irecv_requests);
      
      // Well, did anything happen?
      if (_verbose)
	{
	  for (unsigned int i=0; i<symm_node_cmap_node_ids.size(); ++i)
	    {
	      if (symm_node_cmap_node_ids[i].size())
		{
		  std::cout << "[" << libMesh::processor_id() << "] symm_node_cmap_node_ids["<<i<<"]=";
		  for (unsigned int j=0; j<symm_node_cmap_node_ids[i].size(); ++j)
		    std::cout << symm_node_cmap_node_ids[i][j] << " ";
		  std::cout << std::endl;
		}
	    }
	}


      // Let's now add the border nodes to the Mesh using either our numbering scheme or a lower-ID
      // processor's numbering scheme if that's applicable.
      // eg. Processor zero always numbers his own border nodes.

#ifndef NDEBUG
      // We shall assume that the nemhelper.node_cmap_node_ids arrays are sorted.  This was the
      // case for our test problem, but I'm not sure it will always be the case...  This test
      // will fail in non-optimized runs if the arrays are not sorted.
      for (int i=0; i<nemhelper.num_node_cmaps; ++i)
	{
	  if (!Utility::is_sorted(nemhelper.node_cmap_node_ids[i].begin(),
				  nemhelper.node_cmap_node_ids[i].end()))
	    {
	      std::cerr << "Nemesis_IO: The border node communication map IDs are not sorted.\n"
			<< "Nemesis_IO::read() currently assumes this and can't continue."
			<< std::endl;
	      libmesh_error();
	    }
	}
#endif


      for (int i=0; i<nemhelper.num_border_nodes; ++i)
	{
	  // Caution! Exodus/Nemesis' node numbering is 1-based, so this local index is 1-based!!
	  int local_border_node_index_i = nemhelper.node_mapb[i];

	  // Find the first occurrence of local_border_node_index_i in the
	  // node_cmap_node_ids arrays.  We can use a binary search
	  // because we assumed the arrays were sorted...  We only need the
	  // first such occurrence because it will be the lowest-ID CPU which
	  // matches.
	  int
	    cmap_proc_id_index_match=-1,
	    cmap_node_id_index_match=-1;
	  for (int j=0; j<nemhelper.num_node_cmaps; ++j)
	    {
	      std::vector<int>::iterator it =
		Utility::binary_find(nemhelper.node_cmap_node_ids[j].begin(),
				     nemhelper.node_cmap_node_ids[j].end(),
				     local_border_node_index_i);
// 		std::lower_bound(nemhelper.node_cmap_node_ids[j].begin(),
// 				 nemhelper.node_cmap_node_ids[j].end(),
// 				 local_border_node_index_i);

	      // If binary_find() returns a non-end iterator, the value was found.
	      if (it != nemhelper.node_cmap_node_ids[j].end())
		{
		  cmap_proc_id_index_match = j;
		  cmap_node_id_index_match = std::distance(nemhelper.node_cmap_node_ids[j].begin(), it);
		  break; // out of for loop
		}
	    }

	  // If a shared processor was not found...
	  if (cmap_proc_id_index_match == -1)
	    {
	      std::cerr << "A matching entry in a node communication"
			<< " map was not found for local border node "
			<< local_border_node_index_i << std::endl; 
	      libmesh_error();
	    }

	  // If the shared processor has a higher ID, we use our offset...
	  unsigned int this_node_global_id = my_node_offset + local_border_node_index_i;
	  
	  // ... Otherwise, we use the shared processor's offset and numbering.
	  if (nemhelper.node_cmap_ids[cmap_proc_id_index_match] <
	      static_cast<int>(libMesh::processor_id()))
	    {
	      this_node_global_id =
		global_node_offsets[cmap_proc_id_index_match]+
		symm_node_cmap_node_ids[cmap_proc_id_index_match][cmap_node_id_index_match];
	      
	      if (_verbose)
		{
// 		  std::cout << "[" << libMesh::processor_id() << "] ";
// 		  std::cout << "Adding border node "
// 			    << local_border_node_index_i
// 			    << " with global node index "
// 			    << global_node_offsets[cmap_proc_id_index_match]
// 			    << "+"
// 			    << symm_node_cmap_node_ids[cmap_proc_id_index_match][cmap_node_id_index_match]
// 			    << "="
// 			    << this_node_global_id
// 			    << " (1-based)." << std::endl;
		}
	    }

	  else // Print statement only for case where we use our own numbering.
	    {
	      if (_verbose)
		{
// 		  std::cout << "[" << libMesh::processor_id() << "] ";
// 		  std::cout << "Adding border node "
// 			    << local_border_node_index_i
// 			    << " with global node index "
// 			    << my_node_offset
// 			    << "+"
// 			    << local_border_node_index_i
// 			    << "="
// 			    << this_node_global_id
// 			    << " (1-based)." << std::endl;
		}
	    }


	  // Extract x,y,z values using the local 1-based indexing.
	  const Real x = ex2helper.x[ local_border_node_index_i ];
	  const Real y = ex2helper.y[ local_border_node_index_i ];
	  const Real z = ex2helper.z[ local_border_node_index_i ];

	  // If this_node_global_id==0, the Exodus/Nemesis numbering scheme was
	  // not 1-based!
	  libmesh_assert(this_node_global_id > 0);
	  
	  // Finally, add the border node to the Mesh.  Don't forget: use a zero-based numbering scheme!!
	  mesh.add_point (Point(x,y,z), this_node_global_id-1, libMesh::processor_id());
	} // end for (int i=0; i<nemhelper.num_border_nodes; ++i)


      // See what the node count is up to now.
      if (_verbose)
	{
	  // Report the number of nodes which have been added locally
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "mesh.n_nodes()=" << mesh.n_nodes() << std::endl;

	  // Reports the number of nodes that have been added in total.
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "mesh.parallel_n_nodes()=" << mesh.parallel_n_nodes() << std::endl;
	}



      // --------------------------------------------------------------------------------
      // --------------------------------------------------------------------------------
      // --------------------------------------------------------------------------------

      
      // We can now read in the elements...Exodus stores them in blocks in which all
      // elements have the same geometric type.  This code is adapted directly from exodusII_io.C

      // Assertion: The sum of the border and internal elements on all processors
      // should equal nemhelper.num_elems_global
#ifndef NDEBUG
      int sum_internal_elems=0, sum_border_elems=0;
      for (unsigned int j=3,c=0; c<libMesh::n_processors(); j+=7,++c)
	sum_internal_elems += all_loadbal_data[j];

      for (unsigned int j=4,c=0; c<libMesh::n_processors(); j+=7,++c)
	sum_border_elems += all_loadbal_data[j];

      if (_verbose)
	{
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "sum_internal_elems=" << sum_internal_elems << std::endl;

	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "sum_border_elems=" << sum_border_elems << std::endl;
	}

      libmesh_assert(sum_internal_elems+sum_border_elems == nemhelper.num_elems_global);
#endif

      // Compute my_elem_offset, the amount by which to offset the local elem numbering
      // on my processor.
      unsigned int my_elem_offset = 0;
      for (unsigned int i=0; i<libMesh::processor_id(); ++i)
	my_elem_offset += (all_loadbal_data[7*i + 3]+  // num_internal_elems, proc i
			   all_loadbal_data[7*i + 4]); // num_border_elems, proc i
      std::cout << "[" << libMesh::processor_id() << "] ";
      std::cout << "my_elem_offset=" << my_elem_offset << std::endl;


      // Fills in the: 
      // global_elem_blk_ids[] and
      // global_elem_blk_cnts[] arrays.
      nemhelper.get_eb_info_global();

      // Fills in the vectors
      // elem_mapi[num_internal_elems]
      // elem_mapb[num_border_elems  ]
      // These tell which of the (locally-numbered) elements are internal and which are border elements.
      // In our test example these arrays are sorted (but non-contiguous), which makes it possible to
      // binary search for each element ID... however I don't think we need to distinguish between the
      // two types, since either can have nodes the boundary!
      nemhelper.get_elem_map();
      
      // Fills in the vectors of vectors:
      // elem_cmap_elem_ids[][]
      // elem_cmap_side_ids[][]
      // elem_cmap_proc_ids[][]
      // These arrays are of size num_elem_cmaps * elem_cmap_elem_cnts[i], i = 0..num_elem_cmaps
      nemhelper.get_elem_cmap();
      
      // Get information about the element blocks:
      // (read in the array ex2helper.block_ids[])
      ex2helper.read_block_info();

      // Local indexing offset maps block element indices to local numbering scheme. 
      int nelem_last_block = 0;

      // Instantiate the ElementMaps interface.  This is what translates LibMesh's
      // element numbering scheme to Exodus's.
      ExodusII_IO_Helper::ElementMaps em;
      

      // Read in the element connectivity for each block by
      // looping over all the blocks.
      for (int i=0; i<ex2helper.num_elem_blk; i++)
	{
	  // Read the information for block i:  For ex2helper.block_ids[i], reads
	  // elem_type
	  // num_elem_this_blk
	  // num_nodes_per_elem
	  // num_attr
	  // connect <-- the nodal connectivity array for each element in the block.
	  ex2helper.read_elem_in_block(i);

	  // Set subdomain ID based on the block ID.
	  int subdomain_id = ex2helper.block_ids[i];

	  // Create a type string (this uses the null-terminated string ctor).
	  const std::string type_str ( &(ex2helper.elem_type[0]) ); 

	  // Set any relevant node/edge maps for this element
	  const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(type_str); 

	  if (_verbose)
	    std::cout << "Reading a block of " << type_str << " elements." << std::endl;
      
	  // Loop over all the elements in this block
	  int jmax = nelem_last_block + ex2helper.num_elem_this_blk;
	  for (int j=nelem_last_block; j<jmax; j++)
	    {
	      Elem* elem = Elem::build (conv.get_canonical_type()).release();
	      libmesh_assert (elem);

	      // Assign subdomain and processor ID to the newly-created Elem.
	      // Assigning the processor ID beforehand ensures that the Elem is
	      // not added as an "unpartitioned" element.  Note that the element
	      // numbering in Exodus is also 1-based.
	      elem->subdomain_id() = subdomain_id;
	      elem->processor_id() = libMesh::processor_id();
	      elem->set_id()       = my_elem_offset + j;
		
// 	      if (_verbose)
// 		{
// 		  std::cout << "[" << libMesh::processor_id() << "] ";
// 		  std::cout << "Before mesh.add_elem(), elem->id()==" << elem->id() << std::endl;
// 		}
	      
	      // Add the created Elem to the Mesh, catch the Elem
	      // pointer that the Mesh throws back.
	      elem = mesh.add_elem (elem); 

// 	      if (_verbose)
// 		{
// 		  std::cout << "[" << libMesh::processor_id() << "] ";
// 		  std::cout << "After mesh.add_elem(), elem->id()==" << elem->id() << std::endl;
// 		}
	      
	      // Set all the nodes for this element
	      if (_verbose)
		{
		  std::cout << "[" << libMesh::processor_id() << "] ";
		  // std::cout << "Setting " << ex2helper.num_nodes_per_elem << " nodes per element." << std::endl;
		  std::cout << "Setting nodes for Elem " << elem->id() << std::endl;
		}
	      
	      for (int k=0; k<ex2helper.num_nodes_per_elem; k++)
		{
		  // global index
		  int gi = (j-nelem_last_block)*ex2helper.num_nodes_per_elem + conv.get_node_map(k); 

		  // Exodus (local) node number (1-based)
		  int node_number = ex2helper.connect[gi];

		  // Is this an internal or boundary node?  Try to find in node_mapi[] first...
		  bool node_found = false;
		  std::vector<int>::iterator it =
		    Utility::binary_find(nemhelper.node_mapi.begin(),
					 nemhelper.node_mapi.end(),
					 node_number);
// 		    std::lower_bound(nemhelper.node_mapi.begin(),
// 				     nemhelper.node_mapi.end(),
// 				     node_number);

		  // If binary_find returns a non-end iterator, the value was found!
		  if (it != nemhelper.node_mapi.end())
		    {
		      // Node found was an internal node, use the Exodus number plus
		      // the current processor's node offset.
		      node_number += my_node_offset;

		      if (_verbose)
			{
			  std::cout << "[" << libMesh::processor_id() << "] ";
			  std::cout << "Elem " << elem->id() << ": setting internal node k=" << k
				    << " to " << node_number << std::endl;
			}

		      // We can stop searching
		      node_found = true;
		    }

		  if (!node_found)
		    {
		      // Keep looking for the node ID in the boundary node list.
// 		      it = lower_bound(nemhelper.node_mapb.begin(),
// 				       nemhelper.node_mapb.end(),
// 				       node_number);
		      it = Utility::binary_find(nemhelper.node_mapb.begin(),
						nemhelper.node_mapb.end(),
						node_number);

		      // If binary_find() returns a non-end iterator, the value was found!
		      if (it != nemhelper.node_mapb.end())
			{
			  node_found = true;
			      
			  // Node found was a boundary node.  Determine its offset into
			  // the node_cmap_node_ids vectors and then determine if we will
			  // number it or if we will get its ID from another processor.
			  // Note: this code is copied from above and should be pulled out as
			  // a subroutine...
			  int
			    cmap_proc_id_index_match=-1,
			    cmap_node_id_index_match=-1;
			  for (int m=0; m<nemhelper.num_node_cmaps; ++m)
			    {
			      std::vector<int>::iterator it2 =
				Utility::binary_find(nemhelper.node_cmap_node_ids[m].begin(),
						     nemhelper.node_cmap_node_ids[m].end(),
						     node_number);
// 				std::lower_bound(nemhelper.node_cmap_node_ids[m].begin(),
// 						 nemhelper.node_cmap_node_ids[m].end(),
// 						 node_number);

			      // If binary_find() returns a non-end iterator, the value was found!
			      if (it2 != nemhelper.node_cmap_node_ids[m].end())
				{
				  cmap_proc_id_index_match = m;
				  cmap_node_id_index_match = std::distance(nemhelper.node_cmap_node_ids[m].begin(), it2);
				  break; // out of for (m) loop
				}
			    } // end for (int m=0; m<nemhelper.num_node_cmaps; ++m)

			      // If a shared processor was not found...
			  if (cmap_proc_id_index_match == -1)
			    {
			      std::cerr << "A matching entry in a node communication"
					<< " map was not found for local border node "
					<< node_number << std::endl; 
			      libmesh_error();
			    }

			      
			  // If the shared processor has a lower ID, it numbers the node...
			  if (nemhelper.node_cmap_ids[cmap_proc_id_index_match] <
			      static_cast<int>(libMesh::processor_id()))
			    {
			      node_number =
				global_node_offsets[cmap_proc_id_index_match]+
				symm_node_cmap_node_ids[cmap_proc_id_index_match][cmap_node_id_index_match];
			    }

			  // ... Otherwise, we use our own offset.
			  else
			    {
			      node_number += my_node_offset;
			    }

			  if (_verbose)
			    {
			      std::cout << "[" << libMesh::processor_id() << "] ";
			      std::cout << "Elem " << elem->id() << ": setting boundary node k=" << k
					<< " to " << node_number << std::endl;
			    }
			  
			} // end if (it != nemhelper.node_mapb.end())
		    } // end if (!node_found)
		  
		  // If the node ID was not found in either the internal or boundary node ID
		  // lists, we can't continue!
		  if (!node_found)
		    {
		      std::cerr << "Could not correctly set node k=" << k << " for element " << elem->id() << std::endl;
		      libmesh_error();
		    }
		  
		  
		  // Set node number, subtract 1 since exodus is internally 1-based
		  elem->set_node(k) = mesh.node_ptr((node_number-1)); 
		}
	    }
      
	  // Maintain running sum of # of elements per block.
	  // Should equal total number of elements (for this processor) in the end.
	  nelem_last_block += ex2helper.num_elem_this_blk;
	} // end for (int i=0; i<ex2helper.num_elem_blk; i++)


      // See what the elem count is up to now.
      if (_verbose)
	{
	  // Report the number of elements which have been added locally
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "mesh.n_elem()=" << mesh.n_elem() << std::endl;

	  // Reports the number of elements that have been added in total.
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "mesh.parallel_n_elem()=" << mesh.parallel_n_elem() << std::endl;
	}

      // For ParallelMesh, it seems that _is_serial is true by default.  A hack to
      // make the Mesh think it's parallel might be to call:
      mesh.delete_remote_elements();
      
    } // end if ( libMesh::n_processors() > 1 )

  else
    {
      // Do something different when running in serial??
      libmesh_error();
    }

  
  
}

#else

void Nemesis_IO::read (const std::string& )
{
  std::cerr <<  "ERROR, Nemesis API is not defined!" << std::endl;
  libmesh_error();
}

#endif
