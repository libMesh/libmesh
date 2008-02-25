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
#include "mesh_tools.h"
#include "mesh_communication.h"
#include "parallel.h"
#include "parallel_sort.h"
#include "elem.h"
#include "elem_range.h"
#include "node_range.h"
#ifdef HAVE_LIBHILBERT
#  include "hilbert.h"
#endif
 
#ifdef HAVE_LIBHILBERT 
namespace { // anonymous namespace for helper functions

  // Utility function to map (x,y,z) in [bbox.min, bbox.max]^3 into
  // [0,max_inttype]^3 for computing Hilbert keys
  void get_hilbert_coords (const Point &p,
			   const MeshTools::BoundingBox &bbox,
			   CFixBitVec icoords[3])
  {
    static const Hilbert::inttype max_inttype = static_cast<Hilbert::inttype>(-1);

    const double // put (x,y,z) in [0,1]^3 (don't divide by 0)
      x = ((bbox.first(0) == bbox.second(0)) ? 0. :
	   (p(0)-bbox.first(0))/(bbox.second(0)-bbox.first(0))),
	  
      y = ((bbox.first(1) == bbox.second(1)) ? 0. :
	   (p(1)-bbox.first(1))/(bbox.second(1)-bbox.first(1))),
	  
      z = ((bbox.first(2) == bbox.second(2)) ? 0. :
	   (p(2)-bbox.first(2))/(bbox.second(2)-bbox.first(2)));
	
    // (iccords) in [0,max_inttype]^3
    icoords[0] = static_cast<Hilbert::inttype>(x*max_inttype);
    icoords[1] = static_cast<Hilbert::inttype>(y*max_inttype);
    icoords[2] = static_cast<Hilbert::inttype>(z*max_inttype);
  }


  
  // Compute the hilbert index
  template <typename T>
  Hilbert::HilbertIndices
  get_hilbert_index (const T *p,
		     const MeshTools::BoundingBox &bbox)
  {
    static const unsigned int sizeof_inttype = sizeof(Hilbert::inttype);

    Hilbert::HilbertIndices index;
    CFixBitVec icoords[3];
    Hilbert::BitVecType bv;
    get_hilbert_coords (*p, bbox, icoords);
    Hilbert::coordsToIndex (icoords, 8*sizeof_inttype, 3, bv);
    index = bv;

    return index;
  }

  template <>
  Hilbert::HilbertIndices
  get_hilbert_index (const Elem *e,
		     const MeshTools::BoundingBox &bbox)
  {
    static const unsigned int sizeof_inttype = sizeof(Hilbert::inttype);

    Hilbert::HilbertIndices index;
    CFixBitVec icoords[3];
    Hilbert::BitVecType bv;
    get_hilbert_coords (e->centroid(), bbox, icoords);
    Hilbert::coordsToIndex (icoords, 8*sizeof_inttype, 3, bv);
    index = bv;

    return index;
  }


  
  // Compute the hilbert index
  Hilbert::HilbertIndices
  get_hilbert_index (const Point &p,
		     const MeshTools::BoundingBox &bbox)
  {
    static const unsigned int sizeof_inttype = sizeof(Hilbert::inttype);

    Hilbert::HilbertIndices index;
    CFixBitVec icoords[3];
    Hilbert::BitVecType bv;
    get_hilbert_coords (p, bbox, icoords);
    Hilbert::coordsToIndex (icoords, 8*sizeof_inttype, 3, bv);
    index = bv;

    return index;
  }

  // Helper class for threaded Hilbert key computation
  class ComputeHilbertKeys 
  {
  public:
    ComputeHilbertKeys (const MeshTools::BoundingBox &bbox,
			std::vector<Hilbert::HilbertIndices> &keys) :
      _bbox(bbox),
      _keys(keys)
    {}
							
    // computes the hilbert index for a node
    void operator() (const ConstNodeRange &range) const
    {
      unsigned int pos = range.first_idx();
      for (ConstNodeRange::const_iterator it = range.begin(); it!=range.end(); ++it)
	{
	  const Node* node = (*it);
	  assert (node != NULL);
	  assert (pos < _keys.size());
	  _keys[pos++] = get_hilbert_index (*node, _bbox);
	}	
    }

    // computes the hilbert index for an element
    void operator() (const ConstElemRange &range) const
    { 
      unsigned int pos = range.first_idx();
      for (ConstElemRange::const_iterator it = range.begin(); it!=range.end(); ++it)
	{
	  const Elem* elem = (*it);
	  assert (elem != NULL);
	  assert (pos < _keys.size());
	  _keys[pos++] = get_hilbert_index (elem->centroid(), _bbox);
	}	
    }

  private:					
    const MeshTools::BoundingBox &_bbox;
    std::vector<Hilbert::HilbertIndices> &_keys;
  };

}
#endif



// ------------------------------------------------------------
// MeshCommunication class members
#if defined(HAVE_LIBHILBERT) && defined(HAVE_MPI)
void MeshCommunication::assign_global_indices (MeshBase& mesh) const
{
  START_LOG ("assign_global_indices()", "MeshCommunication");

  // This method determines partition-agnostic global indices
  // for nodes and elements.

  // Algorithm:
  // (1) compute the Hilbert key for each local node/element
  // (2) perform a parallel sort of the Hilbert key
  // (3) get the min/max value on each processor
  // (4) determine the position in the global ranking for
  //     each local object

  // Global bounding box
  MeshTools::BoundingBox bbox =
    MeshTools::bounding_box (mesh);

  // Set up a derived MPI datatype to handle communication of HilbertIndices
  MPI_Datatype hilbert_type;
  MPI_Type_contiguous (3, MPI_UNSIGNED, &hilbert_type);
  MPI_Type_commit     (&hilbert_type);


  //-------------------------------------------------------------
  // (1) compute Hilbert keys
  std::vector<Hilbert::HilbertIndices>
    node_keys, elem_keys;
  
  {
    // Nodes first
    {
      ConstNodeRange nr (mesh.local_nodes_begin(),
			 mesh.local_nodes_end());
      node_keys.resize (nr.size()); 
      Threads::parallel_for (nr, ComputeHilbertKeys (bbox, node_keys));
    }
    
    // Elements next
    {
      ConstElemRange er (mesh.local_elements_begin(),
			 mesh.local_elements_end());
      elem_keys.resize (er.size());
      Threads::parallel_for (er, ComputeHilbertKeys (bbox, elem_keys));
    }    
  } // done computing Hilbert keys

  
 
  //-------------------------------------------------------------
  // (2) parallel sort the Hilbert keys
  Parallel::Sort<Hilbert::HilbertIndices> node_sorter (node_keys);
  node_sorter.sort(); /* done with node_keys */ //node_keys.clear();
    
  const std::vector<Hilbert::HilbertIndices> &my_node_bin = 
    node_sorter.bin();

  Parallel::Sort<Hilbert::HilbertIndices> elem_sorter (elem_keys);
  elem_sorter.sort(); /* done with elem_keys */ //elem_keys.clear();
    
  const std::vector<Hilbert::HilbertIndices> &my_elem_bin = 
    elem_sorter.bin();



  //-------------------------------------------------------------
  // (3) get the max value on each processor
  std::vector<Hilbert::HilbertIndices>    
    node_upper_bounds(libMesh::n_processors()),
    elem_upper_bounds(libMesh::n_processors());

  { // limit scope of temporaries
    std::vector<Hilbert::HilbertIndices> recvbuf(2*libMesh::n_processors());
    std::vector<unsigned short int> /* do not use a vector of bools here since it is not always so! */
      empty_nodes (libMesh::n_processors()),
      empty_elem  (libMesh::n_processors());
    Hilbert::HilbertIndices my_max[2];
    
    Parallel::allgather (static_cast<unsigned short int>(my_node_bin.empty()), empty_nodes);
    Parallel::allgather (static_cast<unsigned short int>(my_elem_bin.empty()),  empty_elem);
     	       
    if (!my_node_bin.empty()) my_max[0] = my_node_bin.back();
    if (!my_elem_bin.empty()) my_max[1] = my_elem_bin.back();

    MPI_Allgather (my_max,      2, hilbert_type,
		   &recvbuf[0], 2, hilbert_type,
		   libMesh::COMM_WORLD);

    // Be cereful here.  The *_upper_bounds will be used to find the processor
    // a given object belongs to.  So, if a processor contains no objects (possible!)
    // then copy the bound from the lower processor id.
    for (unsigned int p=0; p<libMesh::n_processors(); p++)
      {
	node_upper_bounds[p] = recvbuf[2*p+0];
	elem_upper_bounds[p] = recvbuf[2*p+1];

	if (p > 0) // default hilbert index value is the OK upper bound for processor 0.
	  {
	    if (empty_nodes[p]) node_upper_bounds[p] = node_upper_bounds[p-1];
	    if (empty_elem[p])  elem_upper_bounds[p] = elem_upper_bounds[p-1];
	  }
      }
  }



  //-------------------------------------------------------------
  // (4) determine the position in the global ranking for
  //     each local object
  {
    //----------------------------------------------
    // Nodes first -- all nodes, not just local ones
    {
      // Request sets to send to each processor
      std::vector<std::vector<Hilbert::HilbertIndices> > 
	requested_ids (libMesh::n_processors());
      // Results to gather from each processor
      std::vector<std::vector<unsigned int> >
	filled_request (libMesh::n_processors());

      MeshBase::const_node_iterator       it  = mesh.nodes_begin();
      const MeshBase::const_node_iterator end = mesh.nodes_end();

      // build up list of requests
      for (; it != end; ++it)
	{
	  const Node* node = (*it);
	  assert (node != NULL);
	  const Hilbert::HilbertIndices hi = 
	    get_hilbert_index (*node, bbox);
	  const unsigned int pid = 
	    std::distance (node_upper_bounds.begin(), 
			   std::lower_bound(node_upper_bounds.begin(), 
					    node_upper_bounds.end(),
					    hi));

	  assert (pid < libMesh::n_processors());

	  requested_ids[pid].push_back(hi);
	}

      // The number of objects in my_node_bin on each processor
      std::vector<unsigned int> node_bin_sizes(libMesh::n_processors());
      Parallel::allgather (static_cast<unsigned int>(my_node_bin.size()), node_bin_sizes);

      // The offset of my first global index
      unsigned int my_offset = 0;
      for (unsigned int pid=0; pid<libMesh::processor_id(); pid++)
	my_offset += node_bin_sizes[pid];

      // start with pid=0, so that we will trade with ourself
      for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	{
          // Trade my requests with processor procup and procdown
          const unsigned int procup = (libMesh::processor_id() + pid) %
                                       libMesh::n_processors();
          const unsigned int procdown = (libMesh::n_processors() +
                                         libMesh::processor_id() - pid) %
                                         libMesh::n_processors();

          std::vector<Hilbert::HilbertIndices> request_to_fill;
          Parallel::send_receive(procup, requested_ids[procup],
                                 procdown, request_to_fill,
				 hilbert_type);	  

	  // Fill the requests
	  std::vector<unsigned int> global_ids; /**/ global_ids.reserve(request_to_fill.size());
	  for (unsigned int idx=0; idx<request_to_fill.size(); idx++)
	    {
	      const Hilbert::HilbertIndices &hi = request_to_fill[idx];
	      assert (hi <= node_upper_bounds[libMesh::processor_id()]);
	      
	      // find the requested index in my node bin
	      std::vector<Hilbert::HilbertIndices>::const_iterator pos =
		 std::lower_bound (my_node_bin.begin(), my_node_bin.end(), hi);
	      assert (pos != my_node_bin.end());
	      assert (*pos == hi);
	      
	      // Finally, assign the global index based off the position of the index
	      // in my array, properly offset.
	      global_ids.push_back (std::distance(my_node_bin.begin(), pos) + my_offset);
	    }

	  // and trade back
	  Parallel::send_receive (procdown, global_ids,
				  procup,   filled_request[procup]);
	}

      // We now have all the filled requests, so we can loop through our
      // nodes once and assign the global index to each one.
      {
	std::vector<std::vector<unsigned int>::const_iterator>
	  next_obj_on_proc; next_obj_on_proc.reserve(libMesh::n_processors());
	for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	  next_obj_on_proc.push_back(filled_request[pid].begin());

	MeshBase::node_iterator       it  = mesh.nodes_begin();
	const MeshBase::node_iterator end = mesh.nodes_end();

	for (; it != end; ++it)
	  {
	    Node* node = (*it);
	    assert (node != NULL);
	    const Hilbert::HilbertIndices hi = 
	      get_hilbert_index (*node, bbox);
	    const unsigned int pid = 
	      std::distance (node_upper_bounds.begin(), 
			     std::lower_bound(node_upper_bounds.begin(), 
					      node_upper_bounds.end(),
					      hi));

	    assert (pid < libMesh::n_processors());
	    assert (next_obj_on_proc[pid] != filled_request[pid].end());

	    const unsigned int global_index = *next_obj_on_proc[pid];
	    assert (global_index < mesh.n_nodes());
	    node->set_id() = global_index;
	    
	    ++next_obj_on_proc[pid];
	  }
      }
    }

    //---------------------------------------------------
    // elements next -- all elements, not just local ones
    {
      // Request sets to send to each processor
      std::vector<std::vector<Hilbert::HilbertIndices> > 
	requested_ids (libMesh::n_processors());
      // Results to gather from each processor
      std::vector<std::vector<unsigned int> >
	filled_request (libMesh::n_processors());
      
      MeshBase::const_element_iterator       it  = mesh.elements_begin();
      const MeshBase::const_element_iterator end = mesh.elements_end();

      for (; it != end; ++it)
	{
	  const Elem* elem = (*it);
	  assert (elem != NULL);
	  const Hilbert::HilbertIndices hi = 
	    get_hilbert_index (elem->centroid(), bbox);
	  const unsigned int pid = 
	    std::distance (elem_upper_bounds.begin(), 
			   std::lower_bound(elem_upper_bounds.begin(), 
					    elem_upper_bounds.end(),
					    hi));

	  assert (pid < libMesh::n_processors());

	  requested_ids[pid].push_back(hi);
	}

      // The number of objects in my_elem_bin on each processor
      std::vector<unsigned int> elem_bin_sizes(libMesh::n_processors());
      Parallel::allgather (static_cast<unsigned int>(my_elem_bin.size()), elem_bin_sizes);

      // The offset of my first global index
      unsigned int my_offset = 0;
      for (unsigned int pid=0; pid<libMesh::processor_id(); pid++)
	my_offset += elem_bin_sizes[pid];

      // start with pid=0, so that we will trade with ourself
      for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	{
          // Trade my requests with processor procup and procdown
          const unsigned int procup = (libMesh::processor_id() + pid) %
                                       libMesh::n_processors();
          const unsigned int procdown = (libMesh::n_processors() +
                                         libMesh::processor_id() - pid) %
                                         libMesh::n_processors();

          std::vector<Hilbert::HilbertIndices> request_to_fill;
          Parallel::send_receive(procup, requested_ids[procup],
                                 procdown, request_to_fill,
				 hilbert_type);	  

	  // Fill the requests
	  std::vector<unsigned int> global_ids; /**/ global_ids.reserve(request_to_fill.size());
	  for (unsigned int idx=0; idx<request_to_fill.size(); idx++)
	    {
	      const Hilbert::HilbertIndices &hi = request_to_fill[idx];
	      assert (hi <= elem_upper_bounds[libMesh::processor_id()]);
	      
	      // find the requested index in my elem bin
	      std::vector<Hilbert::HilbertIndices>::const_iterator pos =
		std::lower_bound (my_elem_bin.begin(), my_elem_bin.end(), hi);
	      assert (pos != my_elem_bin.end());
	      assert (*pos == hi);
	      
	      // Finally, assign the global index based off the position of the index
	      // in my array, properly offset.
	      global_ids.push_back (std::distance(my_elem_bin.begin(), pos) + my_offset);
	    }

	  // and trade back
	  Parallel::send_receive (procdown, global_ids,
				  procup,   filled_request[procup]);
	}

      // We now have all the filled requests, so we can loop through our
      // elements once and assign the global index to each one.
      {
	std::vector<std::vector<unsigned int>::const_iterator>
	  next_obj_on_proc; next_obj_on_proc.reserve(libMesh::n_processors());
	for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	  next_obj_on_proc.push_back(filled_request[pid].begin());

	MeshBase::element_iterator       it  = mesh.elements_begin();
	const MeshBase::element_iterator end = mesh.elements_end();

	for (; it != end; ++it)
	  {
	    Elem* elem = (*it);
	    assert (elem != NULL);
	    const Hilbert::HilbertIndices hi = 
	      get_hilbert_index (elem->centroid(), bbox);
	    const unsigned int pid = 
	      std::distance (elem_upper_bounds.begin(), 
			     std::lower_bound(elem_upper_bounds.begin(), 
					      elem_upper_bounds.end(),
					      hi));

	    assert (pid < libMesh::n_processors());
	    assert (next_obj_on_proc[pid] != filled_request[pid].end());

	    const unsigned int global_index = *next_obj_on_proc[pid];
	    assert (global_index < mesh.n_elem());
	    elem->set_id() = global_index;
	    
	    ++next_obj_on_proc[pid];
	  }
      }
    }        
  }


  // Clean up
  MPI_Type_free (&hilbert_type);

  STOP_LOG ("assign_global_indices()", "MeshCommunication");
}
#else // HAVE_LIBHILBERT, HAVE_MPI
void MeshCommunication::assign_global_indices (MeshBase&) const
{
}
#endif // HAVE_LIBHILBERT, HAVE_MPI



#if defined(HAVE_LIBHILBERT) && defined(HAVE_MPI)
template <typename ForwardIterator>
void MeshCommunication::find_global_indices (const MeshTools::BoundingBox &bbox,
					     const ForwardIterator &begin,
					     const ForwardIterator &end,
					     std::vector<unsigned int> &index_map) const
{
  START_LOG ("find_global_indices()", "MeshCommunication");

  // This method determines partition-agnostic global indices
  // for nodes and elements.

  // Algorithm:
  // (1) compute the Hilbert key for each local node/element
  // (2) perform a parallel sort of the Hilbert key
  // (3) get the min/max value on each processor
  // (4) determine the position in the global ranking for
  //     each local object

  index_map.clear();
  
  // Set up a derived MPI datatype to handle communication of HilbertIndices
  MPI_Datatype hilbert_type;
  MPI_Type_contiguous (3, MPI_UNSIGNED, &hilbert_type);
  MPI_Type_commit     (&hilbert_type);


  //-------------------------------------------------------------
  // (1) compute Hilbert keys
  std::vector<Hilbert::HilbertIndices> hilbert_keys;
  {
    for (ForwardIterator it=begin; it!=end; ++it)
      if ((*it)->processor_id() == libMesh::processor_id())
	hilbert_keys.push_back(get_hilbert_index (*it, bbox));
    
    // someone needs to take care of unpartitioned objects!
    if (libMesh::processor_id() == 0)
      for (ForwardIterator it=begin; it!=end; ++it)
	if ((*it)->processor_id() == DofObject::invalid_processor_id)
	  hilbert_keys.push_back(get_hilbert_index (*it, bbox));
  }

  //-------------------------------------------------------------
  // (2) parallel sort the Hilbert keys
  Parallel::Sort<Hilbert::HilbertIndices> sorter (hilbert_keys);
  sorter.sort(); 
    
  const std::vector<Hilbert::HilbertIndices> &my_bin = sorter.bin();
  
  //-------------------------------------------------------------
  // (3) get the max value on each processor
  std::vector<Hilbert::HilbertIndices>    
    upper_bounds(libMesh::n_processors());

  { // limit scope of temporaries
    std::vector<Hilbert::HilbertIndices> recvbuf(libMesh::n_processors());
    std::vector<unsigned short int> /* do not use a vector of bools here since it is not always so! */
      empty_bin (libMesh::n_processors());
    Hilbert::HilbertIndices my_max;
    
    Parallel::allgather (static_cast<unsigned short int>(my_bin.empty()), empty_bin);
     	       
    if (!my_bin.empty()) my_max = my_bin.back();

    MPI_Allgather (&my_max,      1, hilbert_type,
		   &recvbuf[0], 1, hilbert_type,
		   libMesh::COMM_WORLD);

    // Be cereful here.  The *_upper_bounds will be used to find the processor
    // a given object belongs to.  So, if a processor contains no objects (possible!)
    // then copy the bound from the lower processor id.
    for (unsigned int p=0; p<libMesh::n_processors(); p++)
      {
	upper_bounds[p] = recvbuf[p];

	if (p > 0) // default hilbert index value is the OK upper bound for processor 0.
	  if (empty_bin[p]) upper_bounds[p] = upper_bounds[p-1];
      }
  }



  //-------------------------------------------------------------
  // (4) determine the position in the global ranking for
  //     each local object
  {
    //----------------------------------------------
    // all objects, not just local ones
    
    // Request sets to send to each processor
    std::vector<std::vector<Hilbert::HilbertIndices> > 
      requested_ids (libMesh::n_processors());
    // Results to gather from each processor
    std::vector<std::vector<unsigned int> >
      filled_request (libMesh::n_processors());

    // build up list of requests
    for (ForwardIterator it = begin; it != end; ++it)
      {
	const Hilbert::HilbertIndices hi = 
	  get_hilbert_index (*it, bbox);
	const unsigned int pid = 
	  std::distance (upper_bounds.begin(), 
			 std::lower_bound(upper_bounds.begin(), 
					  upper_bounds.end(),
					  hi));

	assert (pid < libMesh::n_processors());

	requested_ids[pid].push_back(hi);
      }

    // The number of objects in my_bin on each processor
    std::vector<unsigned int> bin_sizes(libMesh::n_processors());
    Parallel::allgather (static_cast<unsigned int>(my_bin.size()), bin_sizes);
    
    // The offset of my first global index
    unsigned int my_offset = 0;
    for (unsigned int pid=0; pid<libMesh::processor_id(); pid++)
      my_offset += bin_sizes[pid];

    // start with pid=0, so that we will trade with ourself
    for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
      {
	// Trade my requests with processor procup and procdown
	const unsigned int procup = (libMesh::processor_id() + pid) %
                                     libMesh::n_processors();
	const unsigned int procdown = (libMesh::n_processors() +
                                       libMesh::processor_id() - pid) %
                                       libMesh::n_processors();

	std::vector<Hilbert::HilbertIndices> request_to_fill;
	Parallel::send_receive(procup, requested_ids[procup],
			       procdown, request_to_fill,
			       hilbert_type);	  

	// Fill the requests
	std::vector<unsigned int> global_ids; /**/ global_ids.reserve(request_to_fill.size());
	for (unsigned int idx=0; idx<request_to_fill.size(); idx++)
	  {
	    const Hilbert::HilbertIndices &hi = request_to_fill[idx];
	    assert (hi <= upper_bounds[libMesh::processor_id()]);
	    
	    // find the requested index in my node bin
	    std::vector<Hilbert::HilbertIndices>::const_iterator pos =
	      std::lower_bound (my_bin.begin(), my_bin.end(), hi);
	    assert (pos != my_bin.end());
	    assert (*pos == hi);
	    
	    // Finally, assign the global index based off the position of the index
	    // in my array, properly offset.
	    global_ids.push_back (std::distance(my_bin.begin(), pos) + my_offset);
	  }
	
	// and trade back
	Parallel::send_receive (procdown, global_ids,
				procup,   filled_request[procup]);
      }

    // We now have all the filled requests, so we can loop through our
    // nodes once and assign the global index to each one.
    {
      std::vector<std::vector<unsigned int>::const_iterator>
	next_obj_on_proc; next_obj_on_proc.reserve(libMesh::n_processors());
      for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	next_obj_on_proc.push_back(filled_request[pid].begin());
      
      for (ForwardIterator it = begin; it != end; ++it)
	{
	  const Hilbert::HilbertIndices hi = 
	    get_hilbert_index (*it, bbox);
	  const unsigned int pid = 
	    std::distance (upper_bounds.begin(), 
			   std::lower_bound(upper_bounds.begin(), 
					    upper_bounds.end(),
					    hi));
	  
	  assert (pid < libMesh::n_processors());
	  assert (next_obj_on_proc[pid] != filled_request[pid].end());
	  
	  const unsigned int global_index = *next_obj_on_proc[pid];
	  index_map.push_back(global_index);
	  
	  ++next_obj_on_proc[pid];
	}
    }
  }


  // Clean up
  MPI_Type_free (&hilbert_type);

  STOP_LOG ("find_global_indices()", "MeshCommunication");
}
#else // HAVE_LIBHILBERT, HAVE_MPI
template <typename ForwardIterator>
void MeshCommunication::find_global_indices (const MeshTools::BoundingBox &,
					     const ForwardIterator &,
					     const ForwardIterator &,
					     std::vector<unsigned int> &) const
{
}
#endif // HAVE_LIBHILBERT, HAVE_MPI



//------------------------------------------------------------------
template void MeshCommunication::find_global_indices<MeshBase::const_node_iterator> (const MeshTools::BoundingBox &,
										     const MeshBase::const_node_iterator &,
										     const MeshBase::const_node_iterator &,
										     std::vector<unsigned int> &) const;

template void MeshCommunication::find_global_indices<MeshBase::const_element_iterator> (const MeshTools::BoundingBox &,
											const MeshBase::const_element_iterator &,
											const MeshBase::const_element_iterator &,
											std::vector<unsigned int> &) const;
template void MeshCommunication::find_global_indices<MeshBase::node_iterator> (const MeshTools::BoundingBox &,
										     const MeshBase::node_iterator &,
										     const MeshBase::node_iterator &,
										     std::vector<unsigned int> &) const;

template void MeshCommunication::find_global_indices<MeshBase::element_iterator> (const MeshTools::BoundingBox &,
											const MeshBase::element_iterator &,
											const MeshBase::element_iterator &,
											std::vector<unsigned int> &) const;
