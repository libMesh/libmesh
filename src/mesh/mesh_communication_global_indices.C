// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_hilbert.h"
#include "libmesh/parallel_sort.h"
#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/node_range.h"
#ifdef LIBMESH_HAVE_LIBHILBERT
#  include "hilbert.h"
#endif

#ifdef LIBMESH_HAVE_LIBHILBERT
namespace { // anonymous namespace for helper functions

  using namespace libMesh;

  // Utility function to map (x,y,z) in [bbox.min, bbox.max]^3 into
  // [0,max_inttype]^3 for computing Hilbert keys
  void get_hilbert_coords (const Point &p,
			   const MeshTools::BoundingBox &bbox,
			   CFixBitVec icoords[3])
  {
    static const Hilbert::inttype max_inttype = static_cast<Hilbert::inttype>(-1);

    const long double // put (x,y,z) in [0,1]^3 (don't divide by 0)
      x = ((bbox.first(0) == bbox.second(0)) ? 0. :
	   (p(0)-bbox.first(0))/(bbox.second(0)-bbox.first(0))),

#if LIBMESH_DIM > 1
      y = ((bbox.first(1) == bbox.second(1)) ? 0. :
	   (p(1)-bbox.first(1))/(bbox.second(1)-bbox.first(1))),
#else
      y = 0.,
#endif

#if LIBMESH_DIM > 2
      z = ((bbox.first(2) == bbox.second(2)) ? 0. :
	   (p(2)-bbox.first(2))/(bbox.second(2)-bbox.first(2)));
#else
      z = 0.;
#endif

    // (icoords) in [0,max_inttype]^3
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
      dof_id_type pos = range.first_idx();
      for (ConstNodeRange::const_iterator it = range.begin(); it!=range.end(); ++it)
	{
	  const Node* node = (*it);
	  libmesh_assert(node);
	  libmesh_assert_less (pos, _keys.size());
	  _keys[pos++] = get_hilbert_index (*node, _bbox);
	}
    }

    // computes the hilbert index for an element
    void operator() (const ConstElemRange &range) const
    {
      dof_id_type pos = range.first_idx();
      for (ConstElemRange::const_iterator it = range.begin(); it!=range.end(); ++it)
	{
	  const Elem* elem = (*it);
	  libmesh_assert(elem);
	  libmesh_assert_less (pos, _keys.size());
	  _keys[pos++] = get_hilbert_index (elem->centroid(), _bbox);
	}
    }

  private:
    const MeshTools::BoundingBox &_bbox;
    std::vector<Hilbert::HilbertIndices> &_keys;
  };
}
#endif


namespace libMesh
{


// ------------------------------------------------------------
// MeshCommunication class members
#if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)
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

  const Parallel::Communicator &communicator (mesh.communicator());

  // Global bounding box
  MeshTools::BoundingBox bbox =
    MeshTools::bounding_box (mesh);

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

#if 0
    // It's O(N^2) to check that these keys don't duplicate before the
    // sort...
      MeshBase::const_node_iterator nodei = mesh.local_nodes_begin();
      for (std::size_t i = 0; i != node_keys.size(); ++i, ++nodei)
        {
          MeshBase::const_node_iterator nodej = mesh.local_nodes_begin();
          for (std::size_t j = 0; j != i; ++j, ++nodej)
            {
              if (node_keys[i] == node_keys[j])
                {
		  std::cerr << "Error: nodes with duplicate Hilbert keys!" <<
                    std::endl;
                  CFixBitVec icoords[3], jcoords[3];
                  get_hilbert_coords(**nodej, bbox, jcoords);
		  std::cerr <<
		    "node " << (*nodej)->id() << ", " <<
		    *(Point*)(*nodej) << " has HilbertIndices " <<
		    node_keys[j] << std::endl;
                  get_hilbert_coords(**nodei, bbox, icoords);
		  std::cerr <<
		    "node " << (*nodei)->id() << ", " <<
		    *(Point*)(*nodei) << " has HilbertIndices " <<
		    node_keys[i] << std::endl;
                  libmesh_error();
                }
            }
        }
#endif
    }

    // Elements next
    {
      ConstElemRange er (mesh.local_elements_begin(),
			 mesh.local_elements_end());
      elem_keys.resize (er.size());
      Threads::parallel_for (er, ComputeHilbertKeys (bbox, elem_keys));

#if 0
    // For elements, the keys can be (and in the case of TRI, are
    // expected to be) duplicates, but only if the elements are at
    // different levels
      MeshBase::const_element_iterator elemi = mesh.local_elements_begin();
      for (std::size_t i = 0; i != elem_keys.size(); ++i, ++elemi)
        {
          MeshBase::const_element_iterator elemj = mesh.local_elements_begin();
          for (std::size_t j = 0; j != i; ++j, ++elemj)
            {
              if ((elem_keys[i] == elem_keys[j]) &&
                  ((*elemi)->level() == (*elemj)->level()))
                {
		  std::cerr <<
                    "Error: level " << (*elemi)->level() <<
                    " elements with duplicate Hilbert keys!" <<
                    std::endl;
		  std::cerr <<
		    "level " << (*elemj)->level() << " elem\n" <<
                    (**elemj) << " centroid " <<
		    (*elemj)->centroid() << " has HilbertIndices " <<
		    elem_keys[j] << " or " <<
		    get_hilbert_index((*elemj)->centroid(), bbox) <<
                    std::endl;
		  std::cerr <<
		    "level " << (*elemi)->level() << " elem\n" <<
                    (**elemi) << " centroid " <<
		    (*elemi)->centroid() << " has HilbertIndices " <<
		    elem_keys[i] << " or " <<
		    get_hilbert_index((*elemi)->centroid(), bbox) <<
                    std::endl;
                  libmesh_error();
                }
            }
        }
#endif
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
    node_upper_bounds(communicator.size()),
    elem_upper_bounds(communicator.size());

  { // limit scope of temporaries
    std::vector<Hilbert::HilbertIndices> recvbuf(2*communicator.size());
    std::vector<unsigned short int> /* do not use a vector of bools here since it is not always so! */
      empty_nodes (communicator.size()),
      empty_elem  (communicator.size());
    std::vector<Hilbert::HilbertIndices> my_max(2);

    communicator.allgather (static_cast<unsigned short int>(my_node_bin.empty()), empty_nodes);
    communicator.allgather (static_cast<unsigned short int>(my_elem_bin.empty()),  empty_elem);

    if (!my_node_bin.empty()) my_max[0] = my_node_bin.back();
    if (!my_elem_bin.empty()) my_max[1] = my_elem_bin.back();

    communicator.allgather (my_max, /* identical_buffer_sizes = */ true);

    // Be cereful here.  The *_upper_bounds will be used to find the processor
    // a given object belongs to.  So, if a processor contains no objects (possible!)
    // then copy the bound from the lower processor id.
    for (processor_id_type p=0; p<communicator.size(); p++)
      {
	node_upper_bounds[p] = my_max[2*p+0];
	elem_upper_bounds[p] = my_max[2*p+1];

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
	requested_ids (communicator.size());
      // Results to gather from each processor
      std::vector<std::vector<dof_id_type> >
	filled_request (communicator.size());

      {
      MeshBase::const_node_iterator       it  = mesh.nodes_begin();
      const MeshBase::const_node_iterator end = mesh.nodes_end();

      // build up list of requests
      for (; it != end; ++it)
	{
	  const Node* node = (*it);
	  libmesh_assert(node);
	  const Hilbert::HilbertIndices hi =
	    get_hilbert_index (*node, bbox);
	  const processor_id_type pid =
	    libmesh_cast_int<processor_id_type>
	    (std::distance (node_upper_bounds.begin(),
			    std::lower_bound(node_upper_bounds.begin(),
					     node_upper_bounds.end(),
					     hi)));

	  libmesh_assert_less (pid, communicator.size());

	  requested_ids[pid].push_back(hi);
	}
      }

      // The number of objects in my_node_bin on each processor
      std::vector<dof_id_type> node_bin_sizes(communicator.size());
      communicator.allgather (static_cast<dof_id_type>(my_node_bin.size()), node_bin_sizes);

      // The offset of my first global index
      dof_id_type my_offset = 0;
      for (processor_id_type pid=0; pid<communicator.rank(); pid++)
	my_offset += node_bin_sizes[pid];

      // start with pid=0, so that we will trade with ourself
      for (processor_id_type pid=0; pid<communicator.size(); pid++)
	{
          // Trade my requests with processor procup and procdown
          const processor_id_type procup = (communicator.rank() + pid) %
                                            communicator.size();
          const processor_id_type procdown = (communicator.size() +
                                              communicator.rank() - pid) %
                                              communicator.size();

          std::vector<Hilbert::HilbertIndices> request_to_fill;
          communicator.send_receive(procup, requested_ids[procup],
				    procdown, request_to_fill);

	  // Fill the requests
	  std::vector<dof_id_type> global_ids; /**/ global_ids.reserve(request_to_fill.size());
	  for (std::size_t idx=0; idx<request_to_fill.size(); idx++)
	    {
	      const Hilbert::HilbertIndices &hi = request_to_fill[idx];
	      libmesh_assert_less_equal (hi, node_upper_bounds[communicator.rank()]);

	      // find the requested index in my node bin
	      std::vector<Hilbert::HilbertIndices>::const_iterator pos =
		 std::lower_bound (my_node_bin.begin(), my_node_bin.end(), hi);
	      libmesh_assert (pos != my_node_bin.end());
	      libmesh_assert_equal_to (*pos, hi);

	      // Finally, assign the global index based off the position of the index
	      // in my array, properly offset.
	      global_ids.push_back (std::distance(my_node_bin.begin(), pos) + my_offset);
	    }

	  // and trade back
	  communicator.send_receive (procdown, global_ids,
				     procup,   filled_request[procup]);
	}

      // We now have all the filled requests, so we can loop through our
      // nodes once and assign the global index to each one.
      {
	std::vector<std::vector<dof_id_type>::const_iterator>
	  next_obj_on_proc; next_obj_on_proc.reserve(communicator.size());
	for (processor_id_type pid=0; pid<communicator.size(); pid++)
	  next_obj_on_proc.push_back(filled_request[pid].begin());

        {
	MeshBase::node_iterator       it  = mesh.nodes_begin();
	const MeshBase::node_iterator end = mesh.nodes_end();

	for (; it != end; ++it)
	  {
	    Node* node = (*it);
	    libmesh_assert(node);
	    const Hilbert::HilbertIndices hi =
	      get_hilbert_index (*node, bbox);
	    const processor_id_type pid =
	      libmesh_cast_int<processor_id_type>
	      (std::distance (node_upper_bounds.begin(),
			      std::lower_bound(node_upper_bounds.begin(),
					       node_upper_bounds.end(),
					       hi)));

	    libmesh_assert_less (pid, communicator.size());
	    libmesh_assert (next_obj_on_proc[pid] != filled_request[pid].end());

	    const dof_id_type global_index = *next_obj_on_proc[pid];
	    libmesh_assert_less (global_index, mesh.n_nodes());
	    node->set_id() = global_index;

	    ++next_obj_on_proc[pid];
	  }
        }
      }
    }

    //---------------------------------------------------
    // elements next -- all elements, not just local ones
    {
      // Request sets to send to each processor
      std::vector<std::vector<Hilbert::HilbertIndices> >
	requested_ids (communicator.size());
      // Results to gather from each processor
      std::vector<std::vector<dof_id_type> >
	filled_request (communicator.size());

      {
      MeshBase::const_element_iterator       it  = mesh.elements_begin();
      const MeshBase::const_element_iterator end = mesh.elements_end();

      for (; it != end; ++it)
	{
	  const Elem* elem = (*it);
	  libmesh_assert(elem);
	  const Hilbert::HilbertIndices hi =
	    get_hilbert_index (elem->centroid(), bbox);
	  const processor_id_type pid =
	    libmesh_cast_int<processor_id_type>
	    (std::distance (elem_upper_bounds.begin(),
			    std::lower_bound(elem_upper_bounds.begin(),
					     elem_upper_bounds.end(),
					     hi)));

	  libmesh_assert_less (pid, communicator.size());

	  requested_ids[pid].push_back(hi);
	}
      }

      // The number of objects in my_elem_bin on each processor
      std::vector<dof_id_type> elem_bin_sizes(communicator.size());
      communicator.allgather (static_cast<dof_id_type>(my_elem_bin.size()), elem_bin_sizes);

      // The offset of my first global index
      dof_id_type my_offset = 0;
      for (processor_id_type pid=0; pid<communicator.rank(); pid++)
	my_offset += elem_bin_sizes[pid];

      // start with pid=0, so that we will trade with ourself
      for (processor_id_type pid=0; pid<communicator.size(); pid++)
	{
          // Trade my requests with processor procup and procdown
          const processor_id_type procup = (communicator.rank() + pid) %
                                            communicator.size();
          const processor_id_type procdown = (communicator.size() +
                                              communicator.rank() - pid) %
                                              communicator.size();

          std::vector<Hilbert::HilbertIndices> request_to_fill;
          communicator.send_receive(procup, requested_ids[procup],
				    procdown, request_to_fill);

	  // Fill the requests
	  std::vector<dof_id_type> global_ids; /**/ global_ids.reserve(request_to_fill.size());
	  for (std::size_t idx=0; idx<request_to_fill.size(); idx++)
	    {
	      const Hilbert::HilbertIndices &hi = request_to_fill[idx];
	      libmesh_assert_less_equal (hi, elem_upper_bounds[communicator.rank()]);

	      // find the requested index in my elem bin
	      std::vector<Hilbert::HilbertIndices>::const_iterator pos =
		std::lower_bound (my_elem_bin.begin(), my_elem_bin.end(), hi);
	      libmesh_assert (pos != my_elem_bin.end());
	      libmesh_assert_equal_to (*pos, hi);

	      // Finally, assign the global index based off the position of the index
	      // in my array, properly offset.
	      global_ids.push_back (std::distance(my_elem_bin.begin(), pos) + my_offset);
	    }

	  // and trade back
	  communicator.send_receive (procdown, global_ids,
				     procup,   filled_request[procup]);
	}

      // We now have all the filled requests, so we can loop through our
      // elements once and assign the global index to each one.
      {
	std::vector<std::vector<dof_id_type>::const_iterator>
	  next_obj_on_proc; next_obj_on_proc.reserve(communicator.size());
	for (processor_id_type pid=0; pid<communicator.size(); pid++)
	  next_obj_on_proc.push_back(filled_request[pid].begin());

        {
	MeshBase::element_iterator       it  = mesh.elements_begin();
	const MeshBase::element_iterator end = mesh.elements_end();

	for (; it != end; ++it)
	  {
	    Elem* elem = (*it);
	    libmesh_assert(elem);
	    const Hilbert::HilbertIndices hi =
	      get_hilbert_index (elem->centroid(), bbox);
	    const processor_id_type pid =
	      libmesh_cast_int<processor_id_type>
	      (std::distance (elem_upper_bounds.begin(),
			      std::lower_bound(elem_upper_bounds.begin(),
					       elem_upper_bounds.end(),
					       hi)));

	    libmesh_assert_less (pid, communicator.size());
	    libmesh_assert (next_obj_on_proc[pid] != filled_request[pid].end());

	    const dof_id_type global_index = *next_obj_on_proc[pid];
	    libmesh_assert_less (global_index, mesh.n_elem());
	    elem->set_id() = global_index;

	    ++next_obj_on_proc[pid];
	  }
        }
      }
    }
  }

  STOP_LOG ("assign_global_indices()", "MeshCommunication");
}
#else // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI
void MeshCommunication::assign_global_indices (MeshBase&) const
{
}
#endif // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI



#if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)
template <typename ForwardIterator>
void MeshCommunication::find_global_indices (const Parallel::Communicator &communicator,
					     const MeshTools::BoundingBox &bbox,
					     const ForwardIterator &begin,
					     const ForwardIterator &end,
					     std::vector<dof_id_type> &index_map) const
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
  index_map.reserve(std::distance (begin, end));

  //-------------------------------------------------------------
  // (1) compute Hilbert keys
  // These aren't trivial to compute, and we will need them again.
  // But the binsort will sort the input vector, trashing the order
  // that we'd like to rely on.  So, two vectors...
  std::vector<Hilbert::HilbertIndices>
    sorted_hilbert_keys,
    hilbert_keys;
  sorted_hilbert_keys.reserve(index_map.capacity());
  hilbert_keys.reserve(index_map.capacity());
  {
    START_LOG("compute_hilbert_indices()", "MeshCommunication");
    for (ForwardIterator it=begin; it!=end; ++it)
      {
	const Hilbert::HilbertIndices hi(get_hilbert_index (*it, bbox));
	hilbert_keys.push_back(hi);

	if ((*it)->processor_id() == communicator.rank())
	  sorted_hilbert_keys.push_back(hi);

	// someone needs to take care of unpartitioned objects!
	if ((communicator.rank() == 0) &&
	    ((*it)->processor_id() == DofObject::invalid_processor_id))
	  sorted_hilbert_keys.push_back(hi);
      }
    STOP_LOG("compute_hilbert_indices()", "MeshCommunication");
  }

  //-------------------------------------------------------------
  // (2) parallel sort the Hilbert keys
  START_LOG ("parallel_sort()", "MeshCommunication");
  Parallel::Sort<Hilbert::HilbertIndices> sorter (sorted_hilbert_keys);
  sorter.sort();
  STOP_LOG ("parallel_sort()", "MeshCommunication");
  const std::vector<Hilbert::HilbertIndices> &my_bin = sorter.bin();

  // The number of objects in my_bin on each processor
  std::vector<unsigned int> bin_sizes(communicator.size());
  communicator.allgather (static_cast<unsigned int>(my_bin.size()), bin_sizes);

  // The offset of my first global index
  unsigned int my_offset = 0;
  for (unsigned int pid=0; pid<communicator.rank(); pid++)
    my_offset += bin_sizes[pid];

  //-------------------------------------------------------------
  // (3) get the max value on each processor
  std::vector<Hilbert::HilbertIndices>
    upper_bounds(1);

  if (!my_bin.empty())
    upper_bounds[0] = my_bin.back();

  communicator.allgather (upper_bounds, /* identical_buffer_sizes = */ true);

  // Be cereful here.  The *_upper_bounds will be used to find the processor
  // a given object belongs to.  So, if a processor contains no objects (possible!)
  // then copy the bound from the lower processor id.
  for (unsigned int p=1; p<communicator.size(); p++)
    if (!bin_sizes[p]) upper_bounds[p] = upper_bounds[p-1];


  //-------------------------------------------------------------
  // (4) determine the position in the global ranking for
  //     each local object
  {
    //----------------------------------------------
    // all objects, not just local ones

    // Request sets to send to each processor
    std::vector<std::vector<Hilbert::HilbertIndices> >
      requested_ids (communicator.size());
    // Results to gather from each processor
    std::vector<std::vector<unsigned int> >
      filled_request (communicator.size());

    // build up list of requests
    std::vector<Hilbert::HilbertIndices>::const_iterator hi =
      hilbert_keys.begin();

    for (ForwardIterator it = begin; it != end; ++it)
      {
	libmesh_assert (hi != hilbert_keys.end());
	const processor_id_type pid =
	  libmesh_cast_int<processor_id_type>
	  (std::distance (upper_bounds.begin(),
			  std::lower_bound(upper_bounds.begin(),
					   upper_bounds.end(),
					   *hi)));

	libmesh_assert_less (pid, communicator.size());

	requested_ids[pid].push_back(*hi);

	hi++;
	// go ahead and put pid in index_map, that way we
	// don't have to repeat the std::lower_bound()
	index_map.push_back(pid);
      }

    // start with pid=0, so that we will trade with ourself
    std::vector<Hilbert::HilbertIndices> request_to_fill;
    std::vector<unsigned int> global_ids;
    for (unsigned int pid=0; pid<communicator.size(); pid++)
      {
	// Trade my requests with processor procup and procdown
	const unsigned int procup = (communicator.rank() + pid) %
                                     communicator.size();
	const unsigned int procdown = (communicator.size() +
                                       communicator.rank() - pid) %
                                       communicator.size();

	communicator.send_receive(procup, requested_ids[procup],
				  procdown, request_to_fill);

	// Fill the requests
	global_ids.clear(); /**/ global_ids.reserve(request_to_fill.size());
	for (unsigned int idx=0; idx<request_to_fill.size(); idx++)
	  {
	    const Hilbert::HilbertIndices &hilbert_indices = request_to_fill[idx];
	    libmesh_assert_less_equal (hilbert_indices, upper_bounds[communicator.rank()]);

	    // find the requested index in my node bin
	    std::vector<Hilbert::HilbertIndices>::const_iterator pos =
	      std::lower_bound (my_bin.begin(), my_bin.end(), hilbert_indices);
	    libmesh_assert (pos != my_bin.end());
	    libmesh_assert_equal_to (*pos, hilbert_indices);

	    // Finally, assign the global index based off the position of the index
	    // in my array, properly offset.
	    global_ids.push_back (std::distance(my_bin.begin(), pos) + my_offset);
	  }

	// and trade back
	communicator.send_receive (procdown, global_ids,
				   procup,   filled_request[procup]);
      }

    // We now have all the filled requests, so we can loop through our
    // nodes once and assign the global index to each one.
    {
      std::vector<std::vector<unsigned int>::const_iterator>
	next_obj_on_proc; next_obj_on_proc.reserve(communicator.size());
      for (unsigned int pid=0; pid<communicator.size(); pid++)
	next_obj_on_proc.push_back(filled_request[pid].begin());

      unsigned int cnt=0;
      for (ForwardIterator it = begin; it != end; ++it, cnt++)
	{
	  const unsigned int pid = index_map[cnt];

	  libmesh_assert_less (pid, communicator.size());
	  libmesh_assert (next_obj_on_proc[pid] != filled_request[pid].end());

	  const unsigned int global_index = *next_obj_on_proc[pid];
	  index_map[cnt] = global_index;

	  ++next_obj_on_proc[pid];
	}
    }
  }

  STOP_LOG ("find_global_indices()", "MeshCommunication");
}
#else // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI
template <typename ForwardIterator>
void MeshCommunication::find_global_indices (const MeshTools::BoundingBox &,
					     const ForwardIterator &begin,
					     const ForwardIterator &end,
					     std::vector<unsigned int> &index_map) const
{
  index_map.clear();
  index_map.reserve(std::distance (begin, end));
  unsigned int index = 0;
  for (ForwardIterator it=begin; it!=end; ++it)
    index_map.push_back(index++);
}
#endif // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI



//------------------------------------------------------------------
template void MeshCommunication::find_global_indices<MeshBase::const_node_iterator> (const Parallel::Communicator &,
										     const MeshTools::BoundingBox &,
										     const MeshBase::const_node_iterator &,
										     const MeshBase::const_node_iterator &,
										     std::vector<dof_id_type> &) const;

template void MeshCommunication::find_global_indices<MeshBase::const_element_iterator> (const Parallel::Communicator &,
											const MeshTools::BoundingBox &,
											const MeshBase::const_element_iterator &,
											const MeshBase::const_element_iterator &,
											std::vector<dof_id_type> &) const;
template void MeshCommunication::find_global_indices<MeshBase::node_iterator> (const Parallel::Communicator &,
									       const MeshTools::BoundingBox &,
									       const MeshBase::node_iterator &,
									       const MeshBase::node_iterator &,
									       std::vector<dof_id_type> &) const;

template void MeshCommunication::find_global_indices<MeshBase::element_iterator> (const Parallel::Communicator &,
										  const MeshTools::BoundingBox &,
										  const MeshBase::element_iterator &,
										  const MeshBase::element_iterator &,
										  std::vector<dof_id_type> &) const;

} // namespace libMesh
