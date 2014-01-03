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
#include "libmesh/mesh_base.h"
#include "libmesh/metis_partitioner.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/error_vector.h"
#include "libmesh/vectormap.h"
#include "libmesh/metis_csr_graph.h"

#ifdef LIBMESH_HAVE_METIS
// MIPSPro 7.4.2 gets confused about these nested namespaces
# ifdef __sgi
#  include <cstdarg>
# endif
  namespace Metis {
    extern "C" {
#     include "libmesh/ignore_warnings.h"
#     include "metis.h"
#     include "libmesh/restore_warnings.h"
    }
  }
#else
#  include "libmesh/sfc_partitioner.h"
#endif




namespace libMesh
{


// ------------------------------------------------------------
// MetisPartitioner implementation
void MetisPartitioner::_do_partition (MeshBase& mesh,
				      const unsigned int n_pieces)
{
  libmesh_assert_greater (n_pieces, 0);
  libmesh_assert (mesh.is_serial());

  // Check for an easy return
  if (n_pieces == 1)
    {
      this->single_partition (mesh);
      return;
    }

// What to do if the Metis library IS NOT present
#ifndef LIBMESH_HAVE_METIS

  libmesh_here();
  libMesh::err << "ERROR: The library has been built without"    << std::endl
	        << "Metis support.  Using a space-filling curve"  << std::endl
	        << "partitioner instead!"                         << std::endl;

  SFCPartitioner sfcp;

  sfcp.partition (mesh, n_pieces);

// What to do if the Metis library IS present
#else

  START_LOG("partition()", "MetisPartitioner");

  const dof_id_type n_active_elem = mesh.n_active_elem();

  // build the graph
  // std::vector<int> options(5);
  std::vector<int> vwgt(n_active_elem);
  std::vector<int> part(n_active_elem);

  int
    n = static_cast<int>(n_active_elem),  // number of "nodes" (elements)
                                          //   in the graph
//    wgtflag = 2,                          // weights on vertices only,
//                                          //   none on edges
//    numflag = 0,                          // C-style 0-based numbering
    nparts  = static_cast<int>(n_pieces), // number of subdomains to create
    edgecut = 0;                          // the numbers of edges cut by the
                                          //   resulting partition

  // Set the options
  // options[0] = 0; // use default options

  // Metis will only consider the active elements.
  // We need to map the active element ids into a
  // contiguous range.  Further, we want the unique range indexing to be
  // independednt of the element ordering, otherwise a circular dependency
  // can result in which the partitioning depends on the ordering which
  // depends on the partitioning...
  vectormap<dof_id_type, dof_id_type> global_index_map;
  global_index_map.reserve (n_active_elem);

  {
    std::vector<dof_id_type> global_index;

    MeshBase::element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::element_iterator end = mesh.active_elements_end();

    MeshCommunication().find_global_indices (mesh.comm(),
					     MeshTools::bounding_box(mesh),
					     it, end, global_index);

    libmesh_assert_equal_to (global_index.size(), n_active_elem);

    for (std::size_t cnt=0; it != end; ++it)
      {
	const Elem *elem = *it;

	global_index_map.insert (std::make_pair(elem->id(), global_index[cnt++]));
      }
    libmesh_assert_equal_to (global_index_map.size(), n_active_elem);
  }


  // Invoke METIS, but only on processor 0.
  // Then broadcast the resulting decomposition
  if (mesh.processor_id() == 0)
    {
      METIS_CSR_Graph csr_graph;

      csr_graph.offsets.resize(n_active_elem+1, 0);

      // Local scope for these
      {
	// build the graph in CSR format.  Note that
	// the edges in the graph will correspond to
	// face neighbors
	std::vector<const Elem*> neighbors_offspring;

	MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
	const MeshBase::element_iterator elem_end = mesh.active_elements_end();

	std::size_t graph_size=0;

	// (1) first pass - get the row sizes for each element by counting the number
	// of face neighbors.  Also populate the vwght array if necessary
	for (; elem_it != elem_end; ++elem_it)
	  {
	    const Elem* elem = *elem_it;

	    const dof_id_type elem_global_index =
	      global_index_map[elem->id()];

	    libmesh_assert_less (elem_global_index, vwgt.size());

	    // maybe there is a better weight?
	    // The weight is used to define what a balanced graph is
	    if(!_weights)
	      vwgt[elem_global_index] = elem->n_nodes();
	    else
	      vwgt[elem_global_index] = static_cast<int>((*_weights)[elem->id()]);

	    unsigned int num_neighbors = 0;

	    // Loop over the element's neighbors.  An element
	    // adjacency corresponds to a face neighbor
	    for (unsigned int ms=0; ms<elem->n_neighbors(); ms++)
	      {
		const Elem* neighbor = elem->neighbor(ms);

		if (neighbor != NULL)
		  {
		    // If the neighbor is active treat it
		    // as a connection
		    if (neighbor->active())
		      num_neighbors++;

#ifdef LIBMESH_ENABLE_AMR

		    // Otherwise we need to find all of the
		    // neighbor's children that are connected to
		    // us and add them
		    else
		      {
			// The side of the neighbor to which
			// we are connected
			const unsigned int ns =
			  neighbor->which_neighbor_am_i (elem);
			libmesh_assert_less (ns, neighbor->n_neighbors());

			// Get all the active children (& grandchildren, etc...)
			// of the neighbor.
			neighbor->active_family_tree (neighbors_offspring);

			// Get all the neighbor's children that
			// live on that side and are thus connected
			// to us
			for (unsigned int nc=0; nc<neighbors_offspring.size(); nc++)
			  {
			    const Elem* child =
			      neighbors_offspring[nc];

			    // This does not assume a level-1 mesh.
			    // Note that since children have sides numbered
			    // coincident with the parent then this is a sufficient test.
			    if (child->neighbor(ns) == elem)
			      {
				libmesh_assert (child->active());
				num_neighbors++;
			      }
			  }
		      }

#endif /* ifdef LIBMESH_ENABLE_AMR */

		  }
	      }

	    csr_graph.prep_n_nonzeros(elem_global_index, num_neighbors);
	    graph_size += num_neighbors;
	  }

	csr_graph.prepare_for_use();

	// (2) second pass - fill the compressed adjacency array
	elem_it  = mesh.active_elements_begin();

	for (; elem_it != elem_end; ++elem_it)
	  {
	    const Elem* elem = *elem_it;

	    const dof_id_type elem_global_index =
	      global_index_map[elem->id()];

	    unsigned int connection=0;

	    // Loop over the element's neighbors.  An element
	    // adjacency corresponds to a face neighbor
	    for (unsigned int ms=0; ms<elem->n_neighbors(); ms++)
	      {
		const Elem* neighbor = elem->neighbor(ms);

		if (neighbor != NULL)
		  {
		    // If the neighbor is active treat it
		    // as a connection
		    if (neighbor->active())
		      csr_graph(elem_global_index, connection++) = global_index_map[neighbor->id()];

#ifdef LIBMESH_ENABLE_AMR

		    // Otherwise we need to find all of the
		    // neighbor's children that are connected to
		    // us and add them
		    else
		      {
			// The side of the neighbor to which
			// we are connected
			const unsigned int ns =
			  neighbor->which_neighbor_am_i (elem);
			libmesh_assert_less (ns, neighbor->n_neighbors());

			// Get all the active children (& grandchildren, etc...)
			// of the neighbor.
			neighbor->active_family_tree (neighbors_offspring);

			// Get all the neighbor's children that
			// live on that side and are thus connected
			// to us
			for (unsigned int nc=0; nc<neighbors_offspring.size(); nc++)
			  {
			    const Elem* child =
			      neighbors_offspring[nc];

			    // This does not assume a level-1 mesh.
			    // Note that since children have sides numbered
			    // coincident with the parent then this is a sufficient test.
			    if (child->neighbor(ns) == elem)
			      {
				libmesh_assert (child->active());

				csr_graph(elem_global_index, connection++) = global_index_map[child->id()];
			      }
			  }
		      }

#endif /* ifdef LIBMESH_ENABLE_AMR */

		  }
	      }
	  }

        // We create a non-empty vals for a disconnected graph, to
        // work around a segfault from METIS.
	libmesh_assert_equal_to (csr_graph.vals.size(),
                                 std::max(graph_size,std::size_t(1)));
      } // done building the graph

      int ncon = 1;

      // Select which type of partitioning to create

      // Use recursive if the number of partitions is less than or equal to 8
      if (n_pieces <= 8)
	Metis::METIS_PartGraphRecursive(&n, &ncon, &csr_graph.offsets[0], &csr_graph.vals[0], &vwgt[0], NULL,
					NULL, &nparts, NULL, NULL, NULL,
					&edgecut, &part[0]);

      // Otherwise  use kway
      else
	Metis::METIS_PartGraphKway(&n, &ncon, &csr_graph.offsets[0], &csr_graph.vals[0], &vwgt[0], NULL,
				   NULL, &nparts, NULL, NULL, NULL,
				   &edgecut, &part[0]);

    } // end processor 0 part

  // Broadcase the resutling partition
  mesh.comm().broadcast(part);

  // Assign the returned processor ids.  The part array contains
  // the processor id for each active element, but in terms of
  // the contiguous indexing we defined above
  {
    MeshBase::element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::element_iterator end = mesh.active_elements_end();

    for (; it!=end; ++it)
      {
	Elem* elem = *it;

	libmesh_assert (global_index_map.count(elem->id()));

	const dof_id_type elem_global_index =
	  global_index_map[elem->id()];

	libmesh_assert_less (elem_global_index, part.size());
	const processor_id_type elem_procid =
	  static_cast<processor_id_type>(part[elem_global_index]);

	elem->processor_id() = elem_procid;
      }
  }

  STOP_LOG("partition()", "MetisPartitioner");
#endif
}

} // namespace libMesh
