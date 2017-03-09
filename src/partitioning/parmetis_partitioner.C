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



// Local Includes
#include "libmesh/libmesh_config.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"    // also includes mpi.h
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/parmetis_partitioner.h"
#include "libmesh/metis_partitioner.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/parmetis_helper.h"

// Include the ParMETIS header file.
#ifdef LIBMESH_HAVE_PARMETIS

namespace Parmetis {
extern "C" {
#     include "libmesh/ignore_warnings.h"
#     include "parmetis.h"
#     include "libmesh/restore_warnings.h"
}
}

#endif


// Hash maps for interior->boundary element lookups
#include LIBMESH_INCLUDE_UNORDERED_MULTIMAP
#include LIBMESH_INCLUDE_HASH
LIBMESH_DEFINE_HASH_POINTERS


namespace libMesh
{

// Minimum elements on each processor required for us to choose
// Parmetis over Metis.
#ifdef LIBMESH_HAVE_PARMETIS
const unsigned int MIN_ELEM_PER_PROC = 4;
#endif

// ------------------------------------------------------------
// ParmetisPartitioner implementation
ParmetisPartitioner::ParmetisPartitioner()
#ifdef LIBMESH_HAVE_PARMETIS
  :  _pmetis(new ParmetisHelper)
#endif
{}



ParmetisPartitioner::~ParmetisPartitioner()
{
}



void ParmetisPartitioner::_do_partition (MeshBase & mesh,
                                         const unsigned int n_sbdmns)
{
  this->_do_repartition (mesh, n_sbdmns);
}



void ParmetisPartitioner::_do_repartition (MeshBase & mesh,
                                           const unsigned int n_sbdmns)
{
  libmesh_assert_greater (n_sbdmns, 0);

  // Check for an easy return
  if (n_sbdmns == 1)
    {
      this->single_partition(mesh);
      return;
    }

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  // What to do if the Parmetis library IS NOT present
#ifndef LIBMESH_HAVE_PARMETIS

  libmesh_here();
  libMesh::err << "ERROR: The library has been built without" << std::endl
               << "Parmetis support.  Using a Metis"          << std::endl
               << "partitioner instead!"                      << std::endl;

  MetisPartitioner mp;

  mp.partition (mesh, n_sbdmns);

  // What to do if the Parmetis library IS present
#else

  // Revert to METIS on one processor.
  if (mesh.n_processors() == 1)
    {
      MetisPartitioner mp;
      mp.partition (mesh, n_sbdmns);
      return;
    }

  LOG_SCOPE("repartition()", "ParmetisPartitioner");

  // Initialize the data structures required by ParMETIS
  this->initialize (mesh, n_sbdmns);

  // Make sure all processors have enough active local elements.
  // Parmetis tends to crash when it's given only a couple elements
  // per partition.
  {
    bool all_have_enough_elements = true;
    for (std::size_t pid=0; pid<_n_active_elem_on_proc.size(); pid++)
      if (_n_active_elem_on_proc[pid] < MIN_ELEM_PER_PROC)
        all_have_enough_elements = false;

    // Parmetis will not work unless each processor has some
    // elements. Specifically, it will abort when passed a NULL
    // partition array on *any* of the processors.
    if (!all_have_enough_elements)
      {
        // FIXME: revert to METIS, although this requires a serial mesh
        MeshSerializer serialize(mesh);
        MetisPartitioner mp;
        mp.partition (mesh, n_sbdmns);
        return;
      }
  }

  // build the graph corresponding to the mesh
  this->build_graph (mesh);


  // Partition the graph
  std::vector<Parmetis::idx_t> vsize(_pmetis->vwgt.size(), 1);
  Parmetis::real_t itr = 1000000.0;
  MPI_Comm mpi_comm = mesh.comm().get();

  // Call the ParMETIS adaptive repartitioning method.  This respects the
  // original partitioning when computing the new partitioning so as to
  // minimize the required data redistribution.
  Parmetis::ParMETIS_V3_AdaptiveRepart(_pmetis->vtxdist.empty() ? libmesh_nullptr : &_pmetis->vtxdist[0],
                                       _pmetis->xadj.empty()    ? libmesh_nullptr : &_pmetis->xadj[0],
                                       _pmetis->adjncy.empty()  ? libmesh_nullptr : &_pmetis->adjncy[0],
                                       _pmetis->vwgt.empty()    ? libmesh_nullptr : &_pmetis->vwgt[0],
                                       vsize.empty()            ? libmesh_nullptr : &vsize[0],
                                       libmesh_nullptr,
                                       &_pmetis->wgtflag,
                                       &_pmetis->numflag,
                                       &_pmetis->ncon,
                                       &_pmetis->nparts,
                                       _pmetis->tpwgts.empty()  ? libmesh_nullptr : &_pmetis->tpwgts[0],
                                       _pmetis->ubvec.empty()   ? libmesh_nullptr : &_pmetis->ubvec[0],
                                       &itr,
                                       &_pmetis->options[0],
                                       &_pmetis->edgecut,
                                       _pmetis->part.empty()    ? libmesh_nullptr : &_pmetis->part[0],
                                       &mpi_comm);

  // Assign the returned processor ids
  this->assign_partitioning (mesh);

#endif // #ifndef LIBMESH_HAVE_PARMETIS ... else ...

}



// Only need to compile these methods if ParMETIS is present
#ifdef LIBMESH_HAVE_PARMETIS

void ParmetisPartitioner::initialize (const MeshBase & mesh,
                                      const unsigned int n_sbdmns)
{
  LOG_SCOPE("initialize()", "ParmetisPartitioner");

  const dof_id_type n_active_local_elem = mesh.n_active_local_elem();

  // Set parameters.
  _pmetis->wgtflag = 2;                                      // weights on vertices only
  _pmetis->ncon    = 1;                                      // one weight per vertex
  _pmetis->numflag = 0;                                      // C-style 0-based numbering
  _pmetis->nparts  = static_cast<Parmetis::idx_t>(n_sbdmns); // number of subdomains to create
  _pmetis->edgecut = 0;                                      // the numbers of edges cut by the
                                                             // partition

  // Initialize data structures for ParMETIS
  _pmetis->vtxdist.assign (mesh.n_processors()+1, 0);
  _pmetis->tpwgts.assign  (_pmetis->nparts, 1./_pmetis->nparts);
  _pmetis->ubvec.assign   (_pmetis->ncon, 1.05);
  _pmetis->part.assign    (n_active_local_elem, 0);
  _pmetis->options.resize (5);
  _pmetis->vwgt.resize    (n_active_local_elem);

  // Set the options
  _pmetis->options[0] = 1;  // don't use default options
  _pmetis->options[1] = 0;  // default (level of timing)
  _pmetis->options[2] = 15; // random seed (default)
  _pmetis->options[3] = 2;  // processor distribution and subdomain distribution are decoupled

  // Find the number of active elements on each processor.  We cannot use
  // mesh.n_active_elem_on_proc(pid) since that only returns the number of
  // elements assigned to pid which are currently stored on the calling
  // processor. This will not in general be correct for parallel meshes
  // when (pid!=mesh.processor_id()).
  _n_active_elem_on_proc.resize(mesh.n_processors());
  mesh.comm().allgather(n_active_local_elem, _n_active_elem_on_proc);

  // count the total number of active elements in the mesh.  Note we cannot
  // use mesh.n_active_elem() in general since this only returns the number
  // of active elements which are stored on the calling processor.
  // We should not use n_active_elem for any allocation because that will
  // be inheritly unscalable, but it can be useful for libmesh_assertions.
  dof_id_type n_active_elem=0;

  // Set up the vtxdist array.  This will be the same on each processor.
  // ***** Consult the Parmetis documentation. *****
  libmesh_assert_equal_to (_pmetis->vtxdist.size(),
                           cast_int<std::size_t>(mesh.n_processors()+1));
  libmesh_assert_equal_to (_pmetis->vtxdist[0], 0);

  for (processor_id_type pid=0; pid<mesh.n_processors(); pid++)
    {
      _pmetis->vtxdist[pid+1] = _pmetis->vtxdist[pid] + _n_active_elem_on_proc[pid];
      n_active_elem += _n_active_elem_on_proc[pid];
    }
  libmesh_assert_equal_to (_pmetis->vtxdist.back(), static_cast<Parmetis::idx_t>(n_active_elem));

  // ParMetis expects the elements to be numbered in contiguous blocks
  // by processor, i.e. [0, ne0), [ne0, ne0+ne1), ...
  // Since we only partition active elements we should have no expectation
  // that we currently have such a distribution.  So we need to create it.
  // Also, at the same time we are going to map all the active elements into a globally
  // unique range [0,n_active_elem) which is *independent* of the current partitioning.
  // This can be fed to ParMetis as the initial partitioning of the subdomains (decoupled
  // from the partitioning of the objects themselves).  This allows us to get the same
  // resultant partitioning independed of the input partitioning.
  BoundingBox bbox =
    MeshTools::bounding_box(mesh);

  _global_index_by_pid_map.clear();

  // Maps active element ids into a contiguous range independent of partitioning.
  // (only needs local scope)
  vectormap<dof_id_type, dof_id_type> global_index_map;

  {
    std::vector<dof_id_type> global_index;

    // create the mapping which is contiguous by processor
    dof_id_type pid_offset=0;
    for (processor_id_type pid=0; pid<mesh.n_processors(); pid++)
      {
        MeshBase::const_element_iterator       it  = mesh.active_pid_elements_begin(pid);
        const MeshBase::const_element_iterator end = mesh.active_pid_elements_end(pid);

        // note that we may not have all (or any!) the active elements which belong on this processor,
        // but by calling this on all processors a unique range in [0,_n_active_elem_on_proc[pid])
        // is constructed.  Only the indices for the elements we pass in are returned in the array.
        MeshCommunication().find_global_indices (mesh.comm(),
                                                 bbox, it, end,
                                                 global_index);

        for (dof_id_type cnt=0; it != end; ++it)
          {
            const Elem * elem = *it;
            // vectormap::count forces a sort, which is too expensive
            // in a loop
            // libmesh_assert (!_global_index_by_pid_map.count(elem->id()));
            libmesh_assert_less (cnt, global_index.size());
            libmesh_assert_less (global_index[cnt], _n_active_elem_on_proc[pid]);

            _global_index_by_pid_map.insert(std::make_pair(elem->id(), global_index[cnt++] + pid_offset));
          }

        pid_offset += _n_active_elem_on_proc[pid];
      }

    // create the unique mapping for all active elements independent of partitioning
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();

      // Calling this on all processors a unique range in [0,n_active_elem) is constructed.
      // Only the indices for the elements we pass in are returned in the array.
      MeshCommunication().find_global_indices (mesh.comm(),
                                               bbox, it, end,
                                               global_index);

      for (dof_id_type cnt=0; it != end; ++it)
        {
          const Elem * elem = *it;
          // vectormap::count forces a sort, which is too expensive
          // in a loop
          // libmesh_assert (!global_index_map.count(elem->id()));
          libmesh_assert_less (cnt, global_index.size());
          libmesh_assert_less (global_index[cnt], n_active_elem);

          global_index_map.insert(std::make_pair(elem->id(), global_index[cnt++]));
        }
    }
    // really, shouldn't be close!
    libmesh_assert_less_equal (global_index_map.size(), n_active_elem);
    libmesh_assert_less_equal (_global_index_by_pid_map.size(), n_active_elem);

    // At this point the two maps should be the same size.  If they are not
    // then the number of active elements is not the same as the sum over all
    // processors of the number of active elements per processor, which means
    // there must be some unpartitioned objects out there.
    if (global_index_map.size() != _global_index_by_pid_map.size())
      libmesh_error_msg("ERROR:  ParmetisPartitioner cannot handle unpartitioned objects!");
  }

  // Finally, we need to initialize the vertex (partition) weights and the initial subdomain
  // mapping.  The subdomain mapping will be independent of the processor mapping, and is
  // defined by a simple mapping of the global indices we just found.
  {
    std::vector<dof_id_type> subdomain_bounds(mesh.n_processors());

    const dof_id_type first_local_elem = _pmetis->vtxdist[mesh.processor_id()];

    for (processor_id_type pid=0; pid<mesh.n_processors(); pid++)
      {
        dof_id_type tgt_subdomain_size = 0;

        // watch out for the case that n_subdomains < n_processors
        if (pid < static_cast<unsigned int>(_pmetis->nparts))
          {
            tgt_subdomain_size = n_active_elem/std::min
              (cast_int<Parmetis::idx_t>(mesh.n_processors()), _pmetis->nparts);

            if (pid < n_active_elem%_pmetis->nparts)
              tgt_subdomain_size++;
          }
        if (pid == 0)
          subdomain_bounds[0] = tgt_subdomain_size;
        else
          subdomain_bounds[pid] = subdomain_bounds[pid-1] + tgt_subdomain_size;
      }

    libmesh_assert_equal_to (subdomain_bounds.back(), n_active_elem);

    MeshBase::const_element_iterator       elem_it  = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end();

    for (; elem_it != elem_end; ++elem_it)
      {
        const Elem * elem = *elem_it;

        libmesh_assert (_global_index_by_pid_map.count(elem->id()));
        const dof_id_type global_index_by_pid =
          _global_index_by_pid_map[elem->id()];
        libmesh_assert_less (global_index_by_pid, n_active_elem);

        const dof_id_type local_index =
          global_index_by_pid - first_local_elem;

        libmesh_assert_less (local_index, n_active_local_elem);
        libmesh_assert_less (local_index, _pmetis->vwgt.size());

        // TODO:[BSK] maybe there is a better weight?
        _pmetis->vwgt[local_index] = elem->n_nodes();

        // find the subdomain this element belongs in
        libmesh_assert (global_index_map.count(elem->id()));
        const dof_id_type global_index =
          global_index_map[elem->id()];

        libmesh_assert_less (global_index, subdomain_bounds.back());

        const unsigned int subdomain_id =
          std::distance(subdomain_bounds.begin(),
                        std::lower_bound(subdomain_bounds.begin(),
                                         subdomain_bounds.end(),
                                         global_index));
        libmesh_assert_less (subdomain_id, static_cast<unsigned int>(_pmetis->nparts));
        libmesh_assert_less (local_index, _pmetis->part.size());

        _pmetis->part[local_index] = subdomain_id;
      }
  }
}



void ParmetisPartitioner::build_graph (const MeshBase & mesh)
{
  LOG_SCOPE("build_graph()", "ParmetisPartitioner");

  // build the graph in distributed CSR format.  Note that
  // the edges in the graph will correspond to
  // face neighbors
  const dof_id_type n_active_local_elem  = mesh.n_active_local_elem();

  // If we have boundary elements in this mesh, we want to account for
  // the connectivity between them and interior elements.  We can find
  // interior elements from boundary elements, but we need to build up
  // a lookup map to do the reverse.

  typedef LIBMESH_BEST_UNORDERED_MULTIMAP<const Elem *, const Elem *>
    map_type;
  map_type interior_to_boundary_map;

  {
    MeshBase::const_element_iterator       elem_it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator elem_end = mesh.active_elements_end();

    for (; elem_it != elem_end; ++elem_it)
      {
        const Elem * elem = *elem_it;

        // If we don't have an interior_parent then there's nothing to look us
        // up.
        if ((elem->dim() >= LIBMESH_DIM) ||
            !elem->interior_parent())
          continue;

        // get all relevant interior elements
        std::set<const Elem *> neighbor_set;
        elem->find_interior_neighbors(neighbor_set);

        std::set<const Elem *>::iterator n_it = neighbor_set.begin();
        for (; n_it != neighbor_set.end(); ++n_it)
          {
            // FIXME - non-const versions of the Elem set methods
            // would be nice
            Elem * neighbor = const_cast<Elem *>(*n_it);

#if defined(LIBMESH_HAVE_UNORDERED_MULTIMAP) || \
  defined(LIBMESH_HAVE_TR1_UNORDERED_MAP) ||    \
  defined(LIBMESH_HAVE_HASH_MAP) ||             \
  defined(LIBMESH_HAVE_EXT_HASH_MAP)
            interior_to_boundary_map.insert
              (std::make_pair(neighbor, elem));
#else
            interior_to_boundary_map.insert
              (interior_to_boundary_map.begin(),
               std::make_pair(neighbor, elem));
#endif
          }
      }
  }

#ifdef LIBMESH_ENABLE_AMR
  std::vector<const Elem *> neighbors_offspring;
#endif

  std::vector<std::vector<dof_id_type> > graph(n_active_local_elem);
  dof_id_type graph_size=0;

  const dof_id_type first_local_elem = _pmetis->vtxdist[mesh.processor_id()];

  MeshBase::const_element_iterator       elem_it  = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      const Elem * elem = *elem_it;

      libmesh_assert (_global_index_by_pid_map.count(elem->id()));
      const dof_id_type global_index_by_pid =
        _global_index_by_pid_map[elem->id()];

      const dof_id_type local_index =
        global_index_by_pid - first_local_elem;
      libmesh_assert_less (local_index, n_active_local_elem);

      std::vector<dof_id_type> & graph_row = graph[local_index];

      // Loop over the element's neighbors.  An element
      // adjacency corresponds to a face neighbor
      for (unsigned int ms=0; ms<elem->n_neighbors(); ms++)
        {
          const Elem * neighbor = elem->neighbor_ptr(ms);

          if (neighbor != libmesh_nullptr)
            {
              // If the neighbor is active treat it
              // as a connection
              if (neighbor->active())
                {
                  libmesh_assert(_global_index_by_pid_map.count(neighbor->id()));
                  const dof_id_type neighbor_global_index_by_pid =
                    _global_index_by_pid_map[neighbor->id()];

                  graph_row.push_back(neighbor_global_index_by_pid);
                  graph_size++;
                }

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
                  // of the neighbor

                  // FIXME - this is the wrong thing, since we
                  // should be getting the active family tree on
                  // our side only.  But adding too many graph
                  // links may cause hanging nodes to tend to be
                  // on partition interiors, which would reduce
                  // communication overhead for constraint
                  // equations, so we'll leave it.

                  neighbor->active_family_tree (neighbors_offspring);

                  // Get all the neighbor's children that
                  // live on that side and are thus connected
                  // to us
                  for (std::size_t nc=0; nc<neighbors_offspring.size(); nc++)
                    {
                      const Elem * child =
                        neighbors_offspring[nc];

                      // This does not assume a level-1 mesh.
                      // Note that since children have sides numbered
                      // coincident with the parent then this is a sufficient test.
                      if (child->neighbor_ptr(ns) == elem)
                        {
                          libmesh_assert (child->active());
                          libmesh_assert (_global_index_by_pid_map.count(child->id()));
                          const dof_id_type child_global_index_by_pid =
                            _global_index_by_pid_map[child->id()];

                          graph_row.push_back(child_global_index_by_pid);
                          graph_size++;
                        }
                    }
                }

#endif /* ifdef LIBMESH_ENABLE_AMR */


            }
        }

      if ((elem->dim() < LIBMESH_DIM) &&
          elem->interior_parent())
        {
          // get all relevant interior elements
          std::set<const Elem *> neighbor_set;
          elem->find_interior_neighbors(neighbor_set);

          std::set<const Elem *>::iterator n_it = neighbor_set.begin();
          for (; n_it != neighbor_set.end(); ++n_it)
            {
              // FIXME - non-const versions of the Elem set methods
              // would be nice
              Elem * neighbor = const_cast<Elem *>(*n_it);

              const dof_id_type neighbor_global_index_by_pid =
                _global_index_by_pid_map[neighbor->id()];

              graph_row.push_back(neighbor_global_index_by_pid);
              graph_size++;
            }
        }

      // Check for any boundary neighbors
      typedef map_type::iterator map_it_type;
      std::pair<map_it_type, map_it_type>
        bounds = interior_to_boundary_map.equal_range(elem);

      for (map_it_type it = bounds.first; it != bounds.second; ++it)
        {
          const Elem * neighbor = it->second;

          const dof_id_type neighbor_global_index_by_pid =
            _global_index_by_pid_map[neighbor->id()];

          graph_row.push_back(neighbor_global_index_by_pid);
          graph_size++;
        }
    }

  // Reserve space in the adjacency array
  _pmetis->xadj.clear();
  _pmetis->xadj.reserve (n_active_local_elem + 1);
  _pmetis->adjncy.clear();
  _pmetis->adjncy.reserve (graph_size);

  for (std::size_t r=0; r<graph.size(); r++)
    {
      _pmetis->xadj.push_back(_pmetis->adjncy.size());
      std::vector<dof_id_type> graph_row; // build this emtpy
      graph_row.swap(graph[r]); // this will deallocate at the end of scope
      _pmetis->adjncy.insert(_pmetis->adjncy.end(),
                             graph_row.begin(),
                             graph_row.end());
    }

  // The end of the adjacency array for the last elem
  _pmetis->xadj.push_back(_pmetis->adjncy.size());

  libmesh_assert_equal_to (_pmetis->xadj.size(), n_active_local_elem+1);
  libmesh_assert_equal_to (_pmetis->adjncy.size(), graph_size);
}



void ParmetisPartitioner::assign_partitioning (MeshBase & mesh)
{
  LOG_SCOPE("assign_partitioning()", "ParmetisPartitioner");

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  const dof_id_type
    first_local_elem = _pmetis->vtxdist[mesh.processor_id()];

#ifndef NDEBUG
  const dof_id_type n_active_local_elem = mesh.n_active_local_elem();
#endif

  std::vector<std::vector<dof_id_type> >
    requested_ids(mesh.n_processors()),
    requests_to_fill(mesh.n_processors());

  MeshBase::element_iterator elem_it  = mesh.active_elements_begin();
  MeshBase::element_iterator elem_end = mesh.active_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem * elem = *elem_it;

      // we need to get the index from the owning processor
      // (note we cannot assign it now -- we are iterating
      // over elements again and this will be bad!)
      libmesh_assert_less (elem->processor_id(), requested_ids.size());
      requested_ids[elem->processor_id()].push_back(elem->id());
    }

  // Trade with all processors (including self) to get their indices
  for (processor_id_type pid=0; pid<mesh.n_processors(); pid++)
    {
      // Trade my requests with processor procup and procdown
      const processor_id_type procup = (mesh.processor_id() + pid) % mesh.n_processors();
      const processor_id_type procdown = (mesh.n_processors() +
                                          mesh.processor_id() - pid) % mesh.n_processors();

      mesh.comm().send_receive (procup,   requested_ids[procup],
                                procdown, requests_to_fill[procdown]);

      // we can overwrite these requested ids in-place.
      for (std::size_t i=0; i<requests_to_fill[procdown].size(); i++)
        {
          const dof_id_type requested_elem_index =
            requests_to_fill[procdown][i];

          libmesh_assert(_global_index_by_pid_map.count(requested_elem_index));

          const dof_id_type global_index_by_pid =
            _global_index_by_pid_map[requested_elem_index];

          const dof_id_type local_index =
            global_index_by_pid - first_local_elem;

          libmesh_assert_less (local_index, _pmetis->part.size());
          libmesh_assert_less (local_index, n_active_local_elem);

          const unsigned int elem_procid =
            static_cast<unsigned int>(_pmetis->part[local_index]);

          libmesh_assert_less (elem_procid, static_cast<unsigned int>(_pmetis->nparts));

          requests_to_fill[procdown][i] = elem_procid;
        }

      // Trade back
      mesh.comm().send_receive (procdown, requests_to_fill[procdown],
                                procup,   requested_ids[procup]);
    }

  // and finally assign the partitioning.
  // note we are iterating in exactly the same order
  // used to build up the request, so we can expect the
  // required entries to be in the proper sequence.
  elem_it  = mesh.active_elements_begin();
  elem_end = mesh.active_elements_end();

  for (std::vector<unsigned int> counters(mesh.n_processors(), 0);
       elem_it != elem_end; ++elem_it)
    {
      Elem * elem = *elem_it;

      const processor_id_type current_pid = elem->processor_id();

      libmesh_assert_less (counters[current_pid], requested_ids[current_pid].size());

      const processor_id_type elem_procid =
        requested_ids[current_pid][counters[current_pid]++];

      libmesh_assert_less (elem_procid, static_cast<unsigned int>(_pmetis->nparts));
      elem->processor_id() = elem_procid;
    }
}

#endif // #ifdef LIBMESH_HAVE_PARMETIS

} // namespace libMesh
