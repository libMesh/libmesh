// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/parmetis_partitioner.h"

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/metis_partitioner.h"
#include "libmesh/parallel_only.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/parmetis_helper.h"

// TIMPI includes
#include "timpi/communicator.h"    // also includes mpi.h
#include "timpi/parallel_implementation.h"    // for min()

// Include the ParMETIS header file.
#ifdef LIBMESH_HAVE_PARMETIS

// Before we include a header wrapped in a namespace, we'd better make
// sure none of its dependencies end up in that namespace
#include <mpi.h>

namespace Parmetis {
extern "C" {
#     include "libmesh/ignore_warnings.h"
#     include "parmetis.h"
#     include "libmesh/restore_warnings.h"
}
}

#endif


// C++ includes
#include <unordered_map>


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
  :  _pmetis(libmesh_make_unique<ParmetisHelper>())
#endif
{}



ParmetisPartitioner::ParmetisPartitioner (const ParmetisPartitioner & other)
  : Partitioner(other)
#ifdef LIBMESH_HAVE_PARMETIS
  , _pmetis(libmesh_make_unique<ParmetisHelper>(*(other._pmetis)))
#endif
{
}



ParmetisPartitioner::~ParmetisPartitioner() = default;



void ParmetisPartitioner::_do_partition (MeshBase & mesh,
                                         const unsigned int n_sbdmns)
{
  this->_do_repartition (mesh, n_sbdmns);
}



void ParmetisPartitioner::_do_repartition (MeshBase & mesh,
                                           const unsigned int n_sbdmns)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  // Check for easy returns
  if (!mesh.n_elem())
    return;

  if (n_sbdmns == 1)
    {
      this->single_partition(mesh);
      return;
    }

  libmesh_assert_greater (n_sbdmns, 0);

  // What to do if the Parmetis library IS NOT present
#ifndef LIBMESH_HAVE_PARMETIS

  libmesh_do_once(
  libMesh::out << "ERROR: The library has been built without" << std::endl
               << "Parmetis support.  Using a Metis"          << std::endl
               << "partitioner instead!"                      << std::endl;);

  MetisPartitioner mp;

  // Metis and other fallbacks only work in serial, and need to get
  // handed element ranges from an already-serialized mesh.
  mesh.allgather();

  // Don't just call partition() here; that would end up calling
  // post-element-partitioning work redundantly (and at the moment
  // incorrectly)
  mp.partition_range (mesh, mesh.active_elements_begin(),
                      mesh.active_elements_end(), n_sbdmns);

  // What to do if the Parmetis library IS present
#else

  // Revert to METIS on one processor.
  if (mesh.n_processors() == 1)
    {
      // Make sure the mesh knows it's serial
      mesh.allgather();

      MetisPartitioner mp;
      // Don't just call partition() here; that would end up calling
      // post-element-partitioning work redundantly (and at the moment
      // incorrectly)
      mp.partition_range (mesh, mesh.active_elements_begin(),
                          mesh.active_elements_end(), n_sbdmns);
      return;
    }

  LOG_SCOPE("repartition()", "ParmetisPartitioner");

  // Initialize the data structures required by ParMETIS
  this->initialize (mesh, n_sbdmns);

  // build the graph corresponding to the mesh
  this->build_graph (mesh);

  // Make sure all processors have enough active local elements and
  // enough connectivity among them.
  // Parmetis tends to die when it's given only a couple elements
  // per partition or when it can't reach elements from each other.
  {
    bool ready_for_parmetis = true;
    for (const auto & nelem : _n_active_elem_on_proc)
      if (nelem < MIN_ELEM_PER_PROC)
        ready_for_parmetis = false;

    std::size_t my_adjacency = _pmetis->adjncy.size();
    mesh.comm().min(my_adjacency);
    if (!my_adjacency)
      ready_for_parmetis = false;

    // Parmetis will not work unless each processor has some
    // elements. Specifically, it will abort when passed a nullptr
    // partition or adjacency array on *any* of the processors.
    if (!ready_for_parmetis)
      {
        // FIXME: revert to METIS, although this requires a serial mesh
        MeshSerializer serialize(mesh);
        MetisPartitioner mp;
        mp.partition (mesh, n_sbdmns);
        return;
      }
  }


  // Partition the graph
  std::vector<Parmetis::idx_t> vsize(_pmetis->vwgt.size(), 1);
  Parmetis::real_t itr = 1000000.0;
  MPI_Comm mpi_comm = mesh.comm().get();

  // Call the ParMETIS adaptive repartitioning method.  This respects the
  // original partitioning when computing the new partitioning so as to
  // minimize the required data redistribution.
  Parmetis::ParMETIS_V3_AdaptiveRepart(_pmetis->vtxdist.empty() ? nullptr : _pmetis->vtxdist.data(),
                                       _pmetis->xadj.empty()    ? nullptr : _pmetis->xadj.data(),
                                       _pmetis->adjncy.empty()  ? nullptr : _pmetis->adjncy.data(),
                                       _pmetis->vwgt.empty()    ? nullptr : _pmetis->vwgt.data(),
                                       vsize.empty()            ? nullptr : vsize.data(),
                                       nullptr,
                                       &_pmetis->wgtflag,
                                       &_pmetis->numflag,
                                       &_pmetis->ncon,
                                       &_pmetis->nparts,
                                       _pmetis->tpwgts.empty()  ? nullptr : _pmetis->tpwgts.data(),
                                       _pmetis->ubvec.empty()   ? nullptr : _pmetis->ubvec.data(),
                                       &itr,
                                       _pmetis->options.data(),
                                       &_pmetis->edgecut,
                                       _pmetis->part.empty()    ? nullptr : reinterpret_cast<Parmetis::idx_t *>(_pmetis->part.data()),
                                       &mpi_comm);

  // Assign the returned processor ids
  this->assign_partitioning (mesh, _pmetis->part);

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

  // ParMetis expects the elements to be numbered in contiguous blocks
  // by processor, i.e. [0, ne0), [ne0, ne0+ne1), ...
  // Since we only partition active elements we should have no expectation
  // that we currently have such a distribution.  So we need to create it.
  // Also, at the same time we are going to map all the active elements into a globally
  // unique range [0,n_active_elem) which is *independent* of the current partitioning.
  // This can be fed to ParMetis as the initial partitioning of the subdomains (decoupled
  // from the partitioning of the objects themselves).  This allows us to get the same
  // resultant partitioning independent of the input partitioning.
  libMesh::BoundingBox bbox =
    MeshTools::create_bounding_box(mesh);

  _find_global_index_by_pid_map(mesh);


  // count the total number of active elements in the mesh.  Note we cannot
  // use mesh.n_active_elem() in general since this only returns the number
  // of active elements which are stored on the calling processor.
  // We should not use n_active_elem for any allocation because that will
  // be inherently unscalable, but it can be useful for libmesh_assertions.
  dof_id_type n_active_elem=0;

  // Set up the vtxdist array.  This will be the same on each processor.
  // ***** Consult the Parmetis documentation. *****
  libmesh_assert_equal_to (_pmetis->vtxdist.size(),
                           cast_int<std::size_t>(mesh.n_processors()+1));
  libmesh_assert_equal_to (_pmetis->vtxdist[0], 0);

  for (auto pid : make_range(mesh.n_processors()))
    {
      _pmetis->vtxdist[pid+1] = _pmetis->vtxdist[pid] + _n_active_elem_on_proc[pid];
      n_active_elem += _n_active_elem_on_proc[pid];
    }
  libmesh_assert_equal_to (_pmetis->vtxdist.back(), static_cast<Parmetis::idx_t>(n_active_elem));


  // Maps active element ids into a contiguous range independent of partitioning.
  // (only needs local scope)
  std::unordered_map<dof_id_type, dof_id_type> global_index_map;

  {
    std::vector<dof_id_type> global_index;

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

          global_index_map.emplace(elem->id(), global_index[cnt++]);
        }
    }
    // really, shouldn't be close!
    libmesh_assert_less_equal (global_index_map.size(), n_active_elem);
    libmesh_assert_less_equal (_global_index_by_pid_map.size(), n_active_elem);

    // At this point the two maps should be the same size.  If they are not
    // then the number of active elements is not the same as the sum over all
    // processors of the number of active elements per processor, which means
    // there must be some unpartitioned objects out there.
    libmesh_error_msg_if(global_index_map.size() != _global_index_by_pid_map.size(),
                         "ERROR:  ParmetisPartitioner cannot handle unpartitioned objects!");
  }

  // Finally, we need to initialize the vertex (partition) weights and the initial subdomain
  // mapping.  The subdomain mapping will be independent of the processor mapping, and is
  // defined by a simple mapping of the global indices we just found.
  {
    std::vector<dof_id_type> subdomain_bounds(mesh.n_processors());

    const dof_id_type first_local_elem = _pmetis->vtxdist[mesh.processor_id()];

    for (auto pid : make_range(mesh.n_processors()))
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

    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
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
          cast_int<unsigned int>
          (std::distance(subdomain_bounds.begin(),
                         std::lower_bound(subdomain_bounds.begin(),
                                          subdomain_bounds.end(),
                                          global_index)));
        libmesh_assert_less (subdomain_id, _pmetis->nparts);
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

  Partitioner::build_graph(mesh);

  dof_id_type graph_size=0;

  for (auto & row: _dual_graph)
   graph_size += cast_int<dof_id_type>(row.size());

  // Reserve space in the adjacency array
  _pmetis->xadj.clear();
  _pmetis->xadj.reserve (n_active_local_elem + 1);
  _pmetis->adjncy.clear();
  _pmetis->adjncy.reserve (graph_size);

  for (auto & graph_row : _dual_graph)
    {
      _pmetis->xadj.push_back(cast_int<int>(_pmetis->adjncy.size()));
      _pmetis->adjncy.insert(_pmetis->adjncy.end(),
                             graph_row.begin(),
                             graph_row.end());
    }

  // The end of the adjacency array for the last elem
  _pmetis->xadj.push_back(cast_int<int>(_pmetis->adjncy.size()));

  libmesh_assert_equal_to (_pmetis->xadj.size(), n_active_local_elem+1);
  libmesh_assert_equal_to (_pmetis->adjncy.size(), graph_size);
}

#endif // #ifdef LIBMESH_HAVE_PARMETIS

} // namespace libMesh
