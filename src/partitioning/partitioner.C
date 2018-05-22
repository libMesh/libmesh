// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/partitioner.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/parallel_ghost_sync.h"

namespace libMesh
{



// ------------------------------------------------------------
// Partitioner static data
const dof_id_type Partitioner::communication_blocksize = 1000000;



// ------------------------------------------------------------
// Partitioner implementation
void Partitioner::partition (MeshBase & mesh)
{
  this->partition(mesh,mesh.n_processors());
}



void Partitioner::partition (MeshBase & mesh,
                             const unsigned int n)
{
  libmesh_parallel_only(mesh.comm());

  // BSK - temporary fix while redistribution is integrated 6/26/2008
  // Uncomment this to not repartition in parallel
  //   if (!mesh.is_serial())
  //     return;

  // we cannot partition into more pieces than we have
  // active elements!
  const unsigned int n_parts =
    static_cast<unsigned int>
    (std::min(mesh.n_active_elem(), static_cast<dof_id_type>(n)));

  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n_parts;

  if (n_parts == 1)
    {
      this->single_partition (mesh);
      return;
    }

  // First assign a temporary partitioning to any unpartitioned elements
  Partitioner::partition_unpartitioned_elements(mesh, n_parts);

  // Call the partitioning function
  this->_do_partition(mesh,n_parts);

  // Set the parent's processor ids
  Partitioner::set_parent_processor_ids(mesh);

  // Redistribute elements if necessary, before setting node processor
  // ids, to make sure those will be set consistently
  mesh.redistribute();

#ifdef DEBUG
  MeshTools::libmesh_assert_valid_remote_elems(mesh);

  // Messed up elem processor_id()s can leave us without the child
  // elements we need to restrict vectors on a distributed mesh
  MeshTools::libmesh_assert_valid_procids<Elem>(mesh);
#endif

  // Set the node's processor ids
  Partitioner::set_node_processor_ids(mesh);

#ifdef DEBUG
  MeshTools::libmesh_assert_valid_procids<Elem>(mesh);
#endif

  // Give derived Mesh classes a chance to update any cached data to
  // reflect the new partitioning
  mesh.update_post_partitioning();
}



void Partitioner::repartition (MeshBase & mesh)
{
  this->repartition(mesh,mesh.n_processors());
}



void Partitioner::repartition (MeshBase & mesh,
                               const unsigned int n)
{
  // we cannot partition into more pieces than we have
  // active elements!
  const unsigned int n_parts =
    static_cast<unsigned int>
    (std::min(mesh.n_active_elem(), static_cast<dof_id_type>(n)));

  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n_parts;

  if (n_parts == 1)
    {
      this->single_partition (mesh);
      return;
    }

  // First assign a temporary partitioning to any unpartitioned elements
  Partitioner::partition_unpartitioned_elements(mesh, n_parts);

  // Call the partitioning function
  this->_do_repartition(mesh,n_parts);

  // Set the parent's processor ids
  Partitioner::set_parent_processor_ids(mesh);

  // Set the node's processor ids
  Partitioner::set_node_processor_ids(mesh);
}





void Partitioner::single_partition (MeshBase & mesh)
{
  this->single_partition_range(mesh.elements_begin(),
                               mesh.elements_end());
}



void Partitioner::single_partition_range (MeshBase::element_iterator it,
                                          MeshBase::element_iterator end)
{
  LOG_SCOPE("single_partition_range()", "Partitioner");

  for (auto & elem : as_range(it, end))
    {
      elem->processor_id() = 0;

      // Assign all this element's nodes to processor 0 as well.
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        elem->node_ptr(n)->processor_id() = 0;
    }
}

void Partitioner::partition_unpartitioned_elements (MeshBase & mesh)
{
  Partitioner::partition_unpartitioned_elements(mesh, mesh.n_processors());
}



void Partitioner::partition_unpartitioned_elements (MeshBase & mesh,
                                                    const unsigned int n_subdomains)
{
  MeshBase::element_iterator       it  = mesh.unpartitioned_elements_begin();
  const MeshBase::element_iterator end = mesh.unpartitioned_elements_end();

  const dof_id_type n_unpartitioned_elements = MeshTools::n_elem (it, end);

  // the unpartitioned elements must exist on all processors. If the range is empty on one
  // it is empty on all, and we can quit right here.
  if (!n_unpartitioned_elements) return;

  // find the target subdomain sizes
  std::vector<dof_id_type> subdomain_bounds(mesh.n_processors());

  for (processor_id_type pid=0; pid<mesh.n_processors(); pid++)
    {
      dof_id_type tgt_subdomain_size = 0;

      // watch out for the case that n_subdomains < n_processors
      if (pid < n_subdomains)
        {
          tgt_subdomain_size = n_unpartitioned_elements/n_subdomains;

          if (pid < n_unpartitioned_elements%n_subdomains)
            tgt_subdomain_size++;

        }

      //libMesh::out << "pid, #= " << pid << ", " << tgt_subdomain_size << std::endl;
      if (pid == 0)
        subdomain_bounds[0] = tgt_subdomain_size;
      else
        subdomain_bounds[pid] = subdomain_bounds[pid-1] + tgt_subdomain_size;
    }

  libmesh_assert_equal_to (subdomain_bounds.back(), n_unpartitioned_elements);

  // create the unique mapping for all unpartitioned elements independent of partitioning
  // determine the global indexing for all the unpartitioned elements
  std::vector<dof_id_type> global_indices;

  // Calling this on all processors a unique range in [0,n_unpartitioned_elements) is constructed.
  // Only the indices for the elements we pass in are returned in the array.
  MeshCommunication().find_global_indices (mesh.comm(),
                                           MeshTools::create_bounding_box(mesh), it, end,
                                           global_indices);

  dof_id_type cnt=0;
  for (auto & elem : as_range(it, end))
    {
      libmesh_assert_less (cnt, global_indices.size());
      const dof_id_type global_index =
        global_indices[cnt++];

      libmesh_assert_less (global_index, subdomain_bounds.back());
      libmesh_assert_less (global_index, n_unpartitioned_elements);

      const processor_id_type subdomain_id =
        cast_int<processor_id_type>
        (std::distance(subdomain_bounds.begin(),
                       std::upper_bound(subdomain_bounds.begin(),
                                        subdomain_bounds.end(),
                                        global_index)));
      libmesh_assert_less (subdomain_id, n_subdomains);

      elem->processor_id() = subdomain_id;
      //libMesh::out << "assigning " << global_index << " to " << subdomain_id << std::endl;
    }
}



void Partitioner::set_parent_processor_ids(MeshBase & mesh)
{
  // Ignore the parameter when !LIBMESH_ENABLE_AMR
  libmesh_ignore(mesh);

  LOG_SCOPE("set_parent_processor_ids()", "Partitioner");

#ifdef LIBMESH_ENABLE_AMR

  // If the mesh is serial we have access to all the elements,
  // in particular all the active ones.  We can therefore set
  // the parent processor ids indirectly through their children, and
  // set the subactive processor ids while examining their active
  // ancestors.
  // By convention a parent is assigned to the minimum processor
  // of all its children, and a subactive is assigned to the processor
  // of its active ancestor.
  if (mesh.is_serial())
    {
      for (auto & child : mesh.active_element_ptr_range())
        {
          // First set descendents
          std::vector<const Elem *> subactive_family;
          child->total_family_tree(subactive_family);
          for (std::size_t i = 0; i != subactive_family.size(); ++i)
            const_cast<Elem *>(subactive_family[i])->processor_id() = child->processor_id();

          // Then set ancestors
          Elem * parent = child->parent();

          while (parent)
            {
              // invalidate the parent id, otherwise the min below
              // will not work if the current parent id is less
              // than all the children!
              parent->invalidate_processor_id();

              for (auto & child : parent->child_ref_range())
                {
                  libmesh_assert(!child.is_remote());
                  libmesh_assert_not_equal_to (child.processor_id(), DofObject::invalid_processor_id);
                  parent->processor_id() = std::min(parent->processor_id(),
                                                    child.processor_id());
                }
              parent = parent->parent();
            }
        }
    }

  // When the mesh is parallel we cannot guarantee that parents have access to
  // all their children.
  else
    {
      // Setting subactive processor ids is easy: we can guarantee
      // that children have access to all their parents.

      // Loop over all the active elements in the mesh
      for (auto & child : mesh.active_element_ptr_range())
        {
          std::vector<const Elem *> subactive_family;
          child->total_family_tree(subactive_family);
          for (std::size_t i = 0; i != subactive_family.size(); ++i)
            const_cast<Elem *>(subactive_family[i])->processor_id() = child->processor_id();
        }

      // When the mesh is parallel we cannot guarantee that parents have access to
      // all their children.

      // We will use a brute-force approach here.  Each processor finds its parent
      // elements and sets the parent pid to the minimum of its
      // semilocal descendants.
      // A global reduction is then performed to make sure the true minimum is found.
      // As noted, this is required because we cannot guarantee that a parent has
      // access to all its children on any single processor.
      libmesh_parallel_only(mesh.comm());
      libmesh_assert(MeshTools::n_elem(mesh.unpartitioned_elements_begin(),
                                       mesh.unpartitioned_elements_end()) == 0);

      const dof_id_type max_elem_id = mesh.max_elem_id();

      std::vector<processor_id_type>
        parent_processor_ids (std::min(communication_blocksize,
                                       max_elem_id));

      for (dof_id_type blk=0, last_elem_id=0; last_elem_id<max_elem_id; blk++)
        {
          last_elem_id =
            std::min(static_cast<dof_id_type>((blk+1)*communication_blocksize),
                     max_elem_id);
          const dof_id_type first_elem_id = blk*communication_blocksize;

          std::fill (parent_processor_ids.begin(),
                     parent_processor_ids.end(),
                     DofObject::invalid_processor_id);

          // first build up local contributions to parent_processor_ids
          bool have_parent_in_block = false;

          for (auto & parent : as_range(mesh.ancestor_elements_begin(),
                                        mesh.ancestor_elements_end()))
            {
              const dof_id_type parent_idx = parent->id();
              libmesh_assert_less (parent_idx, max_elem_id);

              if ((parent_idx >= first_elem_id) &&
                  (parent_idx <  last_elem_id))
                {
                  have_parent_in_block = true;
                  processor_id_type parent_pid = DofObject::invalid_processor_id;

                  std::vector<const Elem *> active_family;
                  parent->active_family_tree(active_family);
                  for (std::size_t i = 0; i != active_family.size(); ++i)
                    parent_pid = std::min (parent_pid, active_family[i]->processor_id());

                  const dof_id_type packed_idx = parent_idx - first_elem_id;
                  libmesh_assert_less (packed_idx, parent_processor_ids.size());

                  parent_processor_ids[packed_idx] = parent_pid;
                }
            }

          // then find the global minimum
          mesh.comm().min (parent_processor_ids);

          // and assign the ids, if we have a parent in this block.
          if (have_parent_in_block)
            for (auto & parent : as_range(mesh.ancestor_elements_begin(),
                                          mesh.ancestor_elements_end()))
              {
                const dof_id_type parent_idx = parent->id();

                if ((parent_idx >= first_elem_id) &&
                    (parent_idx <  last_elem_id))
                  {
                    const dof_id_type packed_idx = parent_idx - first_elem_id;
                    libmesh_assert_less (packed_idx, parent_processor_ids.size());

                    const processor_id_type parent_pid =
                      parent_processor_ids[packed_idx];

                    libmesh_assert_not_equal_to (parent_pid, DofObject::invalid_processor_id);

                    parent->processor_id() = parent_pid;
                  }
              }
        }
    }

#endif // LIBMESH_ENABLE_AMR
}



void Partitioner::set_node_processor_ids(MeshBase & mesh)
{
  LOG_SCOPE("set_node_processor_ids()","Partitioner");

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  // If we have any unpartitioned elements at this
  // stage there is a problem
  libmesh_assert (MeshTools::n_elem(mesh.unpartitioned_elements_begin(),
                                    mesh.unpartitioned_elements_end()) == 0);


  //   const dof_id_type orig_n_local_nodes = mesh.n_local_nodes();

  //   libMesh::err << "[" << mesh.processor_id() << "]: orig_n_local_nodes="
  //     << orig_n_local_nodes << std::endl;

  // Build up request sets.  Each node is currently owned by a processor because
  // it is connected to an element owned by that processor.  However, during the
  // repartitioning phase that element may have been assigned a new processor id, but
  // it is still resident on the original processor.  We need to know where to look
  // for new ids before assigning new ids, otherwise we may be asking the wrong processors
  // for the wrong information.
  //
  // The only remaining issue is what to do with unpartitioned nodes.  Since they are required
  // to live on all processors we can simply rely on ourselves to number them properly.
  std::map<processor_id_type, std::vector<dof_id_type>>
    requested_node_ids;

  // Loop over all the nodes, count the ones on each processor.  We can skip ourself
  std::vector<dof_id_type> ghost_nodes_from_proc(mesh.n_processors(), 0);

  for (auto & node : mesh.node_ptr_range())
    {
      libmesh_assert(node);
      const processor_id_type current_pid = node->processor_id();
      if (current_pid != mesh.processor_id() &&
          current_pid != DofObject::invalid_processor_id)
        {
          libmesh_assert_less (current_pid, ghost_nodes_from_proc.size());
          ghost_nodes_from_proc[current_pid]++;
        }
    }

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (processor_id_type pid=0; pid != mesh.n_processors(); ++pid)
    if (ghost_nodes_from_proc[pid])
      requested_node_ids[pid].reserve(ghost_nodes_from_proc[pid]);

  // We need to get the new pid for each node from the processor
  // which *currently* owns the node.  We can safely skip ourself
  for (auto & node : mesh.node_ptr_range())
    {
      libmesh_assert(node);
      const processor_id_type current_pid = node->processor_id();
      if (current_pid != mesh.processor_id() &&
          current_pid != DofObject::invalid_processor_id)
        {
          libmesh_assert_less (requested_node_ids[current_pid].size(),
                               ghost_nodes_from_proc[current_pid]);
          requested_node_ids[current_pid].push_back(node->id());
        }

      // Unset any previously-set node processor ids
      node->invalidate_processor_id();
    }

  // Loop over all the active elements
  for (auto & elem : mesh.active_element_ptr_range())
    {
      libmesh_assert(elem);

      libmesh_assert_not_equal_to (elem->processor_id(), DofObject::invalid_processor_id);

      // Consider updating the processor id on this element's nodes
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        {
          Node & node = elem->node_ref(n);
          processor_id_type & pid = node.processor_id();
          pid = node.choose_processor_id(pid, elem->processor_id());
        }
    }

  // And loop over the subactive elements, but don't reassign
  // nodes that are already active on another processor.
  for (auto & elem : as_range(mesh.subactive_elements_begin(),
                              mesh.subactive_elements_end()))
    {
      libmesh_assert(elem);

      libmesh_assert_not_equal_to (elem->processor_id(), DofObject::invalid_processor_id);

      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        if (elem->node_ptr(n)->processor_id() == DofObject::invalid_processor_id)
          elem->node_ptr(n)->processor_id() = elem->processor_id();
    }

  // Same for the inactive elements -- we will have already gotten most of these
  // nodes, *except* for the case of a parent with a subset of children which are
  // ghost elements.  In that case some of the parent nodes will not have been
  // properly handled yet
  for (auto & elem : as_range(mesh.not_active_elements_begin(),
                              mesh.not_active_elements_end()))
    {
      libmesh_assert(elem);

      libmesh_assert_not_equal_to (elem->processor_id(), DofObject::invalid_processor_id);

      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        if (elem->node_ptr(n)->processor_id() == DofObject::invalid_processor_id)
          elem->node_ptr(n)->processor_id() = elem->processor_id();
    }

  // We can't assert that all nodes are connected to elements, because
  // a DistributedMesh with NodeConstraints might have pulled in some
  // remote nodes solely for evaluating those constraints.
  // MeshTools::libmesh_assert_connected_nodes(mesh);

  // For such nodes, we'll do a sanity check later when making sure
  // that we successfully reset their processor ids to something
  // valid.

  auto gather_functor =
    [& mesh]
    (processor_id_type, const std::vector<dof_id_type> & ids,
     std::vector<processor_id_type> & new_pids)
    {
      const std::size_t ids_size = ids.size();
      new_pids.resize(ids_size);

      // Fill those requests in-place
      for (std::size_t i=0; i != ids_size; ++i)
        {
          Node & node = mesh.node_ref(ids[i]);
          const processor_id_type new_pid = node.processor_id();

          // We may have an invalid processor_id() on nodes that have been
          // "detached" from coarsened-away elements but that have not yet
          // themselves been removed.
          // libmesh_assert_not_equal_to (new_pid, DofObject::invalid_processor_id);
          // libmesh_assert_less (new_pid, mesh.n_partitions()); // this is the correct test --
          new_pids[i] = new_pid;                                 //  the number of partitions may
        }                                                        //  not equal the number of processors
    };

  auto action_functor =
    [& mesh]
    (processor_id_type,
     const std::vector<dof_id_type> & ids,
     const std::vector<processor_id_type> & new_pids)
    {
      const std::size_t ids_size = ids.size();
      // Copy the pid changes we've now been informed of
      for (std::size_t i=0; i != ids_size; ++i)
        {
          Node & node = mesh.node_ref(ids[i]);

          // this is the correct test -- the number of partitions may
          // not equal the number of processors

          // But: we may have an invalid processor_id() on nodes that
          // have been "detached" from coarsened-away elements but
          // that have not yet themselves been removed.
          // libmesh_assert_less (filled_request[i], mesh.n_partitions());

          node.processor_id(new_pids[i]);
        }
    };

  const processor_id_type * ex = libmesh_nullptr;
  Parallel::pull_parallel_vector_data
    (mesh.comm(), requested_node_ids, gather_functor, action_functor, ex);

#ifdef DEBUG
  MeshTools::libmesh_assert_valid_procids<Node>(mesh);
  MeshTools::libmesh_assert_canonical_node_procids(mesh);
#endif
}

struct SyncLocalIDs
{
  typedef dof_id_type datum;

  typedef std::unordered_map<dof_id_type, dof_id_type> map_type;

  SyncLocalIDs(map_type & _id_map) : id_map(_id_map) {}

  map_type & id_map;

  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & local_ids) const
  {
    local_ids.resize(ids.size());

    for (std::size_t i=0, imax = ids.size(); i != imax; ++i)
      local_ids[i] = id_map[ids[i]];
  }

  void act_on_data (const std::vector<dof_id_type> & ids,
                    const std::vector<datum> & local_ids)
  {
    for (std::size_t i=0, imax = local_ids.size(); i != imax; ++i)
      id_map[ids[i]] = local_ids[i];
  }
};

void Partitioner::_find_global_index_by_pid_map(const MeshBase & mesh)
{
  const dof_id_type n_active_local_elem = mesh.n_active_local_elem();

  // Find the number of active elements on each processor.  We cannot use
  // mesh.n_active_elem_on_proc(pid) since that only returns the number of
  // elements assigned to pid which are currently stored on the calling
  // processor. This will not in general be correct for parallel meshes
  // when (pid!=mesh.processor_id()).
  _n_active_elem_on_proc.resize(mesh.n_processors());
  mesh.comm().allgather(n_active_local_elem, _n_active_elem_on_proc);

  libMesh::BoundingBox bbox =
    MeshTools::create_bounding_box(mesh);

  _global_index_by_pid_map.clear();

  // create the mapping which is contiguous by processor
  {
    MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_local_elements_end();

    MeshCommunication().find_local_indices (bbox, it, end,
                                              _global_index_by_pid_map);
  }

  SyncLocalIDs sync(_global_index_by_pid_map);

  Parallel::sync_dofobject_data_by_id
      (mesh.comm(), mesh.active_elements_begin(), mesh.active_elements_end(), sync);

  dof_id_type pid_offset=0;
  for (processor_id_type pid=0; pid<mesh.n_processors(); pid++)
    {
      MeshBase::const_element_iterator       it  = mesh.active_pid_elements_begin(pid);
      const MeshBase::const_element_iterator end = mesh.active_pid_elements_end(pid);

      for (; it != end; ++it)
        {
          const Elem * elem = *it;
          libmesh_assert_less (_global_index_by_pid_map[elem->id()], _n_active_elem_on_proc[pid]);

          _global_index_by_pid_map[elem->id()] += pid_offset;
        }

      pid_offset += _n_active_elem_on_proc[pid];
    }
}

void Partitioner::build_graph (const MeshBase & mesh)
{
  LOG_SCOPE("build_graph()", "ParmetisPartitioner");

  const dof_id_type n_active_local_elem  = mesh.n_active_local_elem();
  // If we have boundary elements in this mesh, we want to account for
  // the connectivity between them and interior elements.  We can find
  // interior elements from boundary elements, but we need to build up
  // a lookup map to do the reverse.
  typedef std::unordered_multimap<const Elem *, const Elem *> map_type;
  map_type interior_to_boundary_map;

  for (const auto & elem : mesh.active_element_ptr_range())
    {
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
        interior_to_boundary_map.insert(std::make_pair(*n_it, elem));
    }

#ifdef LIBMESH_ENABLE_AMR
  std::vector<const Elem *> neighbors_offspring;
#endif

   if (!_global_index_by_pid_map.size())
     _find_global_index_by_pid_map(mesh);

   dof_id_type first_local_elem = 0;
   for (processor_id_type pid=0; pid < mesh.processor_id(); pid++)
     first_local_elem += _n_active_elem_on_proc[pid];

  _dual_graph.clear();
  _dual_graph.resize(n_active_local_elem);

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      libmesh_assert (_global_index_by_pid_map.count(elem->id()));
      const dof_id_type global_index_by_pid =
        _global_index_by_pid_map[elem->id()];

      const dof_id_type local_index =
        global_index_by_pid - first_local_elem;
      libmesh_assert_less (local_index, n_active_local_elem);

      std::vector<dof_id_type> & graph_row = _dual_graph[local_index];

      // Loop over the element's neighbors.  An element
      // adjacency corresponds to a face neighbor
      for (auto neighbor : elem->neighbor_ptr_range())
        {
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
            }
        }

      // Check for any boundary neighbors
      for (const auto & pr : as_range(interior_to_boundary_map.equal_range(elem)))
        {
          const Elem * neighbor = pr.second;

          const dof_id_type neighbor_global_index_by_pid =
            _global_index_by_pid_map[neighbor->id()];

          graph_row.push_back(neighbor_global_index_by_pid);
        }
    }

}

void Partitioner::assign_partitioning (const MeshBase & mesh, const std::vector<dof_id_type> & parts)
{
  LOG_SCOPE("assign_partitioning()", "ParmetisPartitioner");

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  dof_id_type first_local_elem = 0;
  for (processor_id_type pid=0; pid < mesh.processor_id(); pid++)
    first_local_elem += _n_active_elem_on_proc[pid];

#ifndef NDEBUG
  const dof_id_type n_active_local_elem = mesh.n_active_local_elem();
#endif

  std::map<processor_id_type, std::vector<dof_id_type>>
    requested_ids;

  // Results to gather from each processor - kept in a map so we
  // do only one loop over elements after all receives are done.
  std::map<processor_id_type, std::vector<processor_id_type>>
    filled_request;

  for (auto & elem : mesh.active_element_ptr_range())
    {
      // we need to get the index from the owning processor
      // (note we cannot assign it now -- we are iterating
      // over elements again and this will be bad!)
      requested_ids[elem->processor_id()].push_back(elem->id());
    }

  auto gather_functor =
    [this,
     & parts,
#ifndef NDEBUG
     & mesh,
     n_active_local_elem,
#endif
     first_local_elem]
    (processor_id_type, const std::vector<dof_id_type> & ids,
     std::vector<processor_id_type> & data)
    {
      const std::size_t ids_size = ids.size();
      data.resize(ids.size());

      for (std::size_t i=0; i != ids_size; i++)
        {
          const dof_id_type requested_elem_index = ids[i];

          libmesh_assert(_global_index_by_pid_map.count(requested_elem_index));

          const dof_id_type global_index_by_pid =
            _global_index_by_pid_map[requested_elem_index];

          const dof_id_type local_index =
            global_index_by_pid - first_local_elem;

          libmesh_assert_less (local_index, parts.size());
          libmesh_assert_less (local_index, n_active_local_elem);

          const processor_id_type elem_procid =
            cast_int<processor_id_type>(parts[local_index]);

          libmesh_assert_less (elem_procid, mesh.n_partitions());

          data[i] = elem_procid;
        }
    };

  auto action_functor =
    [&filled_request]
    (processor_id_type pid,
     const std::vector<dof_id_type> &,
     const std::vector<processor_id_type> & new_procids)
    {
      filled_request[pid] = new_procids;
    };

  // Trade requests with other processors
  const processor_id_type * ex = libmesh_nullptr;
  Parallel::pull_parallel_vector_data
    (mesh.comm(), requested_ids, gather_functor, action_functor, ex);

  // and finally assign the partitioning.
  // note we are iterating in exactly the same order
  // used to build up the request, so we can expect the
  // required entries to be in the proper sequence.
  std::vector<unsigned int> counters(mesh.n_processors(), 0);
  for (auto & elem : mesh.active_element_ptr_range())
    {
      const processor_id_type current_pid = elem->processor_id();

      libmesh_assert_less (counters[current_pid], requested_ids[current_pid].size());

      const processor_id_type elem_procid =
        filled_request[current_pid][counters[current_pid]++];

      libmesh_assert_less (elem_procid, mesh.n_partitions());
      elem->processor_id() = elem_procid;
    }
}


} // namespace libMesh
