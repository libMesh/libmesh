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
#include "libmesh/partitioner.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/parallel_ghost_sync.h"

// TIMPI includes
#include "timpi/parallel_implementation.h"
#include "timpi/parallel_sync.h"

// C/C++ includes
#ifdef LIBMESH_HAVE_PETSC
#include "libmesh/ignore_warnings.h"
#include "petscmat.h"
#include "libmesh/restore_warnings.h"
#endif


namespace {

using namespace libMesh;

struct CorrectProcIds
{
  typedef processor_id_type datum;

  CorrectProcIds(MeshBase & _mesh,
                 std::unordered_set<dof_id_type> & _bad_pids) :
    mesh(_mesh), bad_pids(_bad_pids) {}

  MeshBase & mesh;
  std::unordered_set<dof_id_type> & bad_pids;

  // ------------------------------------------------------------
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & data)
  {
    // Find the processor id of each requested node
    data.resize(ids.size());

    for (auto i : index_range(ids))
      {
        // Look for this point in the mesh and make sure we have a
        // good pid for it before sending that on
        if (ids[i] != DofObject::invalid_id &&
            !bad_pids.count(ids[i]))
          {
            Node & node = mesh.node_ref(ids[i]);

            // Return the node's correct processor id,
            data[i] = node.processor_id();
          }
        else
          data[i] = DofObject::invalid_processor_id;
      }
  }

  // ------------------------------------------------------------
  bool act_on_data (const std::vector<dof_id_type> & ids,
                    const std::vector<datum> proc_ids)
  {
    bool data_changed = false;

    // Set the ghost node processor ids we've now been informed of
    for (auto i : index_range(ids))
      {
        Node & node = mesh.node_ref(ids[i]);

        processor_id_type & proc_id = node.processor_id();
        const processor_id_type new_proc_id = proc_ids[i];

        auto it = bad_pids.find(ids[i]);

        if ((new_proc_id != proc_id ||
             it != bad_pids.end()) &&
            new_proc_id != DofObject::invalid_processor_id)
          {
            if (it != bad_pids.end())
              {
                proc_id = new_proc_id;
                bad_pids.erase(it);
              }
            else
              proc_id = node.choose_processor_id(proc_id, new_proc_id);

            if (proc_id == new_proc_id)
              data_changed = true;
          }
      }

    return data_changed;
  }
};

}

namespace libMesh
{



// ------------------------------------------------------------
// Partitioner static data
const dof_id_type Partitioner::communication_blocksize =
  dof_id_type(1000000);



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

  // Redistribute, in case someone (like our unit tests) is doing
  // something silly (like moving a whole already-distributed mesh
  // back onto rank 0).
  mesh.redistribute();
}



void Partitioner::single_partition_range (MeshBase::element_iterator it,
                                          MeshBase::element_iterator end)
{
  LOG_SCOPE("single_partition_range()", "Partitioner");

  for (auto & elem : as_range(it, end))
    {
      elem->processor_id() = 0;

      // Assign all this element's nodes to processor 0 as well.
      for (Node & node : elem->node_ref_range())
        node.processor_id() = 0;
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
  if (!n_unpartitioned_elements)
    return;

  // find the target subdomain sizes
  std::vector<dof_id_type> subdomain_bounds(mesh.n_processors());

  for (auto pid : make_range(mesh.n_processors()))
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
      for (auto & elem : mesh.active_element_ptr_range())
        {
          // First set descendents
          std::vector<Elem *> subactive_family;
          elem->total_family_tree(subactive_family);
          for (const auto & f : subactive_family)
            f->processor_id() = elem->processor_id();

          // Then set ancestors
          Elem * parent = elem->parent();

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
          std::vector<Elem *> subactive_family;
          child->total_family_tree(subactive_family);
          for (const auto & f : subactive_family)
            f->processor_id() = child->processor_id();
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
                  for (const auto & f : active_family)
                    parent_pid = std::min (parent_pid, f->processor_id());

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

void
Partitioner::processor_pairs_to_interface_nodes(MeshBase & mesh,
                                                std::map<std::pair<processor_id_type, processor_id_type>, std::set<dof_id_type>> & processor_pair_to_nodes)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  processor_pair_to_nodes.clear();

  std::set<dof_id_type> mynodes;
  std::set<dof_id_type> neighbor_nodes;
  std::vector<dof_id_type> common_nodes;

  // Loop over all the active elements
  for (auto & elem : mesh.active_element_ptr_range())
    {
      libmesh_assert(elem);

      libmesh_assert_not_equal_to (elem->processor_id(), DofObject::invalid_processor_id);

      auto n_nodes = elem->n_nodes();

      // prepare data for this element
      mynodes.clear();
      neighbor_nodes.clear();
      common_nodes.clear();

      for (unsigned int inode = 0; inode < n_nodes; inode++)
        mynodes.insert(elem->node_id(inode));

      for (auto i : elem->side_index_range())
        {
          auto neigh = elem->neighbor_ptr(i);
          if (neigh && !neigh->is_remote() && neigh->processor_id() != elem->processor_id())
            {
              neighbor_nodes.clear();
              common_nodes.clear();
              auto neigh_n_nodes = neigh->n_nodes();
              for (unsigned int inode = 0; inode < neigh_n_nodes; inode++)
                neighbor_nodes.insert(neigh->node_id(inode));

              std::set_intersection(mynodes.begin(), mynodes.end(),
                                    neighbor_nodes.begin(), neighbor_nodes.end(),
                                    std::back_inserter(common_nodes));

              auto & map_set = processor_pair_to_nodes[std::make_pair(std::min(elem->processor_id(), neigh->processor_id()),
                                                                      std::max(elem->processor_id(), neigh->processor_id()))];
              for (auto global_node_id : common_nodes)
                map_set.insert(global_node_id);
            }
        }
    }
}

void Partitioner::set_interface_node_processor_ids_linear(MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  std::map<std::pair<processor_id_type, processor_id_type>, std::set<dof_id_type>> processor_pair_to_nodes;

  processor_pairs_to_interface_nodes(mesh, processor_pair_to_nodes);

  for (auto & pmap : processor_pair_to_nodes)
    {
      std::size_t n_own_nodes = pmap.second.size()/2, i = 0;

      for (dof_id_type id : pmap.second)
        {
          auto & node = mesh.node_ref(id);
          if (i <= n_own_nodes)
            node.processor_id() = pmap.first.first;
          else
            node.processor_id() = pmap.first.second;
          i++;
        }
    }
}

void Partitioner::set_interface_node_processor_ids_BFS(MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  std::map<std::pair<processor_id_type, processor_id_type>, std::set<dof_id_type>> processor_pair_to_nodes;

  processor_pairs_to_interface_nodes(mesh, processor_pair_to_nodes);

  std::unordered_map<dof_id_type, std::vector<const Elem *>> nodes_to_elem_map;

  MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

  std::vector<const Node *>  neighbors;
  std::set<dof_id_type> neighbors_order;
  std::vector<dof_id_type> common_nodes;
  std::queue<dof_id_type> nodes_queue;
  std::set<dof_id_type> visted_nodes;

  for (auto & pmap : processor_pair_to_nodes)
    {
      std::size_t n_own_nodes = pmap.second.size()/2;

      // Initialize node assignment
      for (dof_id_type id : pmap.second)
        mesh.node_ref(id).processor_id() = pmap.first.second;

      visted_nodes.clear();
      for (dof_id_type id : pmap.second)
        {
          mesh.node_ref(id).processor_id() = pmap.first.second;

          if (visted_nodes.find(id) != visted_nodes.end())
            continue;
          else
            {
              nodes_queue.push(id);
              visted_nodes.insert(id);
              if (visted_nodes.size() >= n_own_nodes)
                break;
            }

          while (!nodes_queue.empty())
            {
              auto & node = mesh.node_ref(nodes_queue.front());
              nodes_queue.pop();

              neighbors.clear();
              MeshTools::find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbors);
              neighbors_order.clear();
              for (auto & neighbor : neighbors)
                neighbors_order.insert(neighbor->id());

              common_nodes.clear();
              std::set_intersection(pmap.second.begin(), pmap.second.end(),
                                    neighbors_order.begin(), neighbors_order.end(),
                                    std::back_inserter(common_nodes));

              for (auto c_node : common_nodes)
                if (visted_nodes.find(c_node) == visted_nodes.end())
                  {
                    nodes_queue.push(c_node);
                    visted_nodes.insert(c_node);
                    if (visted_nodes.size() >= n_own_nodes)
                      goto queue_done;
                  }

              if (visted_nodes.size() >= n_own_nodes)
                goto queue_done;
            }
        }
    queue_done:
      for (auto node : visted_nodes)
        mesh.node_ref(node).processor_id() = pmap.first.first;
    }
}

void Partitioner::set_interface_node_processor_ids_petscpartitioner(MeshBase & mesh)
{
  libmesh_ignore(mesh); // Only used if LIBMESH_HAVE_PETSC

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

#if LIBMESH_HAVE_PETSC
  std::map<std::pair<processor_id_type, processor_id_type>, std::set<dof_id_type>> processor_pair_to_nodes;

  processor_pairs_to_interface_nodes(mesh, processor_pair_to_nodes);

  std::vector<std::vector<const Elem *>> nodes_to_elem_map;

  MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

  std::vector<const Node *>  neighbors;
  std::set<dof_id_type> neighbors_order;
  std::vector<dof_id_type> common_nodes;

  std::vector<dof_id_type> rows;
  std::vector<dof_id_type> cols;

  std::map<dof_id_type, dof_id_type> global_to_local;

  for (auto & pmap : processor_pair_to_nodes)
    {
      unsigned int i = 0;

      rows.clear();
      rows.resize(pmap.second.size()+1);
      cols.clear();
      for (dof_id_type id : pmap.second)
        global_to_local[id] = i++;

      i = 0;
      for (auto id : pmap.second)
        {
          auto & node = mesh.node_ref(id);
          neighbors.clear();
          MeshTools::find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbors);
          neighbors_order.clear();
          for (auto & neighbor : neighbors)
            neighbors_order.insert(neighbor->id());

          common_nodes.clear();
          std::set_intersection(pmap.second.begin(), pmap.second.end(),
                                neighbors_order.begin(), neighbors_order.end(),
                                std::back_inserter(common_nodes));

          rows[i+1] = rows[i] + cast_int<dof_id_type>(common_nodes.size());

          for (auto c_node : common_nodes)
            cols.push_back(global_to_local[c_node]);

          i++;
        }

      Mat adj;
      MatPartitioning part;
      IS is;
      PetscInt local_size, rows_size, cols_size;
      PetscInt *adj_i, *adj_j;
      const PetscInt *indices;
      PetscCalloc1(rows.size(), &adj_i);
      PetscCalloc1(cols.size(), &adj_j);
      rows_size = cast_int<PetscInt>(rows.size());
      for (PetscInt ii=0; ii<rows_size; ii++)
        adj_i[ii] = rows[ii];

      cols_size = cast_int<PetscInt>(cols.size());
      for (PetscInt ii=0; ii<cols_size; ii++)
        adj_j[ii] = cols[ii];

      const PetscInt sz = cast_int<PetscInt>(pmap.second.size());
      MatCreateMPIAdj(PETSC_COMM_SELF, sz, sz, adj_i, adj_j,nullptr,&adj);
      MatPartitioningCreate(PETSC_COMM_SELF,&part);
      MatPartitioningSetAdjacency(part,adj);
      MatPartitioningSetNParts(part,2);
      PetscObjectSetOptionsPrefix((PetscObject)part, "balance_");
      MatPartitioningSetFromOptions(part);
      MatPartitioningApply(part,&is);

      MatDestroy(&adj);
      MatPartitioningDestroy(&part);

      ISGetLocalSize(is, &local_size);
      ISGetIndices(is, &indices);
      i = 0;
      for (auto id : pmap.second)
        {
          auto & node = mesh.node_ref(id);
          if (indices[i])
            node.processor_id() = pmap.first.second;
          else
            node.processor_id() = pmap.first.first;

          i++;
        }
      ISRestoreIndices(is, &indices);
      ISDestroy(&is);
    }
#else
  libmesh_error_msg("PETSc is required");
#endif
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

  // Start from scratch here: nodes we used to own may not be
  // eligible for us to own any more.
  for (auto & node : mesh.node_ptr_range())
    {
      node->processor_id() = DofObject::invalid_processor_id;
    }

  // Loop over all the active elements
  for (auto & elem : mesh.active_element_ptr_range())
    {
      libmesh_assert(elem);

      libmesh_assert_not_equal_to (elem->processor_id(), DofObject::invalid_processor_id);

      // Consider updating the processor id on this element's nodes
      for (Node & node : elem->node_ref_range())
        {
          processor_id_type & pid = node.processor_id();
          pid = node.choose_processor_id(pid, elem->processor_id());
        }
    }

  // How we finish off the node partitioning depends on our command
  // line options.

  const bool load_balanced_nodes_linear =
      libMesh::on_command_line ("--load-balanced-nodes-linear");

  const bool load_balanced_nodes_bfs =
       libMesh::on_command_line ("--load-balanced-nodes-bfs");

  const bool load_balanced_nodes_petscpartition =
      libMesh::on_command_line ("--load-balanced-nodes-petscpartitioner");

  unsigned int n_load_balance_options = load_balanced_nodes_linear;
  n_load_balance_options += load_balanced_nodes_bfs;
  n_load_balance_options += load_balanced_nodes_petscpartition;
  libmesh_error_msg_if(n_load_balance_options > 1,
                       "Cannot perform more than one load balancing type at a time");

  if (load_balanced_nodes_linear)
    set_interface_node_processor_ids_linear(mesh);
  else if (load_balanced_nodes_bfs)
    set_interface_node_processor_ids_BFS(mesh);
  else if (load_balanced_nodes_petscpartition)
    set_interface_node_processor_ids_petscpartitioner(mesh);

   // Node balancing algorithm will response to assign owned nodes.
   // We still need to sync PIDs
    {
      // For inactive elements, we will have already gotten most of
      // these nodes, *except* for the case of a parent with a subset
      // of active descendants which are remote elements.  In that
      // case some of the parent nodes will not have been properly
      // handled yet on our processor.
      //
      // We don't want to inadvertently give one of them an incorrect
      // processor id, but if we're not in serial then we have to
      // assign them temporary pids to make querying work, so we'll
      // save our *valid* pids before assigning temporaries.
      //
      // Even in serial we'll want to check and make sure we're not
      // overwriting valid active node pids with pids from subactive
      // elements.
      std::unordered_set<dof_id_type> bad_pids;

      for (auto & node : mesh.node_ptr_range())
        if (node->processor_id() == DofObject::invalid_processor_id)
          bad_pids.insert(node->id());

      // If we assign our temporary ids by looping from finer elements
      // to coarser elements, we'll always get an id from the finest
      // ghost element we can see, which will usually be "closer" to
      // the true processor we want to query and so will reduce query
      // cycles that don't reach that processor.

      // But we can still end up with a query cycle that dead-ends, so
      // we need to prepare a "push" communication step here.

      const bool is_serial = mesh.is_serial();
      std::unordered_map
        <processor_id_type,
         std::unordered_map<dof_id_type, processor_id_type>>
        potential_pids;

      const unsigned int n_levels = MeshTools::n_levels(mesh);
      for (unsigned int level = n_levels; level > 0; --level)
        {
          for (auto & elem : as_range(mesh.level_elements_begin(level-1),
                                      mesh.level_elements_end(level-1)))
            {
              libmesh_assert_not_equal_to (elem->processor_id(),
                                           DofObject::invalid_processor_id);

              const processor_id_type elem_pid = elem->processor_id();

              // Consider updating the processor id on this element's nodes
              for (Node & node : elem->node_ref_range())
                {
                  processor_id_type & pid = node.processor_id();
                  if (bad_pids.count(node.id()))
                    pid = node.choose_processor_id(pid, elem_pid);
                  else if (!is_serial)
                    potential_pids[elem_pid][node.id()] = pid;
                }
            }
        }

      if (!is_serial)
        {
          std::unordered_map
            <processor_id_type,
             std::vector<std::pair<dof_id_type, processor_id_type>>>
            potential_pids_vecs;

          for (auto & pair : potential_pids)
            potential_pids_vecs[pair.first].assign(pair.second.begin(), pair.second.end());

          auto pids_action_functor =
            [& mesh, & bad_pids]
            (processor_id_type /* src_pid */,
             const std::vector<std::pair<dof_id_type, processor_id_type>> & data)
            {
              for (auto pair : data)
                {
                  Node & node = mesh.node_ref(pair.first);
                  processor_id_type & pid = node.processor_id();
                  auto it = bad_pids.find(pair.first);
                  if (it != bad_pids.end())
                    {
                      pid = pair.second;
                      bad_pids.erase(it);
                    }
                  else
                    pid = node.choose_processor_id(pid, pair.second);
                }
            };

          Parallel::push_parallel_vector_data
            (mesh.comm(), potential_pids_vecs, pids_action_functor);

          // Using default libMesh options, we'll just need to sync
          // between processors now.  The catch here is that we can't
          // initially trust Node::choose_processor_id() because some
          // of those node processor ids are the temporary ones.
          CorrectProcIds correct_pids(mesh, bad_pids);
          Parallel::sync_node_data_by_element_id
            (mesh, mesh.elements_begin(), mesh.elements_end(),
             Parallel::SyncEverything(), Parallel::SyncEverything(),
             correct_pids);

          // But once we've got all the non-temporary pids synced, we
          // may need to sync again to get any pids on nodes only
          // connected to subactive elements, for which *only*
          // "temporary" pids are possible.
          bad_pids.clear();
          Parallel::sync_node_data_by_element_id
            (mesh,
             mesh.elements_begin(), mesh.elements_end(),
             Parallel::SyncEverything(), Parallel::SyncEverything(),
             correct_pids);
        }
    }

  // We can't assert that all nodes are connected to elements, because
  // a DistributedMesh with NodeConstraints might have pulled in some
  // remote nodes solely for evaluating those constraints.
  // MeshTools::libmesh_assert_connected_nodes(mesh);

#ifdef DEBUG
  MeshTools::libmesh_assert_valid_procids<Node>(mesh);
  //MeshTools::libmesh_assert_canonical_node_procids(mesh);
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

    for (auto i : index_range(ids))
      local_ids[i] = id_map[ids[i]];
  }

  void act_on_data (const std::vector<dof_id_type> & ids,
                    const std::vector<datum> & local_ids)
  {
    for (auto i : index_range(local_ids))
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
  auto n_proc = mesh.n_processors();
  _n_active_elem_on_proc.resize(n_proc);
  mesh.comm().allgather(n_active_local_elem, _n_active_elem_on_proc);

  std::vector<dof_id_type> n_active_elem_before_proc(mesh.n_processors());

  for (auto i : make_range(n_proc-1))
    n_active_elem_before_proc[i+1] =
      n_active_elem_before_proc[i] + _n_active_elem_on_proc[i];

  libMesh::BoundingBox bbox =
    MeshTools::create_bounding_box(mesh);

  _global_index_by_pid_map.clear();

  // create the mapping which is contiguous by processor
  MeshCommunication().find_local_indices (bbox,
                                          mesh.active_local_elements_begin(),
                                          mesh.active_local_elements_end(),
                                          _global_index_by_pid_map);

  SyncLocalIDs sync(_global_index_by_pid_map);

  Parallel::sync_dofobject_data_by_id
      (mesh.comm(), mesh.active_elements_begin(), mesh.active_elements_end(), sync);

  for (const auto & elem : mesh.active_element_ptr_range())
    {
      const processor_id_type pid = elem->processor_id();
      libmesh_assert_less (_global_index_by_pid_map[elem->id()], _n_active_elem_on_proc[pid]);

      _global_index_by_pid_map[elem->id()] += n_active_elem_before_proc[pid];
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

      for (const auto & neighbor : neighbor_set)
        interior_to_boundary_map.emplace(neighbor, elem);
    }

#ifdef LIBMESH_ENABLE_AMR
  std::vector<const Elem *> neighbors_offspring;
#endif

  // This is costly, and we only need to do it if the mesh has
  // changed since we last partitioned... but the mesh probably has
  // changed since we last partitioned, and if it hasn't we don't
  // have a reliable way to be sure of that.
  _find_global_index_by_pid_map(mesh);

  dof_id_type first_local_elem = 0;
  for (auto pid : make_range(mesh.processor_id()))
     first_local_elem += _n_active_elem_on_proc[pid];

  _dual_graph.clear();
  _dual_graph.resize(n_active_local_elem);
  _local_id_to_elem.resize(n_active_local_elem);

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      libmesh_assert (_global_index_by_pid_map.count(elem->id()));
      const dof_id_type global_index_by_pid =
        _global_index_by_pid_map[elem->id()];

      const dof_id_type local_index =
        global_index_by_pid - first_local_elem;
      libmesh_assert_less (local_index, n_active_local_elem);

      std::vector<dof_id_type> & graph_row = _dual_graph[local_index];

      // Save this off to make it easy to index later
      _local_id_to_elem[local_index] = elem;

      // Loop over the element's neighbors.  An element
      // adjacency corresponds to a face neighbor
      for (auto neighbor : elem->neighbor_ptr_range())
        {
          if (neighbor != nullptr)
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
                  for (const auto & child : neighbors_offspring)
                    {
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

          for (const auto & neighbor : neighbor_set)
            {
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
  for (auto pid : make_range(mesh.processor_id()))
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
  const processor_id_type * ex = nullptr;
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
