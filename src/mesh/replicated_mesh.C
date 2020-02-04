// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/replicated_mesh_impl.h"

namespace libMesh
{
template <>
void ReplicatedMeshTempl<Real>::stitching_helper (const ReplicatedMesh * other_mesh,
                                                  boundary_id_type this_mesh_boundary_id,
                                                  boundary_id_type other_mesh_boundary_id,
                                                  Real tol,
                                                  bool clear_stitched_boundary_ids,
                                                  bool verbose,
                                                  bool use_binary_search,
                                                  bool enforce_all_nodes_match_on_boundaries,
                                                  bool skip_find_neighbors)
{
  std::map<dof_id_type, dof_id_type> node_to_node_map, other_to_this_node_map; // The second is the inverse map of the first
  std::map<dof_id_type, std::vector<dof_id_type>> node_to_elems_map;

  typedef dof_id_type                     key_type;
  typedef std::pair<Elem *, unsigned char> val_type;
  typedef std::pair<key_type, val_type>   key_val_pair;
  typedef std::unordered_multimap<key_type, val_type> map_type;
  // Mapping between all side keys in this mesh and elements+side numbers relevant to the boundary in this mesh as well.
  map_type side_to_elem_map;

  // If there is only one mesh (i.e. other_mesh == nullptr), then loop over this mesh twice
  if (!other_mesh)
    {
      other_mesh = this;
    }

  if ((this_mesh_boundary_id  != BoundaryInfo::invalid_id) &&
      (other_mesh_boundary_id != BoundaryInfo::invalid_id))
    {
      LOG_SCOPE("stitch_meshes node merging", "ReplicatedMesh");

      // While finding nodes on the boundary, also find the minimum edge length
      // of all faces on both boundaries.  This will later be used in relative
      // distance checks when stitching nodes.
      Real h_min = std::numeric_limits<Real>::max();
      bool h_min_updated = false;

      // Loop below fills in these sets for the two meshes.
      std::set<dof_id_type> this_boundary_node_ids, other_boundary_node_ids;

      // Pull objects out of the loop to reduce heap operations
      std::unique_ptr<Elem> side;

      {
        // Make temporary fixed-size arrays for loop
        boundary_id_type id_array[2]         = {this_mesh_boundary_id, other_mesh_boundary_id};
        std::set<dof_id_type> * set_array[2] = {&this_boundary_node_ids, &other_boundary_node_ids};
        const ReplicatedMesh * mesh_array[2] = {this, other_mesh};

        for (unsigned i=0; i<2; ++i)
          {
            // First we deal with node boundary IDs.
            // We only enter this loop if we have at least one
            // nodeset.
            if (mesh_array[i]->get_boundary_info().n_nodeset_conds() > 0)
              {
                // build_node_list() returns a vector of (node-id, bc-id) tuples
                for (const auto & t : mesh_array[i]->get_boundary_info().build_node_list())
                  {
                    boundary_id_type node_bc_id = std::get<1>(t);
                    if (node_bc_id == id_array[i])
                      {
                        dof_id_type this_node_id = std::get<0>(t);
                        set_array[i]->insert( this_node_id );

                        // We need to set h_min to some value. It's too expensive to
                        // search for the element that actually contains this node,
                        // since that would require a PointLocator. As a result, we
                        // just use the first (non-NodeElem!) element in the mesh to
                        // give us hmin if it's never been set before.
                        if (!h_min_updated)
                          {
                            for (const auto & elem : mesh_array[i]->active_element_ptr_range())
                              {
                                Real current_h_min = elem->hmin();
                                if (current_h_min > 0.)
                                  {
                                    h_min = current_h_min;
                                    h_min_updated = true;
                                    break;
                                  }
                              }

                            // If, after searching all the active elements, we did not update
                            // h_min, give up and set h_min to 1 so that we don't repeat this
                            // fruitless search
                            if (!h_min_updated)
                              {
                                h_min_updated = true;
                                h_min = 1.0;
                              }
                          }
                      }
                  }
              }

            // Container to catch boundary IDs passed back from BoundaryInfo.
            std::vector<boundary_id_type> bc_ids;

            for (auto & el : mesh_array[i]->element_ptr_range())
              {
                // Now check whether elem has a face on the specified boundary
                for (auto side_id : el->side_index_range())
                  if (el->neighbor_ptr(side_id) == nullptr)
                    {
                      // Get *all* boundary IDs on this side, not just the first one!
                      mesh_array[i]->get_boundary_info().boundary_ids (el, side_id, bc_ids);

                      if (std::find(bc_ids.begin(), bc_ids.end(), id_array[i]) != bc_ids.end())
                        {
                          el->build_side_ptr(side, side_id);
                          for (auto & n : side->node_ref_range())
                            set_array[i]->insert(n.id());

                          h_min = std::min(h_min, side->hmin());
                          h_min_updated = true;

                          // This side is on the boundary, add its information to side_to_elem
                          if (skip_find_neighbors && (i==0))
                            {
                              key_type key = el->key(side_id);
                              val_type val;
                              val.first = el;
                              val.second = cast_int<unsigned char>(side_id);

                              key_val_pair kvp;
                              kvp.first = key;
                              kvp.second = val;
                              side_to_elem_map.insert (kvp);
                            }
                        }

                      // Also, check the edges on this side. We don't have to worry about
                      // updating neighbor info in this case since elements don't store
                      // neighbor info on edges.
                      for (auto edge_id : el->edge_index_range())
                        {
                          if (el->is_edge_on_side(edge_id, side_id))
                            {
                              // Get *all* boundary IDs on this edge, not just the first one!
                              mesh_array[i]->get_boundary_info().edge_boundary_ids (el, edge_id, bc_ids);

                              if (std::find(bc_ids.begin(), bc_ids.end(), id_array[i]) != bc_ids.end())
                                {
                                  std::unique_ptr<Elem> edge (el->build_edge_ptr(edge_id));
                                  for (auto & n : edge->node_ref_range())
                                    set_array[i]->insert( n.id() );

                                  h_min = std::min(h_min, edge->hmin());
                                  h_min_updated = true;
                                }
                            }
                        }
                    }
              }
          }
      }

      if (verbose)
        {
          libMesh::out << "In ReplicatedMesh::stitch_meshes:\n"
                       << "This mesh has "  << this_boundary_node_ids.size()
                       << " nodes on boundary " << this_mesh_boundary_id  << ".\n"
                       << "Other mesh has " << other_boundary_node_ids.size()
                       << " nodes on boundary " << other_mesh_boundary_id << ".\n";

          if (h_min_updated)
            {
              libMesh::out << "Minimum edge length on both surfaces is " << h_min << ".\n";
            }
          else
            {
              libMesh::out << "No elements on specified surfaces." << std::endl;
            }
        }

      // We require nanoflann for the "binary search" (really kd-tree)
      // option to work. If it's not available, turn that option off,
      // warn the user, and fall back on the N^2 search algorithm.
      if (use_binary_search)
        {
#ifndef LIBMESH_HAVE_NANOFLANN
          use_binary_search = false;
          libmesh_warning("The use_binary_search option in the "
                          "ReplicatedMesh stitching algorithms requires nanoflann "
                          "support. Falling back on N^2 search algorithm.");
#endif
        }

      if (this_boundary_node_ids.size())
      {
        if (use_binary_search)
        {
#ifdef LIBMESH_HAVE_NANOFLANN
          typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<Real, VectorOfNodesAdaptor<Real>>, VectorOfNodesAdaptor<Real>, 3> kd_tree_t;

          // Create the dataset needed to build the kd tree with nanoflann
          std::vector<std::pair<Point, dof_id_type>> this_mesh_nodes(this_boundary_node_ids.size());
          std::set<dof_id_type>::iterator current_node = this_boundary_node_ids.begin(),
                                          node_ids_end = this_boundary_node_ids.end();
          for (unsigned int ctr = 0; current_node != node_ids_end; ++current_node, ++ctr)
          {
            this_mesh_nodes[ctr].first = this->point(*current_node);
            this_mesh_nodes[ctr].second = *current_node;
          }

          VectorOfNodesAdaptor<Real> vec_nodes_adaptor(this_mesh_nodes);

          kd_tree_t this_kd_tree(3, vec_nodes_adaptor, 10);
          this_kd_tree.buildIndex();

          // Storage for nearest neighbor in the loop below
          std::vector<size_t> ret_index(1);
          std::vector<Real> ret_dist_sqr(1);

          // Loop over other mesh. For each node, find its nearest neighbor in this mesh, and fill in the maps.
          for (auto node : other_boundary_node_ids)
          {
            const Real query_pt[] = {other_mesh->point(node)(0), other_mesh->point(node)(1), other_mesh->point(node)(2)};
            this_kd_tree.knnSearch(&query_pt[0], 1, &ret_index[0], &ret_dist_sqr[0]);
            if (ret_dist_sqr[0] < TOLERANCE*TOLERANCE)
            {
              node_to_node_map[this_mesh_nodes[ret_index[0]].second] = node;
              other_to_this_node_map[node] = this_mesh_nodes[ret_index[0]].second;
            }
          }

          // If the 2 maps don't have the same size, it means we have overwritten a value in node_to_node_map
          // It means one node in this mesh is the nearest neighbor of several nodes in other mesh.
          // Not possible !
          if (node_to_node_map.size() != other_to_this_node_map.size())
            libmesh_error_msg("Error: Found multiple matching nodes in stitch_meshes");
#endif
        }
        else
        {
          // In the unlikely event that two meshes composed entirely of
          // NodeElems are being stitched together, we will not have
          // selected a valid h_min value yet, and the distance
          // comparison below will be true for essentially any two
          // nodes. In this case we simply fall back on an absolute
          // distance check.
          if (!h_min_updated)
            {
              libmesh_warning("No valid h_min value was found, falling back on "
                              "absolute distance check in the N^2 search algorithm.");
              h_min = 1.;
            }

          // Otherwise, use a simple N^2 search to find the closest matching points. This can be helpful
          // in the case that we have tolerance issues which cause mismatch between the two surfaces
          // that are being stitched.
          for (const auto & this_node_id : this_boundary_node_ids)
          {
            Node & this_node = this->node_ref(this_node_id);

            bool found_matching_nodes = false;

            for (const auto & other_node_id : other_boundary_node_ids)
            {
              const Node & other_node = other_mesh->node_ref(other_node_id);

              Real node_distance = (this_node - other_node).norm();

              if (node_distance < tol*h_min)
              {
                // Make sure we didn't already find a matching node!
                if (found_matching_nodes)
                  libmesh_error_msg("Error: Found multiple matching nodes in stitch_meshes");

                node_to_node_map[this_node_id] = other_node_id;
                other_to_this_node_map[other_node_id] = this_node_id;

                found_matching_nodes = true;
              }
            }
          }
        }
      }

      // Build up the node_to_elems_map, using only one loop over other_mesh
      for (auto & el : other_mesh->element_ptr_range())
        {
          // For each node on the element, find the corresponding node
          // on "this" Mesh, 'this_node_id', if it exists, and push
          // the current element ID back onto node_to_elems_map[this_node_id].
          // For that we will use the reverse mapping we created at
          // the same time as the forward mapping.
          for (auto & n : el->node_ref_range())
            {
              dof_id_type other_node_id = n.id();
              std::map<dof_id_type, dof_id_type>::iterator it =
                other_to_this_node_map.find(other_node_id);

              if (it != other_to_this_node_map.end())
                {
                  dof_id_type this_node_id = it->second;
                  node_to_elems_map[this_node_id].push_back( el->id() );
                }
            }
        }

      if (verbose)
        {
          libMesh::out << "In ReplicatedMesh::stitch_meshes:\n"
                       << "Found " << node_to_node_map.size()
                       << " matching nodes.\n"
                       << std::endl;
        }

      if (enforce_all_nodes_match_on_boundaries)
        {
          std::size_t n_matching_nodes = node_to_node_map.size();
          std::size_t this_mesh_n_nodes = this_boundary_node_ids.size();
          std::size_t other_mesh_n_nodes = other_boundary_node_ids.size();
          if ((n_matching_nodes != this_mesh_n_nodes) || (n_matching_nodes != other_mesh_n_nodes))
            libmesh_error_msg("Error: We expected the number of nodes to match.");
        }
    }
  else
    {
      if (verbose)
        {
          libMesh::out << "Skip node merging in ReplicatedMesh::stitch_meshes:" << std::endl;
        }
    }

  dof_id_type node_delta = this->max_node_id();
  dof_id_type elem_delta = this->max_elem_id();

  unique_id_type unique_delta =
#ifdef LIBMESH_ENABLE_UNIQUE_ID
    this->parallel_max_unique_id();
#else
    0;
#endif

  // If other_mesh != nullptr, then we have to do a bunch of work
  // in order to copy it to this mesh
  if (this!=other_mesh)
    {
      LOG_SCOPE("stitch_meshes copying", "ReplicatedMesh");

      // Increment the node_to_node_map and node_to_elems_map
      // to account for id offsets
      for (auto & pr : node_to_node_map)
        pr.second += node_delta;

      for (auto & pr : node_to_elems_map)
        for (auto & entry : pr.second)
          entry += elem_delta;

      // Copy mesh data. If we skip the call to find_neighbors(), the lists
      // of neighbors will be copied verbatim from the other mesh
      this->copy_nodes_and_elements(*other_mesh, skip_find_neighbors,
                                    elem_delta, node_delta,
                                    unique_delta);

      // Copy BoundaryInfo from other_mesh too.  We do this via the
      // list APIs rather than element-by-element for speed.
      BoundaryInfo & boundary = this->get_boundary_info();
      const BoundaryInfo & other_boundary = other_mesh->get_boundary_info();

      for (const auto & t : other_boundary.build_node_list())
        boundary.add_node(std::get<0>(t) + node_delta,
                          std::get<1>(t));

      for (const auto & t : other_boundary.build_side_list())
        boundary.add_side(std::get<0>(t) + elem_delta,
                          std::get<1>(t),
                          std::get<2>(t));

      for (const auto & t : other_boundary.build_edge_list())
        boundary.add_edge(std::get<0>(t) + elem_delta,
                          std::get<1>(t),
                          std::get<2>(t));

      for (const auto & t : other_boundary.build_shellface_list())
        boundary.add_shellface(std::get<0>(t) + elem_delta,
                               std::get<1>(t),
                               std::get<2>(t));

    } // end if (other_mesh)

  // Finally, we need to "merge" the overlapping nodes
  // We do this by iterating over node_to_elems_map and updating
  // the elements so that they "point" to the nodes that came
  // from this mesh, rather than from other_mesh.
  // Then we iterate over node_to_node_map and delete the
  // duplicate nodes that came from other_mesh.

  {
    LOG_SCOPE("stitch_meshes node updates", "ReplicatedMesh");

    // Container to catch boundary IDs passed back from BoundaryInfo.
    std::vector<boundary_id_type> bc_ids;

    for (const auto & pr : node_to_elems_map)
      {
        dof_id_type target_node_id = pr.first;
        dof_id_type other_node_id = node_to_node_map[target_node_id];
        Node & target_node = this->node_ref(target_node_id);

        std::size_t n_elems = pr.second.size();
        for (std::size_t i=0; i<n_elems; i++)
          {
            dof_id_type elem_id = pr.second[i];
            Elem * el = this->elem_ptr(elem_id);

            // find the local node index that we want to update
            unsigned int local_node_index = el->local_node(other_node_id);
            libmesh_assert_not_equal_to(local_node_index, libMesh::invalid_uint);

            // We also need to copy over the nodeset info here,
            // because the node will get deleted below
            this->get_boundary_info().boundary_ids(el->node_ptr(local_node_index), bc_ids);
            el->set_node(local_node_index) = &target_node;
            this->get_boundary_info().add_node(&target_node, bc_ids);
          }
      }
  }

  {
    LOG_SCOPE("stitch_meshes node deletion", "ReplicatedMesh");
    for (const auto & pr : node_to_node_map)
      {
        // In the case that this==other_mesh, the two nodes might be the same (e.g. if
        // we're stitching a "sliver"), hence we need to skip node deletion in that case.
        if ((this == other_mesh) && (pr.second == pr.first))
          continue;

        dof_id_type this_node_id = pr.second;
        this->delete_node( this->node_ptr(this_node_id) );
      }
  }

  // If find_neighbors() wasn't called in prepare_for_use(), we need to
  // manually loop once more over all elements adjacent to the stitched boundary
  // and fix their lists of neighbors.
  // This is done according to the following steps:
  //   1. Loop over all copied elements adjacent to the boundary using node_to_elems_map (trying to avoid duplicates)
  //   2. Look at all their sides with a nullptr neighbor and update them using side_to_elem_map if necessary
  //   3. Update the corresponding side in side_to_elem_map as well
  if (skip_find_neighbors)
    {
      LOG_SCOPE("stitch_meshes neighbor fixes", "ReplicatedMesh");

      // Pull objects out of the loop to reduce heap operations
      std::unique_ptr<Elem> my_side, their_side;

      std::set<dof_id_type> fixed_elems;
      for (const auto & pr : node_to_elems_map)
        {
          std::size_t n_elems = pr.second.size();
          for (std::size_t i=0; i<n_elems; i++)
            {
              dof_id_type elem_id = pr.second[i];
              if (fixed_elems.find(elem_id) == fixed_elems.end())
                {
                  Elem * el = this->elem_ptr(elem_id);
                  fixed_elems.insert(elem_id);
                  for (auto s : el->side_index_range())
                    {
                      if (el->neighbor_ptr(s) == nullptr)
                        {
                          key_type key = el->key(s);
                          auto bounds = side_to_elem_map.equal_range(key);

                          if (bounds.first != bounds.second)
                            {
                              // Get the side for this element
                              el->side_ptr(my_side, s);

                              // Look at all the entries with an equivalent key
                              while (bounds.first != bounds.second)
                                {
                                  // Get the potential element
                                  Elem * neighbor = bounds.first->second.first;

                                  // Get the side for the neighboring element
                                  const unsigned int ns = bounds.first->second.second;
                                  neighbor->side_ptr(their_side, ns);
                                  //libmesh_assert(my_side.get());
                                  //libmesh_assert(their_side.get());

                                  // If found a match with my side
                                  //
                                  // We need special tests here for 1D:
                                  // since parents and children have an equal
                                  // side (i.e. a node), we need to check
                                  // ns != ms, and we also check level() to
                                  // avoid setting our neighbor pointer to
                                  // any of our neighbor's descendants
                                  if ((*my_side == *their_side) &&
                                      (el->level() == neighbor->level()) &&
                                      ((el->dim() != 1) || (ns != s)))
                                    {
                                      // So share a side.  Is this a mixed pair
                                      // of subactive and active/ancestor
                                      // elements?
                                      // If not, then we're neighbors.
                                      // If so, then the subactive's neighbor is

                                      if (el->subactive() ==
                                          neighbor->subactive())
                                        {
                                          // an element is only subactive if it has
                                          // been coarsened but not deleted
                                          el->set_neighbor (s,neighbor);
                                          neighbor->set_neighbor(ns,el);
                                        }
                                      else if (el->subactive())
                                        {
                                          el->set_neighbor(s,neighbor);
                                        }
                                      else if (neighbor->subactive())
                                        {
                                          neighbor->set_neighbor(ns,el);
                                        }
                                      // It's OK to invalidate the
                                      // bounds.first iterator here,
                                      // as we are immediately going
                                      // to break out of this while
                                      // loop. bounds.first will
                                      // therefore not be used for
                                      // anything else.
                                      side_to_elem_map.erase (bounds.first);
                                      break;
                                    }

                                  ++bounds.first;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  this->prepare_for_use( /*skip_renumber_nodes_and_elements= */ false, skip_find_neighbors);

  // After the stitching, we may want to clear boundary IDs from element
  // faces that are now internal to the mesh
  if (clear_stitched_boundary_ids)
    {
      LOG_SCOPE("stitch_meshes clear bcids", "ReplicatedMesh");

      // Container to catch boundary IDs passed back from BoundaryInfo.
      std::vector<boundary_id_type> bc_ids;

      for (auto & el : element_ptr_range())
        for (auto side_id : el->side_index_range())
          if (el->neighbor_ptr(side_id) != nullptr)
            {
              // Completely remove the side from the boundary_info object if it has either
              // this_mesh_boundary_id or other_mesh_boundary_id.
              this->get_boundary_info().boundary_ids (el, side_id, bc_ids);

              if (std::find(bc_ids.begin(), bc_ids.end(), this_mesh_boundary_id) != bc_ids.end() ||
                  std::find(bc_ids.begin(), bc_ids.end(), other_mesh_boundary_id) != bc_ids.end())
                this->get_boundary_info().remove_side(el, side_id);
            }

      // Removing stitched-away boundary ids might have removed an id
      // *entirely*, so we need to recompute boundary id sets to check
      // for that.
      this->get_boundary_info().regenerate_id_sets();
    }
}

template class ReplicatedMeshTempl<Real>;
}
