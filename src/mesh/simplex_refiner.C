// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/simplex_refiner.h"

#include "libmesh/boundary_info.h"
#include "libmesh/face_tri3.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/node.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/remote_elem.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/utility.h"

// TIMPI includes
#include "timpi/parallel_implementation.h"
#include "timpi/parallel_sync.h"


namespace libMesh
{

//-----------------------------------------------------------------
// Mesh refinement methods
SimplexRefiner::SimplexRefiner (UnstructuredMesh & mesh) :
  _mesh(mesh),
  _desired_volume(0.1)
{}



bool SimplexRefiner::refine_elements()
{
  if (!_mesh.is_prepared())
    _mesh.prepare_for_use();

  // We should add more options later: to refine via circumcenter
  // insertion instead of edge center insertion, to do Delaunay swaps,
  // etc.
  return this->refine_via_edges();
}



bool SimplexRefiner::should_refine_elem(Elem & elem)
{
  const Real min_volume_target = this->desired_volume();
  FunctionBase<Real> *volume_func = this->get_desired_volume_function();

  // If this isn't a question, why are we here?
  libmesh_assert(min_volume_target > 0 ||
                 volume_func != nullptr);

  // If we don't have position-dependent volume targets we can make a
  // decision quickly
  const Real volume = elem.volume();

  if (volume < TOLERANCE*TOLERANCE)
    libmesh_warning
      ("Warning: trying to refine sliver element with volume " <<
       volume << "!\n");

  if (!volume_func)
    return (volume > min_volume_target);

  // If we do?
  //
  // See if we're meeting the local volume target at all the elem
  // vertices first
  for (auto v : make_range(elem.n_vertices()))
    {
      // If we have an auto volume function, we'll use it and override other volume options
      const Real local_volume_target = (*volume_func)(elem.point(v));
      libmesh_error_msg_if
        (local_volume_target <= 0,
         "Non-positive desired element volumes are unachievable");
      if (volume > local_volume_target)
        return true;
    }

  // If our vertices are happy, it's still possible that our interior
  // isn't.  Are we allowed not to bother checking it?
  if (!min_volume_target)
    return false;

  libmesh_not_implemented_msg
    ("Combining a minimum desired_volume with an volume function isn't yet supported.");
}


std::size_t SimplexRefiner::refine_via_edges(Elem & elem,
                                             dof_id_type coarse_id)
{
  // Do we need to be refined for our own sake, or only if our
  // neighbors were refined?
  const bool should_refine = this->should_refine_elem(elem);

  // We only currently handle coarse conforming meshes here.
  // Uniformly refined meshes should be flattened before using this.
  libmesh_assert_equal_to(elem.level(), 0u);

  // Adding Tets or Tri6/Tri7 wouldn't be too much harder.  Quads or
  // other elements will get tricky.
  libmesh_assert_equal_to(elem.type(), TRI3);

  const int n_sides = elem.n_sides();
  int side_to_refine = 0;
  Real length_to_refine = 0;
  std::pair<Node *, Node *> vertices_to_refine;

  // Find the longest side and refine it first.  We might change that
  // heuristic at some point, to enable anisotropic refinement.
  for (int side : make_range(n_sides))
    {
      std::pair<Node *, Node *> vertices
        {elem.node_ptr((side+1)%n_sides),
         elem.node_ptr(side)};

      if ((Point&)(*vertices.first) >
          (Point&)(*vertices.second))
        std::swap(vertices.first, vertices.second);

      const auto sidevec = *(Point*)vertices.first -
                           *(Point*)vertices.second;
      Real length = sidevec.norm();

      // We might need to ask the owner of a coarse ghost element
      // about whether it saw any refinement forced by a coarse
      // neighbor with a smaller desired area.
      if (const processor_id_type pid = elem.processor_id();
          pid != _mesh.processor_id() && !_mesh.is_serial() &&
          !_desired_volume_func.get() &&
          elem.id() == coarse_id)
        edge_queries[pid].emplace_back(vertices.first->id(),
                                       vertices.second->id());

      // Test should_refine, but also test for any edges that were
      // already refined by neighbors, which might happen even if
      // !should_refine in cases where we're refining non-uniformly.
      if (length > length_to_refine &&
          (should_refine ||
           new_nodes.find(vertices) != new_nodes.end()))
        {
          side_to_refine = side;
          length_to_refine = length;
          vertices_to_refine = vertices;
        }
    }

  // We might not refine anything.
  if (!length_to_refine)
    return 0;

  // But if we must, then do it.
  //
  // Find all the edge neighbors we can.  For a triangle there'll just
  // be the one, but we want to extend to tets later.
  std::vector<Elem *> neighbors;
  if (Elem * neigh = elem.neighbor_ptr(side_to_refine); neigh)
    neighbors.push_back(neigh);

  Node * midedge_node;
  if (auto it = new_nodes.find(vertices_to_refine);
      it != new_nodes.end())
    {
      midedge_node = it->second;
    }
  else
    {
      // Add with our own processor id first, so we can request a temporary
      // id that won't conflict with anyone else's canonical ids later
      midedge_node =
        _mesh.add_point((*(Point*)vertices_to_refine.first +
                         *(Point*)vertices_to_refine.second)/2,
                        DofObject::invalid_id,
                        _mesh.processor_id());

      // Then give the new node the best pid we can.
      midedge_node->processor_id() = elem.processor_id();
      for (const Elem * neigh : neighbors)
        if (neigh != remote_elem)
          midedge_node->processor_id() = std::min(midedge_node->processor_id(),
                                                  neigh->processor_id());

      new_nodes[vertices_to_refine] = midedge_node;
    }

  // Good for Tris and Tets; we'll need more for Quads.
  constexpr int max_subelems = 2;
  std::unique_ptr<Tri3> subelem[max_subelems];

  // We'll switch(elem->type()) here eventually
  subelem[0] = std::make_unique<Tri3>();
  subelem[0]->set_node(0) = elem.node_ptr(side_to_refine);
  subelem[0]->set_node(1) = midedge_node;
  subelem[0]->set_node(2) = elem.node_ptr((side_to_refine+2)%3);

  subelem[1] = std::make_unique<Tri3>();
  subelem[1]->set_node(0) = elem.node_ptr((side_to_refine+2)%3);
  subelem[1]->set_node(1) = midedge_node;
  subelem[1]->set_node(2) = elem.node_ptr((side_to_refine+1)%3);

  // Preserve boundary and sort-of-preserve neighbor information.  We
  // need remote_elem neighbors to avoid breaking distributed meshes
  // (which can't reconstruct them from scratch), and we need at least
  // placeholder neighbors to make sure we'll have good processor_id()
  // options for new nodes in future splits of the same edge.
  subelem[0]->set_neighbor(1, subelem[1].get());
  subelem[1]->set_neighbor(0, subelem[0].get());

  BoundaryInfo & boundary_info = _mesh.get_boundary_info();
  std::vector<boundary_id_type> bc_ids;
  boundary_info.boundary_ids(&elem, side_to_refine, bc_ids);
  boundary_info.add_side(subelem[0].get(), 0, bc_ids);
  boundary_info.add_side(subelem[1].get(), 1, bc_ids);
  subelem[0]->set_neighbor(0, elem.neighbor_ptr(side_to_refine));
  subelem[1]->set_neighbor(1, elem.neighbor_ptr(side_to_refine));

  boundary_info.boundary_ids(&elem, (side_to_refine+1)%3, bc_ids);
  boundary_info.add_side(subelem[1].get(), 2, bc_ids);
  subelem[1]->set_neighbor(2, elem.neighbor_ptr((side_to_refine+1)%3));

  boundary_info.boundary_ids(&elem, (side_to_refine+2)%3, bc_ids);
  boundary_info.add_side(subelem[0].get(), 2, bc_ids);
  subelem[0]->set_neighbor(2, elem.neighbor_ptr((side_to_refine+2)%3));

  // Be sure the correct data is set for all subelems.
  const unsigned int nei = elem.n_extra_integers();
  for (unsigned int i=0; i != max_subelems; ++i)
    if (subelem[i])
      {
        // We'd love to assert that we don't have a sliver here,
        // but maybe our input had a sliver and we can't fix it.
        libmesh_assert_greater_equal(subelem[i]->volume() + TOLERANCE,
                                     elem.volume() / 2);
        subelem[i]->processor_id() = elem.processor_id();
        subelem[i]->subdomain_id() = elem.subdomain_id();

        // Copy any extra element data.  Since the subelements
        // haven't been added to the mesh yet any allocation has
        // to be done manually.
        subelem[i]->add_extra_integers(nei);
        for (unsigned int ei=0; ei != nei; ++ei)
          subelem[ei]->set_extra_integer(ei, elem.get_extra_integer(ei));

        // Copy any mapping data.
        subelem[i]->set_mapping_type(elem.mapping_type());
        subelem[i]->set_mapping_data(elem.mapping_data());
      }

  boundary_info.remove(&elem);

  if (auto it = this->added_elements.find(coarse_id);
      it != this->added_elements.end())
    {
      auto & added = it->second;
      if (auto sub_it = added.find(elem.vertex_average());
          sub_it != added.end())
        {
          if (&elem == sub_it->second)
            added.erase(sub_it);
        }
    }

  // We only add an element to new_elements *right* before
  // checking it for further refinement, so we only need to check the
  // back() to see where the elem here came from.
  if (!this->new_elements.empty() &&
      &elem == this->new_elements.back().get())
    {
      this->new_elements.pop_back();
      this->coarse_parent.erase(&elem);
    }
  else
    _mesh.delete_elem(&elem);

  // We refined at least ourselves
  std::size_t n_refined_elements = 1;

  for (unsigned int i=0; i != max_subelems; ++i)
    {
      Elem * add_elem = subelem[i].get();
      added_elements[coarse_id].emplace(add_elem->vertex_average(),
                                        add_elem);
      this->new_elements.push_back(std::move(subelem[i]));
      this->coarse_parent[add_elem] = coarse_id;
      n_refined_elements +=
        this->refine_via_edges(*add_elem, coarse_id);
    }

  return n_refined_elements;
}



std::size_t SimplexRefiner::refine_via_edges()
{
  this->new_nodes.clear();
  this->new_elements.clear();

  // We need to look at neighbors' pids to determine new node pids,
  // but we might be deleting elements and leaving dangling neighbor
  // pointers.  Put proxy neighbors in those pointers places instead.
  for (auto & elem : _mesh.element_ptr_range())
    for (auto n : make_range(elem->n_neighbors()))
      {
        Elem * neigh = elem->neighbor_ptr(n);
        if (!neigh || neigh == remote_elem)
          continue;

        const processor_id_type p = neigh->processor_id();
        std::unique_ptr<Elem> & proxy_neighbor = proxy_elements[p];
        if (!proxy_neighbor.get())
          proxy_neighbor = Elem::build(NODEELEM);

        elem->set_neighbor(n, proxy_neighbor.get());
      }

  // If we have a _desired_volume_func then we might still have
  // some unsplit edges that could have been split by their
  // remote_elem neighbors, and we'll need to query those.
  auto edge_gather_functor =
    [this]
    (processor_id_type,
     const std::vector<std::pair<dof_id_type, dof_id_type>> & edges,
     std::vector<refinement_datum> & edge_refinements)
    {
      // Fill those requests
      const std::size_t query_size = edges.size();
      edge_refinements.resize(query_size);
      for (std::size_t i=0; i != query_size; ++i)
        {
          auto vertex_ids = edges[i];
          std::pair<Node *, Node *> vertices
            {_mesh.node_ptr(vertex_ids.first),
             _mesh.node_ptr(vertex_ids.second)};
          fill_refinement_datum(vertices, edge_refinements[i]);
        }
    };

  auto edge_action_functor =
    [this]
    (processor_id_type,
     const std::vector<std::pair<dof_id_type, dof_id_type>> &,
     const std::vector<refinement_datum> & edge_refinements
    )
    {
      for (const auto & one_edges_refinements : edge_refinements)
        for (const auto & refinement : one_edges_refinements)
          {
            std::pair<Node *, Node *> vertices
              {_mesh.node_ptr(std::get<0>(refinement)),
               _mesh.node_ptr(std::get<1>(refinement))};

            if ((Point&)(*vertices.first) >
                (Point&)(*vertices.second))
              std::swap(vertices.first, vertices.second);

            if (auto it = new_nodes.find(vertices);
                it != new_nodes.end())
              {
                _mesh.renumber_node(it->second->id(),
                                    std::get<2>(refinement));
                it->second->processor_id() = std::get<3>(refinement);
              }
            else
              {
                Node * new_node =
                  _mesh.add_point((*(Point*)vertices.first +
                                   *(Point*)vertices.second)/2,
                                  std::get<2>(refinement),
                                  std::get<3>(refinement));
                new_nodes[vertices] = new_node;
              }
          }
    };

  std::size_t refined_elements = 0;
  std::size_t newly_refined_elements = 0;
  do
    {
      newly_refined_elements = 0;

      // Refinement results should agree on ghost elements, except for
      // id and unique_id assignment
      for (auto & elem : _mesh.element_ptr_range())
        {
          auto it = coarse_parent.find(elem);
          const dof_id_type coarse_id =
            (it == coarse_parent.end()) ?
            elem->id() : it->second;

          newly_refined_elements +=
            this->refine_via_edges(*elem, coarse_id);
        }
      refined_elements += newly_refined_elements;
      for (auto & elem : this->new_elements)
        _mesh.add_elem(std::move(elem));
      this->new_elements.clear();

      // Yeah, this is just a lower bound in parallel
      _mesh.comm().max(newly_refined_elements);

      if (newly_refined_elements && !_mesh.is_replicated())
        {
          bool have_edge_queries = !edge_queries.empty();
          _mesh.comm().max(have_edge_queries);
          refinement_datum * refinement_data_ex = nullptr;
          if (have_edge_queries)
            Parallel::pull_parallel_vector_data
              (_mesh.comm(), edge_queries, edge_gather_functor,
               edge_action_functor, refinement_data_ex);
        }
    }
  while (newly_refined_elements);

  // If we're replicated then we should be done here.
  if (!_mesh.is_replicated())
    {
      // If we're not replicated then we might have a bunch of nodes
      // and elements that still have temporary id and unique_id
      // values.  Give them permanent ones.
      std::unordered_map<processor_id_type, std::vector<dof_id_type>> elems_to_query;
      for (const auto & [coarse_id, added_elem_map] : added_elements)
        {
          if (added_elem_map.empty())
            continue;

          const processor_id_type pid =
            added_elem_map.begin()->second->processor_id();
          if (pid == _mesh.processor_id())
            continue;

          elems_to_query[pid].push_back(coarse_id);
        }

      // Return fine element data based on vertex_average()
      typedef
#ifdef LIBMESH_ENABLE_UNIQUE_ID
        std::vector<std::tuple<Point, dof_id_type, unique_id_type>>
#else
        std::vector<std::tuple<Point, dof_id_type>>
#endif
        elem_refinement_datum;

      auto added_gather_functor =
        [this]
        (processor_id_type,
         const std::vector<dof_id_type> & coarse_elems,
         std::vector<elem_refinement_datum> & coarse_refinements)
        {
          // Fill those requests
          const std::size_t query_size = coarse_elems.size();
          coarse_refinements.resize(query_size);
          for (auto i : make_range(query_size))
            {
              const dof_id_type coarse_id = coarse_elems[i];
              const auto & added =
                libmesh_map_find(added_elements, coarse_id);

              for (auto [vertex_avg, elem] : added)
                {
                  coarse_refinements[i].emplace_back
                    (vertex_avg,
                     elem->id()
#ifdef LIBMESH_ENABLE_UNIQUE_ID
                     , elem->unique_id()
#endif
                    );
                }
            }
        };

      auto added_action_functor =
        [this]
        (processor_id_type,
         const std::vector<dof_id_type> & coarse_elems,
         const std::vector<elem_refinement_datum> & coarse_refinements)
        {
          const std::size_t query_size = coarse_elems.size();
          for (auto i : make_range(query_size))
            {
              const dof_id_type coarse_id = coarse_elems[i];
              const auto & refinement_data = coarse_refinements[i];
              const auto & our_added =
                libmesh_map_find(added_elements, coarse_id);
              for (auto [vertex_avg, id
#ifdef LIBMESH_ENABLE_UNIQUE_ID
                         , unique_id
#endif
                   ] : refinement_data)
                {
                  Elem & our_elem = *libmesh_map_find(our_added,
                                                      vertex_avg);
                  libmesh_assert_equal_to(our_elem.vertex_average(),
                                          vertex_avg);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
                  our_elem.set_unique_id(unique_id);
#endif
                  _mesh.renumber_elem(our_elem.id(), id);
                }
            }
        };

      elem_refinement_datum * refinement_data_ex = nullptr;
      Parallel::pull_parallel_vector_data
        (_mesh.comm(), elems_to_query, added_gather_functor,
         added_action_functor, refinement_data_ex);

      // That took care of our element ids; now get node ids.
      MeshCommunication mc;
      mc.make_node_proc_ids_parallel_consistent (_mesh);
      mc.make_node_ids_parallel_consistent (_mesh);
      mc.make_node_unique_ids_parallel_consistent (_mesh);
    }

  // In theory the operations on ghost elements are, after remote edge
  // splittings are known, embarrassingly parallel.  In practice ...
  // let's verify that we haven't just desynchronized our mesh.
#ifdef DEBUG
  MeshTools::libmesh_assert_equal_points(_mesh);
  MeshTools::libmesh_assert_equal_connectivity (_mesh);
#  ifdef LIBMESH_ENABLE_UNIQUE_ID
  MeshTools::libmesh_assert_valid_unique_ids(_mesh);
#  endif
#endif

  return refined_elements;
}



void SimplexRefiner::set_desired_volume_function
  (FunctionBase<Real> * desired)
{
  if (desired)
    _desired_volume_func = desired->clone();
  else
    _desired_volume_func.reset();
}



FunctionBase<Real> * SimplexRefiner::get_desired_volume_function ()
{
  return _desired_volume_func.get();
}


void SimplexRefiner::fill_refinement_datum(std::pair<Node *, Node *> vertices,
                                           refinement_datum & vec)
{
  auto it = new_nodes.find(vertices);
  if (it == new_nodes.end())
    return;

  vec.emplace_back(vertices.first->id(),
                   vertices.second->id(),
                   it->second->id(),
                   it->second->processor_id());

  std::pair<Node *, Node *> subedge {vertices.first, it->second};
  if ((Point&)(*subedge.first) >
      (Point&)(*subedge.second))
    std::swap(subedge.first, subedge.second);
  fill_refinement_datum(subedge, vec);

  subedge = std::make_pair(it->second, vertices.second);
  if ((Point&)(*subedge.first) >
      (Point&)(*subedge.second))
    std::swap(subedge.first, subedge.second);
  fill_refinement_datum(subedge, vec);
}


} // namespace libMesh

