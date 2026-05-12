// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_config.h"


// C++ includes
#include <sstream>

// Local includes
#include "libmesh/mesh_tet_interface.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/remote_elem.h"
#include "libmesh/unstructured_mesh.h"

namespace {
  using namespace libMesh;

  std::unordered_set<Elem *>
  flood_component (std::unordered_set<Elem *> & all_components,
                        Elem * elem)
  {
    libmesh_assert(!all_components.count(elem));

    std::unordered_set<Elem *> current_component;

    std::unordered_set<Elem *> elements_to_consider = {elem};

    while (!elements_to_consider.empty())
      {
        std::unordered_set<Elem *> next_elements_to_consider;

        for (Elem * considering : elements_to_consider)
          {
            all_components.insert(considering);
            current_component.insert(considering);

            for (auto s : make_range(considering->n_sides()))
              {
                Elem * neigh = considering->neighbor_ptr(s);

                libmesh_error_msg_if
                  (!neigh,
                   "Tet generation encountered a 2D element with a null neighbor, but a\n"
                   "boundary must be a 2D closed manifold (surface).\n");

                if (all_components.find(neigh) == all_components.end())
                  next_elements_to_consider.insert(neigh);
              }
          }
        elements_to_consider = next_elements_to_consider;
      }

    return current_component;
  }

  // Returns six times the signed volume of a tet formed by the given
  // 3 points and the origin
  Real six_times_signed_tet_volume (const Point & p1,
                                    const Point & p2,
                                    const Point & p3)
  {
    return p1(0)*p2(1)*p3(2)
         - p1(0)*p3(1)*p2(2)
         - p2(0)*p1(1)*p3(2)
         + p2(0)*p3(1)*p1(2)
         + p3(0)*p1(1)*p2(2)
         - p3(0)*p2(1)*p1(2);
  }

  // Returns six times the signed volume of the space defined by the
  // manifold of surface elements in component
  Real six_times_signed_volume (const std::unordered_set<Elem *> component)
  {
    Real six_vol = 0;

    for (const Elem * elem: component)
      {
        libmesh_assert_equal_to(elem->dim(), 2);
        for (auto n : make_range(elem->n_vertices()-2))
          six_vol += six_times_signed_tet_volume(elem->point(0),
                                                 elem->point(n+1),
                                                 elem->point(n+2));
      }

    return six_vol;
  }
}

namespace libMesh
{

//----------------------------------------------------------------------
// MeshTetInterface class members
MeshTetInterface::MeshTetInterface (UnstructuredMesh & mesh) :
  _verbosity(0), _desired_volume(0), _smooth_after_generating(false),
  _elem_type(TET4), _mesh(mesh)
{
}


MeshTetInterface::~MeshTetInterface() = default;


void MeshTetInterface::attach_hole_list
  (std::unique_ptr<std::vector<std::unique_ptr<UnstructuredMesh>>> holes)
{
  _holes = std::move(holes);
}


BoundingBox MeshTetInterface::volume_to_surface_mesh(UnstructuredMesh & mesh)
{
  // If we've been handed an unprepared mesh then we need to be made
  // aware of that and fix that; we're relying on neighbor pointers.
  MeshTools::libmesh_assert_valid_is_prepared(mesh);

  if (!mesh.is_prepared())
    mesh.prepare_for_use();

  // We'll return a bounding box for use by subclasses in basic sanity checks.
  BoundingBox surface_bb;

  // First convert all volume boundaries to surface elements; this
  // gives us a manifold bounding the mesh, though it may not be a
  // connected manifold even if the volume mesh was connected.
  {
    // Make sure ids are in sync and valid on a DistributedMesh
    const dof_id_type max_orig_id = mesh.max_elem_id();
#ifdef LIBMESH_ENABLE_UNIQUE_ID
    const unique_id_type max_unique_id = mesh.parallel_max_unique_id();
#endif

    // Change this if we add arbitrary polyhedra...
    const dof_id_type max_sides = 6;

    std::unordered_set<Elem *> elems_to_delete;

    std::vector<std::unique_ptr<Elem>> elems_to_add;

    // Convert all faces to surface elements
    for (auto * elem : mesh.active_element_ptr_range())
      {
        libmesh_error_msg_if (elem->dim() < 2,
          "Cannot use meshes with 0D or 1D elements to define a volume");

        // If we've already got 2D elements then those are (part of)
        // our surface.
        if (elem->dim() == 2)
          continue;

        // 3D elements will be removed after we've extracted their
        // surface faces.
        elems_to_delete.insert(elem);

        for (auto s : make_range(elem->n_sides()))
          {
            // If there's a neighbor on this side then there's not a
            // boundary
            if (elem->neighbor_ptr(s))
              {
                // We're not supporting AMR meshes here yet
                if (elem->level() != elem->neighbor_ptr(s)->level())
                  libmesh_not_implemented_msg
                    ("Tetrahedralizaton of adapted meshes is not currently supported");
                continue;
              }

            elems_to_add.push_back(elem->build_side_ptr(s));
            Elem * side_elem = elems_to_add.back().get();

            // Wipe the interior_parent before it can become a
            // dangling pointer later
            side_elem->set_interior_parent(nullptr);

            // If the mesh is replicated then its automatic id
            // setting is fine.  If not, then we need unambiguous ids
            // independent of element traversal.
            if (!mesh.is_replicated())
              {
                side_elem->set_id(max_orig_id + max_sides*elem->id() + s);
#ifdef LIBMESH_ENABLE_UNIQUE_ID
                side_elem->set_unique_id(max_unique_id + max_sides*elem->id() + s);
#endif
              }
          }
      }

    // If the mesh is replicated then its automatic neighbor finding
    // is fine.  If not, then we need to insert them ourselves, but
    // it's easy because we can use the fact (from our implementation
    // above) that our new elements have no parents or children, plus
    // the fact (from the tiny fraction of homology I understand) that
    // a manifold boundary is a manifold with no boundary.
    //
    // See UnstructuredMesh::find_neighbors() for more explanation of
    // (a more complicated version of) the algorithm here.
    if (!mesh.is_replicated())
      {
        typedef dof_id_type                     key_type;
        typedef std::pair<Elem *, unsigned char> val_type;
        typedef std::unordered_multimap<key_type, val_type> map_type;
        map_type side_to_elem_map;

        std::unique_ptr<Elem> my_side, their_side;

        for (auto & elem : elems_to_add)
          {
            for (auto s : elem->side_index_range())
              {
                if (elem->neighbor_ptr(s))
                  continue;
                const dof_id_type key = elem->low_order_key(s);
                auto bounds = side_to_elem_map.equal_range(key);
                if (bounds.first != bounds.second)
                  {
                    elem->side_ptr(my_side, s);
                    while (bounds.first != bounds.second)
                      {
                        Elem * potential_neighbor = bounds.first->second.first;
                        const unsigned int ns = bounds.first->second.second;
                        potential_neighbor->side_ptr(their_side, ns);
                        if (*my_side == *their_side)
                          {
                            elem->set_neighbor(s, potential_neighbor);
                            potential_neighbor->set_neighbor(ns, elem.get());
                            side_to_elem_map.erase (bounds.first);
                            break;
                          }
                        ++bounds.first;
                      }

                    if (!elem->neighbor_ptr(s))
                      side_to_elem_map.emplace
                        (key, std::make_pair(elem.get(), cast_int<unsigned char>(s)));
                  }
              }
          }

        // At this point we *should* have a match for everything, so
        // anything we don't have a match for is remote.
        for (auto & elem : elems_to_add)
          for (auto s : elem->side_index_range())
            if (!elem->neighbor_ptr(s))
              elem->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));
      }

    // Remove volume and edge elements
    for (Elem * elem : elems_to_delete)
      mesh.delete_elem(elem);

    // Add the new elements outside the loop so we don't risk
    // invalidating iterators.
    for (auto & elem : elems_to_add)
      mesh.add_elem(std::move(elem));
  }

  // Fix up neighbor pointers, element counts, etc.
  mesh.prepare_for_use();

  // We're making tets; we need to start with tris
  MeshTools::Modification::all_tri(mesh);

  // Partition surface into connected components.  At this point I'm
  // finally going to give up and serialize, because at least we got
  // from 3D down to 2D first, and because I don't want to have to
  // turn flood_component into a while loop with a parallel sync in
  // the middle, and because we do have to serialize *eventually*
  // anyways unless we get a parallel tetrahedralizer backend someday.
  MeshSerializer mesh_serializer(mesh);

  std::vector<std::unordered_set<Elem *>> components;
  std::unordered_set<Elem *> in_component;

  for (auto * elem : mesh.element_ptr_range())
    if (!in_component.count(elem))
      components.emplace_back(flood_component(in_component, elem));

  const std::unordered_set<Elem *> * biggest_component = nullptr;
  Real biggest_six_vol = 0;
  for (const auto & component : components)
    {
      Real six_vol = six_times_signed_volume(component);
      if (std::abs(six_vol) > std::abs(biggest_six_vol))
        {
          biggest_six_vol = six_vol;
          biggest_component = &component;
        }
    }

  if (!biggest_component)
    libmesh_error_msg("No non-zero-volume component found among " <<
                      components.size() << " boundary components");

  for (const auto & component : components)
    if (&component != biggest_component)
      {
        for (Elem * elem: component)
          mesh.delete_elem(elem);
      }
    else
      {
        for (Elem * elem: component)
          {
            if (biggest_six_vol < 0)
              elem->flip(&mesh.get_boundary_info());

            for (auto & node : elem->node_ref_range())
              surface_bb.union_with(node);
          }
      }

  mesh.prepare_for_use();

  return surface_bb;
}


std::set<MeshTetInterface::SurfaceIntegrity> MeshTetInterface::check_hull_integrity() const
{
  // Check for easy return: if the Mesh is empty (i.e. if
  // somebody called triangulate_conformingDelaunayMesh on
  // a Mesh with no elements, then hull integrity check must
  // fail...
  if (_mesh.n_elem() == 0)
    return {EMPTY_MESH};

  std::set<MeshTetInterface::SurfaceIntegrity> returnval;

  const BoundingBox bb = MeshTools::create_bounding_box(this->_mesh);
  const Point extents = bb.max() - bb.min();
  if (extents(0) == 0 ||
      extents(1) == 0 ||
      extents(2) == 0)
    returnval.insert(DEGENERATE_MESH);

  // Figure a area to use for relative tolerances when detecting
  // degenerate elements
  const Real ref_area = std::abs(extents(0) * extents(1)) +
                        std::abs(extents(0) * extents(2)) +
                        std::abs(extents(1) * extents(2));

  for (auto & elem : this->_mesh.element_ptr_range())
    {
      // Check for proper element type
      if (elem->type() != TRI3)
        {
          if (this->_verbosity >= 50)
            std::cerr << "Non-Tri3: " << elem->get_info() << std::endl;
          returnval.insert(NON_TRI3);
        }

      // Make sure it's a decent element.
      if (elem->volume() < ref_area * TOLERANCE * TOLERANCE)
        {
          if (this->_verbosity >= 50)
            std::cerr << "Degenerate element: " << elem->get_info() << std::endl;
          returnval.insert(DEGENERATE_ELEMENT);
        }

      for (auto s : elem->side_index_range())
        {
          const Elem * const neigh = elem->neighbor_ptr(s);

          if (neigh == nullptr)
            {
              if (this->_verbosity >= 50)
                std::cerr << "Element missing neighbor " << s << ": " << elem->get_info() << std::endl;
              returnval.insert(MISSING_NEIGHBOR);
              continue;
            }

          // Make sure our neighbor points back to us
          const unsigned int nn = neigh->which_neighbor_am_i(elem);

          if (nn >= 3)
            {
              if (this->_verbosity >= 50)
                std::cerr << "Element missing backlink " << s << ": " << elem->get_info() << std::endl;
              returnval.insert(MISSING_BACKLINK);
              continue;
            }

          // Our neighbor should have the same the edge nodes we do on
          // the neighboring edgei
          const Node * const n1 = elem->node_ptr(s);
          const Node * const n2 = elem->node_ptr((s+1)%3);

          const unsigned int i1 = neigh->local_node(n1->id());
          const unsigned int i2 = neigh->local_node(n2->id());
          if (i1 >= 3 || i2 >= 3)
            {
              if (this->_verbosity >= 50)
                std::cerr << "Element with bad neighbor " << s << " nodes: " << elem->get_info() << std::endl;
              returnval.insert(BAD_NEIGHBOR_NODES);
              continue;
            }

          // It should have those edge nodes in the opposite order
          // (because they have the same orientation we do)
          if ((i2 + 1)%3 != i1)
            {
              if (this->_verbosity >= 50)
                std::cerr << "Element orientation mismatch with neighbor " << s << ": " << elem->get_info() << std::endl;
              returnval.insert(NON_ORIENTED);
              continue;
            }

          // And it should have those edge nodes in the expected
          // places relative to its neighbor link
          if (i2 != nn)
            {
              if (this->_verbosity >= 50)
                std::cerr << "Element with bad links on neighbor " << s << ": " << elem->get_info() << std::endl;
              returnval.insert(BAD_NEIGHBOR_LINKS);
              continue;
            }
        }
    }

  // Return anything and everything we found
  return returnval;
}



std::set<MeshTetInterface::SurfaceIntegrity> MeshTetInterface::improve_hull_integrity()
{
  // We don't really do anything parallel here, but we aspire to.
  libmesh_parallel_only(this->_mesh.comm());

  std::set<MeshTetInterface::SurfaceIntegrity> integrityproblems =
    this->check_hull_integrity();

  // If we have no problem, or a problem we can't fix, we're done.
  if (integrityproblems.empty() ||
      integrityproblems.count(NON_TRI3) ||
      integrityproblems.count(EMPTY_MESH))
    return integrityproblems;

  // Possibly the user gave us an unprepared mesh with missing or bad
  // neighbor links?
  if (integrityproblems.count(MISSING_NEIGHBOR) ||
      integrityproblems.count(MISSING_BACKLINK) ||
      integrityproblems.count(BAD_NEIGHBOR_LINKS))
  {
    this->_mesh.find_neighbors();
    integrityproblems = this->check_hull_integrity();
  }

  // If find_neighbors() doesn't fix these, I give up.
  if (integrityproblems.count(MISSING_NEIGHBOR) ||
      integrityproblems.count(MISSING_BACKLINK) ||
      integrityproblems.count(BAD_NEIGHBOR_LINKS))
    return integrityproblems;

  // find_neighbors() might have fixed everything
  if (integrityproblems.empty())
    return integrityproblems;

  // A non-oriented (but orientable!) surface is the only thing we
  // shouldn't have fixed or given up on by now.
  libmesh_assert_equal_to(integrityproblems.size(), 1);
  libmesh_assert_equal_to(integrityproblems.count(NON_ORIENTED), 1);

  // We need one known-good triangle to start from.  We'll pick the
  // most-negative-x normal among the triangles on the most-negative-x
  // point.

  // We'll just implement this in serial for now.
  MeshSerializer mesh_serializer(this->_mesh);

  // I don't see why we'd need boundary info here, but maybe we'll
  // want to preserve edge/node conditions eventually?
  BoundaryInfo & bi = this->_mesh.get_boundary_info();

  const Node * lowest_point = (*this->_mesh.elements_begin())->node_ptr(0);

  // Index by ids, not pointers, for consistency in parallel
  std::unordered_set<dof_id_type> attached_elements;

  for (Elem * elem : this->_mesh.element_ptr_range())
    {
      for (const Node & node : elem->node_ref_range())
        {
          if (node(0) < (*lowest_point)(0))
            {
              lowest_point = &node;
              attached_elements.clear();
            }
          if (&node == lowest_point)
            attached_elements.insert(elem->id());
        }
    }

  Elem * best_elem = nullptr;
  Real best_abs_normal_0 = 0;

  for (dof_id_type id : attached_elements)
    {
      Elem * elem = this->_mesh.elem_ptr(id);
      const Point e01 = elem->point(1) - elem->point(0);
      const Point e02 = elem->point(2) - elem->point(0);
      const Point normal = e01.cross(e02).unit();
      const Real abs_normal_0 = std::abs(normal(0));

      if (!best_elem || abs_normal_0 > best_abs_normal_0)
        {
          best_elem = elem;
          best_abs_normal_0 = abs_normal_0;

          // Make sure that element is actually a good one, by
          // flipping it if it's not.
          if (abs_normal_0 == normal(0))
            elem->flip(&bi);
        }
    }

  // Now flood-fill from that element to get a consistent orientation
  // for the others.
  std::unordered_set<dof_id_type> frontier_elements{best_elem->id()},
                                  finished_elements{};

  while (!frontier_elements.empty())
    {
      const dof_id_type elem_id = *frontier_elements.begin();
      Elem & elem = this->_mesh.elem_ref(elem_id);
      for (auto s : elem.side_index_range())
        {
          Elem * neigh = elem.neighbor_ptr(s);
          libmesh_assert(neigh);
          libmesh_assert_less(neigh->which_neighbor_am_i(&elem), 3);

          const Node * const n1 = elem.node_ptr(s);
          const Node * const n2 = elem.node_ptr((s+1)%3);
          const unsigned int i1 = neigh->local_node(n1->id());
          const unsigned int i2 = neigh->local_node(n2->id());
          libmesh_assert_less(i1, 3);
          libmesh_assert_less(i2, 3);

          const dof_id_type neigh_id = neigh->id();

          const bool frontier_neigh = frontier_elements.count(neigh_id);
          const bool finished_neigh = finished_elements.count(neigh_id);

          // Are we flipped?
          if ((i2 + 1)%3 != i1)
            {
              // Are we a Moebius strip???  We give up.
              if (frontier_neigh || finished_neigh)
                return integrityproblems;

              neigh->flip(&bi);
            }

          if (!frontier_neigh && !finished_neigh)
            frontier_elements.insert(neigh_id);
        }

      finished_elements.insert(elem_id);
      frontier_elements.erase(elem_id);
    }

  this->_mesh.find_neighbors();

  libmesh_assert(this->check_hull_integrity().empty());

  return {};
}


void MeshTetInterface::process_hull_integrity_result
  (const std::set<MeshTetInterface::SurfaceIntegrity> & result) const
{
  std::ostringstream err_msg;

  if (result.empty()) // success
    return;

  err_msg << "Error! Conforming Delaunay mesh tetrahedralization requires a convex hull." << std::endl;

  if (result.count(NON_TRI3))
    {
      err_msg << "At least one non-Tri3 element was found in the input boundary mesh.  ";
      err_msg << "Our constrained Delaunay tetrahedralization boundary must be a triangulation of Tri3 elements." << std::endl;
    }
  if (result.count(MISSING_NEIGHBOR))
    {
      err_msg << "At least one triangle without three neighbors was found in the input boundary mesh.  ";
      err_msg << "A constrained Delaunay tetrahedralization boundary must be a triangular manifold without boundary." << std::endl;
    }
  if (result.count(EMPTY_MESH))
    {
      err_msg << "The input boundary mesh was empty!" << std::endl;
      err_msg << "Our constrained Delaunay tetrahedralization boundary must be a triangulation of Tri3 elements." << std::endl;
    }
  if (result.count(MISSING_BACKLINK))
    {
      err_msg << "At least one triangle neighbor without a return neighbor link was found in the input boundary mesh.  ";
      err_msg << "A constrained Delaunay tetrahedralization boundary must be a conforming and non-adaptively-refined mesh." << std::endl;
    }
  if (result.count(BAD_NEIGHBOR_NODES))
    {
      err_msg << "At least one triangle neighbor without expected node links was found in the input boundary mesh.  ";
      err_msg << "A constrained Delaunay tetrahedralization boundary must be a conforming and non-adaptively-refined mesh." << std::endl;
    }
  if (result.count(NON_ORIENTED))
    {
      err_msg << "At least one triangle neighbor with an inconsistent orientation was found in the input boundary mesh.  ";
      err_msg << "A constrained Delaunay tetrahedralization boundary must be an oriented Tri3 mesh." << std::endl;
    }
  if (result.count(BAD_NEIGHBOR_LINKS))
    err_msg << "At least one triangle neighbor with inconsistent node and neighbor links was found in the input boundary mesh." << std::endl;
  if (result.count(DEGENERATE_ELEMENT))
    err_msg << "At least one input triangle is degenerate, with near-zero area relative to the manifold." << std::endl;
  if (result.count(DEGENERATE_MESH))
    err_msg << "The input mesh is degenerate, with zero thickness in at least one direction." << std::endl;

  libmesh_error_msg(err_msg.str());
}



void MeshTetInterface::delete_2D_hull_elements()
{
  for (auto & elem : this->_mesh.element_ptr_range())
    {
      // Check for proper element type. Yes, we legally delete elements while
      // iterating over them because no entries from the underlying container
      // are actually erased.
      if (elem->type() == TRI3)
        _mesh.delete_elem(elem);
    }

  // We just removed any boundary info associated with hull element
  // edges, so let's update the boundary id caches.
  this->_mesh.get_boundary_info().regenerate_id_sets();
}



} // namespace libMesh
