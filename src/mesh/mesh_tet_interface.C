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
  _desired_volume(0), _smooth_after_generating(false),
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


unsigned MeshTetInterface::check_hull_integrity()
{
  // Check for easy return: if the Mesh is empty (i.e. if
  // somebody called triangulate_conformingDelaunayMesh on
  // a Mesh with no elements, then hull integrity check must
  // fail...
  if (_mesh.n_elem() == 0)
    return 3;

  for (auto & elem : this->_mesh.element_ptr_range())
    {
      // Check for proper element type
      if (elem->type() != TRI3)
        {
          //libmesh_error_msg("ERROR: Some of the elements in the original mesh were not TRI3!");
          return 1;
        }

      for (auto neigh : elem->neighbor_ptr_range())
        {
          if (neigh == nullptr)
            {
              // libmesh_error_msg("ERROR: Non-convex hull, cannot be tetrahedralized.");
              return 2;
            }
        }
    }

  // If we made it here, return success!
  return 0;
}



void MeshTetInterface::process_hull_integrity_result(unsigned result)
{
  if (result != 0)
    {
      libMesh::err << "Error! Conforming Delaunay mesh tetrahedralization requires a convex hull." << std::endl;

      if (result==1)
        {
          libMesh::err << "Non-TRI3 elements were found in the input Mesh.  ";
          libMesh::err << "A constrained Delaunay triangulation requires a convex hull of TRI3 elements." << std::endl;
        }

      libmesh_error_msg("Consider calling TetGenMeshInterface::pointset_convexhull() followed by Mesh::find_neighbors() first.");
    }
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
