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
#include "libmesh/unstructured_mesh.h"

namespace {
  using namespace libMesh;

  void flood_component (std::unordered_set<Elem *> & all_components,
                        std::unordered_set<Elem *> & current_component,
                        Elem * elem)
  {
    libmesh_assert(elem);

    if (current_component.count(elem))
      return;

    libmesh_assert(!all_components.count(elem));

    all_components.insert(elem);
    current_component.insert(elem);

    for (auto s : make_range(elem->n_sides()))
      flood_component(all_components, current_component,
                      elem->neighbor_ptr(s));
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


void MeshTetInterface::volume_to_surface_mesh()
{
  {
    std::unordered_set<Elem *> elems_to_delete;

    std::vector<std::unique_ptr<Elem>> elems_to_add;

    // Convert all faces to surface elements
    for (auto * elem : this->_mesh.element_ptr_range())
      {
        if (elem->dim() != 2)
          elems_to_delete.insert(elem);
        for (auto s : make_range(elem->n_sides()))
          {
            if (elem->neighbor_ptr(s))
              continue;

            elems_to_add.push_back(elem->build_side_ptr(s));
          }
      }

    // Add the new elements outside the loop so we don't risk
    // invalidating iterators.  Wipe their interior_parent pointers so
    // those aren't dangling later.
    for (auto & elem : elems_to_add)
    {
      elem->set_interior_parent(nullptr);
      this->_mesh.add_elem(std::move(elem));
    }

    // Remove volume and edge elements
    for (Elem * elem : elems_to_delete)
      this->_mesh.delete_elem(elem);
  }

  // Fix up neighbor pointers, element counts, etc.
  this->_mesh.prepare_for_use();

  // Partition surface into connected components
  std::vector<std::unordered_set<Elem *>> components;
  std::unordered_set<Elem *> in_component;

  for (auto * elem : this->_mesh.element_ptr_range())
    if (!in_component.count(elem))
      {
        components.push_back({});
        flood_component(in_component, components.back(), elem);
      }

  const std::unordered_set<Elem *> * biggest_component = nullptr;
  Real biggest_six_vol = 0;
  for (const auto & component : components)
    {
      Real six_vol = six_times_signed_volume(component);
      if (six_vol > biggest_six_vol)
        {
          six_vol = biggest_six_vol;
          biggest_component = &component;
        }
    }

  if (!biggest_component)
    libmesh_error_msg("No positive-volume component found among " <<
                      components.size() << " boundary components");

  for (const auto & component : components)
    if (&component != biggest_component)
      for (Elem * elem: component)
        this->_mesh.delete_elem(elem);

  this->_mesh.prepare_for_use();
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
