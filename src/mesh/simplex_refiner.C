// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/node.h"
#include "libmesh/partitioner.h"
#include "libmesh/remote_elem.h"
#include "libmesh/sync_refinement_flags.h"
#include "libmesh/unstructured_mesh.h"



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

std::size_t SimplexRefiner::refine_via_edges(Elem & elem)
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

  int side_to_refine = 0;
  Real length_to_refine = 0;
  std::pair<Point *, Point *> vertices_to_refine;
  for (int side : make_range(3))
    {
      std::pair<Point *, Point *> vertices
        {elem.node_ptr((side+1)%3),
         elem.node_ptr(side)};

      if (vertices.first > vertices.second)
        std::swap(vertices.first, vertices.second);

      const auto sidevec = *vertices.first - *vertices.second;
      Real length = sidevec.norm();
      if (length > length_to_refine &&
          (should_refine ||
           new_nodes.find(vertices) != new_nodes.end()))
        {
          side_to_refine = side;
          length_to_refine = length;
          vertices_to_refine = vertices;
        }
    }

  if (length_to_refine)
    {
      Node * midedge_node;
      if (auto it = new_nodes.find(vertices_to_refine);
          it != new_nodes.end())
        {
          midedge_node = it->second;
        }
      else
        {
          midedge_node =
            _mesh.add_point((*vertices_to_refine.first +
                             *vertices_to_refine.second)/2);
          new_nodes[vertices_to_refine] = midedge_node;
        }

      constexpr int max_subelems = 2;
      std::unique_ptr<Tri3> subelem[max_subelems];

      std::vector<boundary_id_type> bc_ids;

      // We'll switch(elem->type()) here eventually
      subelem[0] = std::make_unique<Tri3>();
      subelem[0]->set_node(0) = elem.node_ptr(side_to_refine);
      subelem[0]->set_node(1) = midedge_node;
      subelem[0]->set_node(2) = elem.node_ptr((side_to_refine+2)%3);

      subelem[1] = std::make_unique<Tri3>();
      subelem[1]->set_node(0) = elem.node_ptr((side_to_refine+2)%3);
      subelem[1]->set_node(1) = midedge_node;
      subelem[1]->set_node(2) = elem.node_ptr((side_to_refine+1)%3);

      BoundaryInfo & boundary_info = _mesh.get_boundary_info();
      boundary_info.boundary_ids(&elem, side_to_refine, bc_ids);
      boundary_info.add_side(subelem[0].get(), 0, bc_ids);
      boundary_info.add_side(subelem[1].get(), 1, bc_ids);
      boundary_info.boundary_ids(&elem, (side_to_refine+1)%3, bc_ids);
      boundary_info.add_side(subelem[1].get(), 2, bc_ids);
      boundary_info.boundary_ids(&elem, (side_to_refine+2)%3, bc_ids);
      boundary_info.add_side(subelem[0].get(), 2, bc_ids);

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
      if (!this->new_elements.empty() &&
          &elem == this->new_elements.back().get())
        this->new_elements.pop_back();
      else
        _mesh.delete_elem(&elem);

      // We refined at least ourselves
      std::size_t n_refined_elements = 1;

      for (unsigned int i=0; i != max_subelems; ++i)
        {
          this->new_elements.push_back(std::move(subelem[i]));
          n_refined_elements +=
            this->refine_via_edges(*(this->new_elements.back()));
        }

      return n_refined_elements;
    }

  // We didn't refine anything
  return 0;
}



std::size_t SimplexRefiner::refine_via_edges()
{
  this->new_nodes.clear();
  this->new_elements.clear();

  std::size_t refined_elements = 0;
  std::size_t newly_refined_elements = 0;
  do
    {
      newly_refined_elements = 0;
      for (auto & elem : _mesh.local_element_ptr_range())
        {
          newly_refined_elements +=
            this->refine_via_edges(*elem);
        }
      refined_elements += newly_refined_elements;
      for (auto & elem : this->new_elements)
        _mesh.add_elem(std::move(elem));
      this->new_elements.clear();
    }
  while (newly_refined_elements);

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


} // namespace libMesh

