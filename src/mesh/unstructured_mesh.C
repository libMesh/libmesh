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

#include "libmesh/unstructured_mesh_impl.h"

namespace libMesh
{
template <>
void UnstructuredMeshTempl<Real>::read (const std::string & name,
                                        void *,
                                        bool skip_renumber_nodes_and_elements,
                                        bool skip_find_neighbors)
{
  // Set the skip_renumber_nodes_and_elements flag on all processors
  // if necessary.
  // This ensures that renumber_nodes_and_elements is *not* called
  // during prepare_for_use() for certain types of mesh files.
  // This is required in cases where there is an associated solution
  // file which expects a certain ordering of the nodes.
  if (name.rfind(".gmv") + 4 == name.size())
    this->allow_renumbering(false);

  NameBasedIO(*this).read(name);

  if (skip_renumber_nodes_and_elements)
    {
      // Use MeshBase::allow_renumbering() yourself instead.
      libmesh_deprecated();
      this->allow_renumbering(false);
    }

  // Done reading the mesh.  Now prepare it for use.
  this->prepare_for_use(/*skip_renumber (deprecated)*/ false,
                        skip_find_neighbors);
}

template <>
void UnstructuredMeshTempl<Real>::write (const std::string & name)
{
  LOG_SCOPE("write()", "Mesh");

  NameBasedIO(*this).write(name);
}

template <>
void UnstructuredMeshTempl<Real>::write (const std::string & name,
                                         const std::vector<Number> & v,
                                         const std::vector<std::string> & vn)
{
  LOG_SCOPE("write()", "Mesh");

  NameBasedIO(*this).write_nodal_data(name, v, vn);
}

template <>
void UnstructuredMeshTempl<Real>::all_first_order ()
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  START_LOG("all_first_order()", "Mesh");

  /**
   * Prepare to identify (and then delete) a bunch of no-longer-used nodes.
   */
  std::vector<bool> node_touched_by_me(this->max_node_id(), false);

  // Loop over the high-ordered elements.
  // First make sure they _are_ indeed high-order, and then replace
  // them with an equivalent first-order element.
  for (auto & so_elem : this->element_ptr_range())
    {
      libmesh_assert(so_elem);

      /*
       * build the first-order equivalent, add to
       * the new_elements list.
       */
      Elem * lo_elem = Elem::build
        (Elem::first_order_equivalent_type
         (so_elem->type()), so_elem->parent()).release();

      const unsigned short n_sides = so_elem->n_sides();

      for (unsigned short s=0; s != n_sides; ++s)
        if (so_elem->neighbor_ptr(s) == RemoteElem::get_instance())
          lo_elem->set_neighbor(s, const_cast<RemoteElem *>(RemoteElem::get_instance()));

#ifdef LIBMESH_ENABLE_AMR
      /*
       * Reset the parent links of any child elements
       */
      if (so_elem->has_children())
        for (unsigned int c = 0, nc = so_elem->n_children(); c != nc; ++c)
          {
            Elem * child = so_elem->child_ptr(c);
            if (child != RemoteElem::get_instance())
              child->set_parent(lo_elem);
            lo_elem->add_child(child, c);
          }

      /*
       * Reset the child link of any parent element
       */
      if (so_elem->parent())
        {
          unsigned int c =
            so_elem->parent()->which_child_am_i(so_elem);
          lo_elem->parent()->replace_child(lo_elem, c);
        }

      /*
       * Copy as much data to the new element as makes sense
       */
      lo_elem->set_p_level(so_elem->p_level());
      lo_elem->set_refinement_flag(so_elem->refinement_flag());
      lo_elem->set_p_refinement_flag(so_elem->p_refinement_flag());
#endif

      libmesh_assert_equal_to (lo_elem->n_vertices(), so_elem->n_vertices());

      /*
       * By definition the vertices of the linear and
       * second order element are identically numbered.
       * transfer these.
       */
      for (unsigned int v=0, snv=so_elem->n_vertices(); v < snv; v++)
        {
          lo_elem->set_node(v) = so_elem->node_ptr(v);
          node_touched_by_me[lo_elem->node_id(v)] = true;
        }

      /*
       * find_neighbors relies on remote_elem neighbor links being
       * properly maintained.
       */
      for (unsigned short s=0; s != n_sides; s++)
        {
          if (so_elem->neighbor_ptr(s) == RemoteElem::get_instance())
            lo_elem->set_neighbor(s, const_cast<RemoteElem*>(RemoteElem::get_instance()));
        }

      /**
       * If the second order element had any boundary conditions they
       * should be transferred to the first-order element.  The old
       * boundary conditions will be removed from the BoundaryInfo
       * data structure by insert_elem.
       */
      this->get_boundary_info().copy_boundary_ids
        (this->get_boundary_info(), so_elem, lo_elem);

      /*
       * The new first-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the second-order element.
       */
      lo_elem->set_id(so_elem->id());
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      lo_elem->set_unique_id() = so_elem->unique_id();
#endif
      lo_elem->processor_id() = so_elem->processor_id();
      lo_elem->subdomain_id() = so_elem->subdomain_id();
      this->insert_elem(lo_elem);
    }

  // Deleting nodes does not invalidate iterators, so this is safe.
  for (const auto & node : this->node_ptr_range())
    if (!node_touched_by_me[node->id()])
      this->delete_node(node);

  // If crazy people applied boundary info to non-vertices and then
  // deleted those non-vertices, we should make sure their boundary id
  // caches are correct.
  this->get_boundary_info().regenerate_id_sets();

  STOP_LOG("all_first_order()", "Mesh");

  // On hanging nodes that used to also be second order nodes, we
  // might now have an invalid nodal processor_id()
  Partitioner::set_node_processor_ids(*this);

  // delete or renumber nodes if desired
  this->prepare_for_use();
}

template class UnstructuredMeshTempl<Real>;
}
