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
#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/node_range.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/sphere.h"
#include "libmesh/threads.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/int_range.h"
#include "libmesh/utility.h"
#include "libmesh/boundary_info.h"

#ifndef NDEBUG
#  include "libmesh/remote_elem.h"
#endif

// C++ includes
#include <limits>
#include <numeric> // for std::accumulate
#include <set>
#include <unordered_map>
#include <unordered_set>



// ------------------------------------------------------------
// anonymous namespace for helper classes
namespace {

using namespace libMesh;

/**
 * SumElemWeight(Range) sums the number of nodes per element
 * for each element in the provided range. The join() method
 * defines how to combine the reduction operation from two
 * distinct instances of this class which may be executed on
 * separate threads.
 */
class SumElemWeight
{
public:
  SumElemWeight () :
    _weight(0)
  {}

  SumElemWeight (SumElemWeight &, Threads::split) :
    _weight(0)
  {}

  void operator()(const ConstElemRange & range)
  {
    for (const auto & elem : range)
      _weight += elem->n_nodes();
  }

  dof_id_type weight() const
  { return _weight; }

  // If we don't have threads we never need a join, and icpc yells a
  // warning if it sees an anonymous function that's never used
#if LIBMESH_USING_THREADS
  void join (const SumElemWeight & other)
  { _weight += other.weight(); }
#endif

private:
  dof_id_type _weight;
};


/**
 * FindBBox(Range) computes the bounding box for the objects
 * in the specified range.  This class may be split and subranges
 * can be executed on separate threads.  The join() method
 * defines how the results from two separate threads are combined.
 */
class FindBBox
{
public:
  FindBBox () : _bbox()
  {}

  FindBBox (FindBBox & other, Threads::split) :
    _bbox(other._bbox)
  {}

  void operator()(const ConstNodeRange & range)
  {
    for (const auto & node : range)
      {
        libmesh_assert(node);
        _bbox.union_with(*node);
      }
  }

  void operator()(const ConstElemRange & range)
  {
    for (const auto & elem : range)
      {
        libmesh_assert(elem);
        _bbox.union_with(elem->loose_bounding_box());
      }
  }

  Point & min() { return _bbox.min(); }

  Point & max() { return _bbox.max(); }

  // If we don't have threads we never need a join, and icpc yells a
  // warning if it sees an anonymous function that's never used
#if LIBMESH_USING_THREADS
  void join (const FindBBox & other)
  {
    _bbox.union_with(other._bbox);
  }
#endif

  libMesh::BoundingBox & bbox ()
  {
    return _bbox;
  }

private:
  BoundingBox _bbox;
};

#ifdef DEBUG
void assert_semiverify_dofobj(const Parallel::Communicator & communicator,
                              const DofObject * d,
                              unsigned int sysnum = libMesh::invalid_uint)
{
  if (d)
    {
      const unsigned int n_sys = d->n_systems();

      std::vector<unsigned int> n_vars (n_sys, 0);
      for (unsigned int s = 0; s != n_sys; ++s)
        if (sysnum == s ||
            sysnum == libMesh::invalid_uint)
          n_vars[s] = d->n_vars(s);

      const unsigned int tot_n_vars =
        std::accumulate(n_vars.begin(), n_vars.end(), 0);

      std::vector<unsigned int> n_comp (tot_n_vars, 0);
      std::vector<dof_id_type> first_dof (tot_n_vars, 0);

      for (unsigned int s = 0, i=0; s != n_sys; ++s)
        {
          if (sysnum != s &&
              sysnum != libMesh::invalid_uint)
            continue;

          for (unsigned int v = 0; v != n_vars[s]; ++v, ++i)
            {
              n_comp[i] = d->n_comp(s,v);
              first_dof[i] = n_comp[i] ? d->dof_number(s,v,0) : DofObject::invalid_id;
            }
        }

      libmesh_assert(communicator.semiverify(&n_sys));
      libmesh_assert(communicator.semiverify(&n_vars));
      libmesh_assert(communicator.semiverify(&n_comp));
      libmesh_assert(communicator.semiverify(&first_dof));
    }
  else
    {
      const unsigned int * p_ui = nullptr;
      const std::vector<unsigned int> * p_vui = nullptr;
      const std::vector<dof_id_type> * p_vdid = nullptr;

      libmesh_assert(communicator.semiverify(p_ui));
      libmesh_assert(communicator.semiverify(p_vui));
      libmesh_assert(communicator.semiverify(p_vui));
      libmesh_assert(communicator.semiverify(p_vdid));
    }
}



#ifdef LIBMESH_ENABLE_UNIQUE_ID
void assert_dofobj_unique_id(const Parallel::Communicator & comm,
                             const DofObject * d,
                             const std::unordered_set<unique_id_type> & unique_ids)
{
  // Duplicating some semiverify code here so we can reuse
  // tempmin,tempmax afterward

  unique_id_type tempmin, tempmax;
  if (d)
    {
      tempmin = tempmax = d->unique_id();
    }
  else
    {
      TIMPI::Attributes<unique_id_type>::set_highest(tempmin);
      TIMPI::Attributes<unique_id_type>::set_lowest(tempmax);
    }
  comm.min(tempmin);
  comm.max(tempmax);
  bool invalid = d && ((d->unique_id() != tempmin) ||
                       (d->unique_id() != tempmax));
  comm.max(invalid);

  // First verify that everything is in sync
  libmesh_assert(!invalid);

  // Then verify that any remote id doesn't duplicate a local one.
  if (!d && tempmin == tempmax)
    libmesh_assert(!unique_ids.count(tempmin));
}
#endif // LIBMESH_ENABLE_UNIQUE_ID
#endif // DEBUG

void find_nodal_neighbors_helper(const dof_id_type global_id,
                                 const std::vector<const Elem *> & node_to_elem_vec,
                                 std::vector<const Node *> & neighbors)
{
  // We'll construct a std::set<const Node *> for more efficient
  // searching while finding the nodal neighbors, and return it to the
  // user in a std::vector.
  std::set<const Node *> neighbor_set;

  // Look through the elements that contain this node
  // find the local node id... then find the side that
  // node lives on in the element
  // next, look for the _other_ node on that side
  // That other node is a "nodal_neighbor"... save it
  for (const auto & elem : node_to_elem_vec)
    {
      // We only care about active elements...
      if (elem->active())
        {
          // Which local node number is global_id?
          unsigned local_node_number = elem->local_node(global_id);

          // Make sure it was found
          libmesh_assert_not_equal_to(local_node_number, libMesh::invalid_uint);

          const unsigned short n_edges = elem->n_edges();

          // If this element has no edges, the edge-based algorithm below doesn't make sense.
          if (!n_edges)
            {
              switch (elem->type())
                {
                case EDGE2:
                  {
                    switch (local_node_number)
                      {
                      case 0:
                        // The other node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(1));
                        break;

                      case 1:
                        // The other node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(0));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                case EDGE3:
                  {
                    switch (local_node_number)
                      {
                        // The outside nodes have node 2 as a neighbor
                      case 0:
                      case 1:
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                        // The middle node has the outer nodes as neighbors
                      case 2:
                        neighbor_set.insert(elem->node_ptr(0));
                        neighbor_set.insert(elem->node_ptr(1));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                case EDGE4:
                  {
                    switch (local_node_number)
                      {
                      case 0:
                        // The left-middle node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                      case 1:
                        // The right-middle node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(3));
                        break;

                        // The left-middle node
                      case 2:
                        neighbor_set.insert(elem->node_ptr(0));
                        neighbor_set.insert(elem->node_ptr(3));
                        break;

                        // The right-middle node
                      case 3:
                        neighbor_set.insert(elem->node_ptr(1));
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                default:
                  libmesh_error_msg("Unrecognized ElemType: " << Utility::enum_to_string(elem->type()) << std::endl);
                }
            }

          const auto elem_order = Elem::type_to_default_order_map[elem->type()];

          // Index of the current edge
          unsigned current_edge = 0;

          const unsigned short n_nodes = elem->n_nodes();

          while (current_edge < n_edges)
            {
              // Find the edge the node is on
              bool found_edge = false;
              for (; current_edge<n_edges; ++current_edge)
                if (elem->is_node_on_edge(local_node_number, current_edge))
                  {
                    found_edge = true;
                    break;
                  }

              // Did we find one?
              if (found_edge)
                {
                  const Node * node_to_save = nullptr;

                  // Find another node in this element on this edge
                  for (unsigned other_node_this_edge = 0; other_node_this_edge != n_nodes; other_node_this_edge++)
                    {
                      const bool both_vertices = elem->is_vertex(local_node_number) &&
                                                 elem->is_vertex(other_node_this_edge);
                      if ( elem->is_node_on_edge(other_node_this_edge, current_edge) && // On the current edge
                           elem->node_id(other_node_this_edge) != global_id          && // But not the original node
                            // vertex nodes on the same edge of higher order elements are not nodal neighbors
                          (elem_order == 1 || !both_vertices))
                        {
                          // We've found a nodal neighbor!  Save a pointer to it..
                          node_to_save = elem->node_ptr(other_node_this_edge);

                          // Make sure we found something
                          libmesh_assert(node_to_save != nullptr);

                          neighbor_set.insert(node_to_save);
                        }
                    }
                }

              // Keep looking for edges, node may be on more than one edge
              current_edge++;
            }
        } // if (elem->active())
    } // for

  // Assign the entries from the set to the vector.  Note: this
  // replaces any existing contents in neighbors and modifies its size
  // accordingly.
  neighbors.assign(neighbor_set.begin(), neighbor_set.end());
}

}


namespace libMesh
{

// ------------------------------------------------------------
// MeshTools functions

namespace MeshTools
{

dof_id_type total_weight(const MeshBase & mesh)
{
  if (!mesh.is_serial())
    {
      libmesh_parallel_only(mesh.comm());
      dof_id_type weight = MeshTools::weight (mesh, mesh.processor_id());
      mesh.comm().sum(weight);
      dof_id_type unpartitioned_weight =
        MeshTools::weight (mesh, DofObject::invalid_processor_id);
      return weight + unpartitioned_weight;
    }

  SumElemWeight sew;

  Threads::parallel_reduce (ConstElemRange (mesh.elements_begin(),
                                            mesh.elements_end()),
                            sew);
  return sew.weight();

}



dof_id_type weight(const MeshBase & mesh, const processor_id_type pid)
{
  SumElemWeight sew;

  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(pid),
                                            mesh.pid_elements_end(pid)),
                            sew);
  return sew.weight();
}



void build_nodes_to_elem_map (const MeshBase & mesh,
                              std::vector<std::vector<dof_id_type>> & nodes_to_elem_map)
{
  // A vector indexed over all nodes is too inefficient to use for a
  // distributed mesh.  Use the unordered_map API instead.
  if (!mesh.is_serial())
    libmesh_deprecated();

  nodes_to_elem_map.resize (mesh.max_node_id());

  for (const auto & elem : mesh.element_ptr_range())
    for (auto & node : elem->node_ref_range())
      {
        libmesh_assert_less (node.id(), nodes_to_elem_map.size());
        libmesh_assert_less (elem->id(), mesh.n_elem());

        nodes_to_elem_map[node.id()].push_back(elem->id());
      }
}



void build_nodes_to_elem_map (const MeshBase & mesh,
                              std::vector<std::vector<const Elem *>> & nodes_to_elem_map)
{
  // A vector indexed over all nodes is too inefficient to use for a
  // distributed mesh.  Use the unordered_map API instead.
  if (!mesh.is_serial())
    libmesh_deprecated();

  nodes_to_elem_map.resize (mesh.max_node_id());

  for (const auto & elem : mesh.element_ptr_range())
    for (auto & node : elem->node_ref_range())
      {
        libmesh_assert_less (node.id(), nodes_to_elem_map.size());

        nodes_to_elem_map[node.id()].push_back(elem);
      }
}



void build_nodes_to_elem_map (const MeshBase & mesh,
                              std::unordered_map<dof_id_type, std::vector<dof_id_type>> & nodes_to_elem_map)
{
  nodes_to_elem_map.clear();

  for (const auto & elem : mesh.element_ptr_range())
    for (auto & node : elem->node_ref_range())
      nodes_to_elem_map[node.id()].push_back(elem->id());
}



void build_nodes_to_elem_map (const MeshBase & mesh,
                              std::unordered_map<dof_id_type, std::vector<const Elem *>> & nodes_to_elem_map)
{
  nodes_to_elem_map.clear();

  for (const auto & elem : mesh.element_ptr_range())
    for (auto & node : elem->node_ref_range())
      nodes_to_elem_map[node.id()].push_back(elem);
}



std::unordered_set<dof_id_type>
find_boundary_nodes(const MeshBase & mesh)
{
  std::unordered_set<dof_id_type> boundary_nodes;

  // Loop over elements, find those on boundary, and
  // mark them as true in on_boundary.
  for (const auto & elem : mesh.active_element_ptr_range())
    for (auto s : elem->side_index_range())
      if (elem->neighbor_ptr(s) == nullptr) // on the boundary
        {
          auto nodes_on_side = elem->nodes_on_side(s);

          for (auto & local_id : nodes_on_side)
            boundary_nodes.insert(elem->node_ptr(local_id)->id());
        }

  return boundary_nodes;
}

std::unordered_set<dof_id_type>
find_block_boundary_nodes(const MeshBase & mesh)
{
  std::unordered_set<dof_id_type> block_boundary_nodes;

  // Loop over elements, find those on boundary, and
  // mark them as true in on_boundary.
  for (const auto & elem : mesh.active_element_ptr_range())
    for (auto s : elem->side_index_range())
      if (elem->neighbor_ptr(s) && (elem->neighbor_ptr(s)->subdomain_id() != elem->subdomain_id()))
        {
          auto nodes_on_side = elem->nodes_on_side(s);

          for (auto & local_id : nodes_on_side)
            block_boundary_nodes.insert(elem->node_ptr(local_id)->id());
        }

  return block_boundary_nodes;
}



libMesh::BoundingBox
create_bounding_box (const MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  FindBBox find_bbox;

  // Start with any unpartitioned elements we know about locally
  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(DofObject::invalid_processor_id),
                                            mesh.pid_elements_end(DofObject::invalid_processor_id)),
                            find_bbox);

  // And combine with our local elements
  find_bbox.bbox().union_with(create_local_bounding_box(mesh));

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



libMesh::BoundingBox
create_nodal_bounding_box (const MeshBase & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  FindBBox find_bbox;

  // Start with any unpartitioned nodes we know about locally
  Threads::parallel_reduce (ConstNodeRange (mesh.pid_nodes_begin(DofObject::invalid_processor_id),
                                            mesh.pid_nodes_end(DofObject::invalid_processor_id)),
                            find_bbox);

  // Add our local nodes
  Threads::parallel_reduce (ConstNodeRange (mesh.local_nodes_begin(),
                                            mesh.local_nodes_end()),
                            find_bbox);

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



Sphere
bounding_sphere(const MeshBase & mesh)
{
  libMesh::BoundingBox bbox = create_bounding_box(mesh);

  const Real  diag = (bbox.second - bbox.first).norm();
  const Point cent = (bbox.second + bbox.first)/2;

  return Sphere (cent, .5*diag);
}



libMesh::BoundingBox
create_local_bounding_box (const MeshBase & mesh)
{
  FindBBox find_bbox;

  Threads::parallel_reduce (ConstElemRange (mesh.local_elements_begin(),
                                            mesh.local_elements_end()),
                            find_bbox);

  return find_bbox.bbox();
}



libMesh::BoundingBox
create_processor_bounding_box (const MeshBase & mesh,
                                          const processor_id_type pid)
{
  // This can only be run in parallel, with consistent arguments.
  libmesh_parallel_only(mesh.comm());
  libmesh_assert(mesh.comm().verify(pid));

  libmesh_assert_less (pid, mesh.n_processors());

  FindBBox find_bbox;

  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(pid),
                                            mesh.pid_elements_end(pid)),
                            find_bbox);

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



Sphere
processor_bounding_sphere (const MeshBase & mesh,
                                      const processor_id_type pid)
{
  libMesh::BoundingBox bbox =
    create_processor_bounding_box(mesh, pid);

  const Real  diag = (bbox.second - bbox.first).norm();
  const Point cent = (bbox.second + bbox.first)/2;

  return Sphere (cent, .5*diag);
}



libMesh::BoundingBox
create_subdomain_bounding_box (const MeshBase & mesh,
                               const subdomain_id_type sid)
{
  // This can only be run in parallel, with consistent arguments.
  libmesh_parallel_only(mesh.comm());
  libmesh_assert(mesh.comm().verify(sid));

  FindBBox find_bbox;

  Threads::parallel_reduce
    (ConstElemRange (mesh.active_local_subdomain_elements_begin(sid),
                     mesh.active_local_subdomain_elements_end(sid)),
     find_bbox);

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



Sphere
subdomain_bounding_sphere (const MeshBase & mesh,
                           const subdomain_id_type sid)
{
  libMesh::BoundingBox bbox =
    create_subdomain_bounding_box(mesh, sid);

  const Real  diag = (bbox.second - bbox.first).norm();
  const Point cent = (bbox.second + bbox.first)/2;

  return Sphere (cent, .5*diag);
}



void elem_types (const MeshBase & mesh,
                 std::vector<ElemType> & et)
{
  // Loop over the the elements.  If the current element type isn't in
  // the vector, insert it.
  for (const auto & elem : mesh.element_ptr_range())
    if (!std::count(et.begin(), et.end(), elem->type()))
      et.push_back(elem->type());
}



dof_id_type n_elem_of_type (const MeshBase & mesh,
                            const ElemType type)
{
  return static_cast<dof_id_type>(std::distance(mesh.type_elements_begin(type),
                                                mesh.type_elements_end  (type)));
}



dof_id_type n_active_elem_of_type (const MeshBase & mesh,
                                   const ElemType type)
{
  return static_cast<dof_id_type>(std::distance(mesh.active_type_elements_begin(type),
                                                mesh.active_type_elements_end  (type)));
}

dof_id_type n_non_subactive_elem_of_type_at_level(const MeshBase & mesh,
                                                  const ElemType type,
                                                  const unsigned int level)
{
  dof_id_type cnt = 0;

  // iterate over the elements of the specified type
  for (const auto & elem : as_range(mesh.type_elements_begin(type),
                                    mesh.type_elements_end(type)))
    if ((elem->level() == level) && !elem->subactive())
      cnt++;

  return cnt;
}


unsigned int n_active_local_levels(const MeshBase & mesh)
{
  unsigned int nl = 0;

  for (auto & elem : mesh.active_local_element_ptr_range())
    nl = std::max(elem->level() + 1, nl);

  return nl;
}



unsigned int n_active_levels(const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int nl = n_active_local_levels(mesh);

  for (const auto & elem : as_range(mesh.unpartitioned_elements_begin(),
                                    mesh.unpartitioned_elements_end()))
    if (elem->active())
      nl = std::max(elem->level() + 1, nl);

  mesh.comm().max(nl);
  return nl;
}



unsigned int n_local_levels(const MeshBase & mesh)
{
  unsigned int nl = 0;

  for (const auto & elem : as_range(mesh.local_elements_begin(),
                                    mesh.local_elements_end()))
    nl = std::max(elem->level() + 1, nl);

  return nl;
}



unsigned int n_levels(const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int nl = n_local_levels(mesh);

  for (const auto & elem : as_range(mesh.unpartitioned_elements_begin(),
                                    mesh.unpartitioned_elements_end()))
    nl = std::max(elem->level() + 1, nl);

  mesh.comm().max(nl);

  // n_levels() is only valid and should only be called in cases where
  // the mesh is validly distributed (or serialized).  Let's run an
  // expensive test in debug mode to make sure this is such a case.
#ifdef DEBUG
  const unsigned int paranoid_nl = paranoid_n_levels(mesh);
  libmesh_assert_equal_to(nl, paranoid_nl);
#endif
  return nl;
}



unsigned int paranoid_n_levels(const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int nl = 0;
  for (const auto & elem : mesh.element_ptr_range())
    nl = std::max(elem->level() + 1, nl);

  mesh.comm().max(nl);
  return nl;
}



dof_id_type n_connected_components(const MeshBase & mesh,
                                   Real constraint_tol)
{
  LOG_SCOPE("n_connected_components()", "MeshTools");

  // Yes, I'm being lazy.  This is for mesh analysis before a
  // simulation, not anything going in any loops.
  if (!mesh.is_serial_on_zero())
    libmesh_not_implemented();

  dof_id_type n_components = 0;

  if (mesh.processor_id())
  {
    mesh.comm().broadcast(n_components);
    return n_components;
  }

  // All nodes in a set here are connected (at least indirectly) to
  // all other nodes in the same set, but have not yet been discovered
  // to be connected to nodes in other sets.
  std::vector<std::unordered_set<const Node *>> components;

  // With a typical mesh with few components and somewhat-contiguous
  // ordering, vector performance should be fine.  With a mesh with
  // many components or completely scrambled ordering, performance
  // can be a disaster.
  auto find_component = [&components](const Node * n) {
    std::unordered_set<const Node *> * component = nullptr;

    for (auto & c: components)
      if (c.find(n) != c.end())
        {
          libmesh_assert(component == nullptr);
          component = &c;
        }

    return component;
  };

  auto add_to_component =
    [&find_component]
    (std::unordered_set<const Node *> & component, const Node * n) {

    auto current_component = find_component(n);
    // We may already know we're in the desired component
    if (&component == current_component)
      return;

    // If we're unknown, we should be in the desired component
    else if (!current_component)
      component.insert(n);

    // If we think we're in another component, it should actually be
    // part of the desired component
    else
      {
        component.merge(*current_component);
        libmesh_assert(current_component->empty());
      }
  };

  auto & constraint_rows = mesh.get_constraint_rows();

  for (const auto & elem : mesh.element_ptr_range())
    {
      const Node * first_node = elem->node_ptr(0);

      auto component = find_component(first_node);

      // If we didn't find one, make a new one, reusing an existing
      // slot if possible or growing our vector if necessary
      if (!component)
        for (auto & c: components)
          if (c.empty())
            component = &c;

      if (!component)
        component = &components.emplace_back();

      for (const Node & n : elem->node_ref_range())
        {
          add_to_component(*component, &n);

          auto it = constraint_rows.find(&n);
          if (it == constraint_rows.end())
            continue;

          for (const auto & [pr, val] : it->second)
            {
              // Ignore too-trivial constraint coefficients if
              // we get a non-default-0 constraint_tol
              if (std::abs(val) < constraint_tol)
                continue;

              const Elem * spline_elem = pr.first;
              libmesh_assert(spline_elem == mesh.elem_ptr(spline_elem->id()));

              const Node * spline_node =
                spline_elem->node_ptr(pr.second);

              add_to_component(*component, spline_node);
            }
        }
    }

  for (auto & component : components)
    if (!component.empty())
      ++n_components;

  // We calculated this on proc 0; now let everyone else know too
  mesh.comm().broadcast(n_components);

  return n_components;
}



void get_not_subactive_node_ids(const MeshBase & mesh,
                                std::set<dof_id_type> & not_subactive_node_ids)
{
  for (const auto & elem : mesh.element_ptr_range())
    if (!elem->subactive())
      for (auto & n : elem->node_ref_range())
        not_subactive_node_ids.insert(n.id());
}



dof_id_type n_elem (const MeshBase::const_element_iterator & begin,
                    const MeshBase::const_element_iterator & end)
{
  return cast_int<dof_id_type>(std::distance(begin, end));
}



dof_id_type n_nodes (const MeshBase::const_node_iterator & begin,
                     const MeshBase::const_node_iterator & end)
{
  return cast_int<dof_id_type>(std::distance(begin, end));
}



Real volume (const MeshBase & mesh,
             unsigned int dim)
{
  libmesh_parallel_only(mesh.comm());

  if (dim == libMesh::invalid_uint)
    dim = mesh.mesh_dimension();

  Real vol = 0;

  // first my local elements
  for (const auto & elem : as_range(mesh.local_elements_begin(),
                                    mesh.local_elements_end()))
    if (elem->dim() == dim)
      vol += elem->volume();

  // then count any unpartitioned objects, once
  if (mesh.processor_id() == 0)
    for (const auto & elem : as_range(mesh.unpartitioned_elements_begin(),
                                      mesh.unpartitioned_elements_end()))
      if (elem->dim() == dim)
        vol += elem->volume();

  mesh.comm().sum(vol);
  return vol;
}



unsigned int n_p_levels (const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int max_p_level = 0;

  // first my local elements
  for (const auto & elem : as_range(mesh.local_elements_begin(),
                                    mesh.local_elements_end()))
    max_p_level = std::max(elem->p_level(), max_p_level);

  // then any unpartitioned objects
  for (const auto & elem : as_range(mesh.unpartitioned_elements_begin(),
                                    mesh.unpartitioned_elements_end()))
    max_p_level = std::max(elem->p_level(), max_p_level);

  mesh.comm().max(max_p_level);
  return max_p_level + 1;
}



void find_nodal_neighbors(const MeshBase &,
                          const Node & node,
                          const std::vector<std::vector<const Elem *>> & nodes_to_elem_map,
                          std::vector<const Node *> & neighbors)
{
  find_nodal_neighbors_helper(node.id(), nodes_to_elem_map[node.id()],
                              neighbors);
}



void find_nodal_neighbors(const MeshBase &,
                          const Node & node,
                          const std::unordered_map<dof_id_type, std::vector<const Elem *>> & nodes_to_elem_map,
                          std::vector<const Node *> & neighbors)
{
  const std::vector<const Elem *> node_to_elem_vec =
    libmesh_map_find(nodes_to_elem_map, node.id());
  find_nodal_neighbors_helper(node.id(), node_to_elem_vec, neighbors);
}

void find_nodal_or_face_neighbors(
    const MeshBase & mesh,
    const Node & node,
    const std::unordered_map<dof_id_type, std::vector<const Elem *>> & nodes_to_elem_map,
    std::vector<const Node *> & neighbors)
{
  // Find all the nodal neighbors... that is the nodes directly connected
  // to this node through one edge.
  find_nodal_neighbors(mesh, node, nodes_to_elem_map, neighbors);

  // If no neighbors are found, use all nodes on the containing side as
  // neighbors.
  if (!neighbors.size())
    {
      // Grab the element containing node
      const auto * elem = libmesh_map_find(nodes_to_elem_map, node.id()).front();
      // Find the element side containing node
      for (const auto &side : elem->side_index_range())
        {
          const auto &nodes_on_side = elem->nodes_on_side(side);
          const auto it =
              std::find_if(nodes_on_side.begin(), nodes_on_side.end(), [&](auto local_node_id) {
                return elem->node_id(local_node_id) == node.id();
              });

          if (it != nodes_on_side.end())
            {
              for (const auto &local_node_id : nodes_on_side)
                // No need to add node itself as a neighbor
                if (const auto *node_ptr = elem->node_ptr(local_node_id);
                    *node_ptr != node)
                  neighbors.push_back(node_ptr);
              break;
            }
        }
    }
  libmesh_assert(neighbors.size());
}



void find_hanging_nodes_and_parents(const MeshBase & mesh,
                                    std::map<dof_id_type, std::vector<dof_id_type>> & hanging_nodes)
{
  // Loop through all the elements
  for (auto & elem : mesh.active_local_element_ptr_range())
    if (elem->type() == QUAD4)
      for (auto s : elem->side_index_range())
        {
          // Loop over the sides looking for sides that have hanging nodes
          // This code is inspired by compute_proj_constraints()
          const Elem * neigh = elem->neighbor_ptr(s);

          // If not a boundary side
          if (neigh != nullptr)
            {
              // Is there a coarser element next to this one?
              if (neigh->level() < elem->level())
                {
                  const Elem * ancestor = elem;
                  while (neigh->level() < ancestor->level())
                    ancestor = ancestor->parent();
                  unsigned int s_neigh = neigh->which_neighbor_am_i(ancestor);
                  libmesh_assert_less (s_neigh, neigh->n_neighbors());

                  // Couple of helper uints...
                  unsigned int local_node1=0;
                  unsigned int local_node2=0;

                  bool found_in_neighbor = false;

                  // Find the two vertices that make up this side
                  while (!elem->is_node_on_side(local_node1++,s)) { }
                  local_node1--;

                  // Start looking for the second one with the next node
                  local_node2=local_node1+1;

                  // Find the other one
                  while (!elem->is_node_on_side(local_node2++,s)) { }
                  local_node2--;

                  //Pull out their global ids:
                  dof_id_type node1 = elem->node_id(local_node1);
                  dof_id_type node2 = elem->node_id(local_node2);

                  // Now find which node is present in the neighbor
                  // FIXME This assumes a level one rule!
                  // The _other_ one is the hanging node

                  // First look for the first one
                  // FIXME could be streamlined a bit
                  for (unsigned int n=0;n<neigh->n_sides();n++)
                    if (neigh->node_id(n) == node1)
                      found_in_neighbor=true;

                  dof_id_type hanging_node=0;

                  if (!found_in_neighbor)
                    hanging_node=node1;
                  else // If it wasn't node1 then it must be node2!
                    hanging_node=node2;

                  // Reset for reuse
                  local_node1=0;

                  // Find the first node that makes up the side in the neighbor (these should be the parent nodes)
                  while (!neigh->is_node_on_side(local_node1++,s_neigh)) { }
                  local_node1--;

                  local_node2=local_node1+1;

                  // Find the second node...
                  while (!neigh->is_node_on_side(local_node2++,s_neigh)) { }
                  local_node2--;

                  // Save them if we haven't already found the parents for this one
                  if (hanging_nodes[hanging_node].size()<2)
                    {
                      hanging_nodes[hanging_node].push_back(neigh->node_id(local_node1));
                      hanging_nodes[hanging_node].push_back(neigh->node_id(local_node2));
                    }
                }
            }
        }
}



void clear_spline_nodes(MeshBase & mesh)
{
  std::vector<Elem *> nodeelem_to_delete;

  for (auto & elem : mesh.element_ptr_range())
    if (elem->type() == NODEELEM &&
        elem->mapping_type() == RATIONAL_BERNSTEIN_MAP)
      nodeelem_to_delete.push_back(elem);

  auto & constraint_rows = mesh.get_constraint_rows();

  // All our constraint_rows ought to be for spline constraints we're
  // about to get rid of.
#ifndef NDEBUG
  for (auto & node_row : constraint_rows)
    for (auto pr : node_row.second)
      {
        const Elem * elem = pr.first.first;
        libmesh_assert(elem->type() == NODEELEM);
        libmesh_assert(elem->mapping_type() == RATIONAL_BERNSTEIN_MAP);
      }
#endif

  constraint_rows.clear();

  for (Elem * elem : nodeelem_to_delete)
    {
      Node * node = elem->node_ptr(0);
      mesh.delete_elem(elem);
      mesh.delete_node(node);
    }
}



bool valid_is_prepared (const MeshBase & mesh)
{
  LOG_SCOPE("valid_is_prepared()", "MeshTools");

  if (!mesh.is_prepared())
    return true;

  std::unique_ptr<MeshBase> mesh_clone = mesh.clone();

  // Try preparing (without allowing repartitioning or renumbering, to
  // avoid false assertion failures)
  bool old_allow_renumbering = mesh_clone->allow_renumbering();
  mesh_clone->allow_renumbering(false);
  bool old_allow_remote_element_removal =
    mesh_clone->allow_remote_element_removal();
  bool old_skip_partitioning = mesh_clone->skip_partitioning();
  mesh_clone->skip_partitioning(true);
  mesh_clone->allow_remote_element_removal(false);
  mesh_clone->prepare_for_use();
  mesh_clone->allow_renumbering(old_allow_renumbering);
  mesh_clone->allow_remote_element_removal(old_allow_remote_element_removal);
  mesh_clone->skip_partitioning(old_skip_partitioning);

  return (mesh == *mesh_clone);
}



#ifndef NDEBUG

void libmesh_assert_equal_n_systems (const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_equal_n_systems()", "MeshTools");

  unsigned int n_sys = libMesh::invalid_uint;

  for (const auto & elem : mesh.element_ptr_range())
    {
      if (n_sys == libMesh::invalid_uint)
        n_sys = elem->n_systems();
      else
        libmesh_assert_equal_to (elem->n_systems(), n_sys);
    }

  for (const auto & node : mesh.node_ptr_range())
    {
      if (n_sys == libMesh::invalid_uint)
        n_sys = node->n_systems();
      else
        libmesh_assert_equal_to (node->n_systems(), n_sys);
    }
}



#ifdef LIBMESH_ENABLE_AMR
void libmesh_assert_old_dof_objects (const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_old_dof_objects()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      if (elem->refinement_flag() == Elem::JUST_REFINED ||
          elem->refinement_flag() == Elem::INACTIVE)
        continue;

      if (elem->has_dofs())
        libmesh_assert(elem->get_old_dof_object());

      for (auto & node : elem->node_ref_range())
        if (node.has_dofs())
          libmesh_assert(node.get_old_dof_object());
    }
}
#else
void libmesh_assert_old_dof_objects (const MeshBase &) {}
#endif // LIBMESH_ENABLE_AMR



void libmesh_assert_valid_node_pointers(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_node_pointers()", "MeshTools");

  // Here we specifically do not want "auto &" because we need to
  // reseat the (temporary) pointer variable in the loop below,
  // without modifying the original.
  for (const Elem * elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);
      while (elem)
        {
          elem->libmesh_assert_valid_node_pointers();
          for (auto n : elem->neighbor_ptr_range())
            if (n && n != remote_elem)
              n->libmesh_assert_valid_node_pointers();

          libmesh_assert_not_equal_to (elem->parent(), remote_elem);
          elem = elem->parent();
        }
    }
}



void libmesh_assert_valid_remote_elems(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_remote_elems()", "MeshTools");

  for (const auto & elem : as_range(mesh.local_elements_begin(),
                                    mesh.local_elements_end()))
    {
      libmesh_assert (elem);

      // We currently don't allow active_local_elements to have
      // remote_elem neighbors
      if (elem->active())
        for (auto n : elem->neighbor_ptr_range())
          libmesh_assert_not_equal_to (n, remote_elem);

#ifdef LIBMESH_ENABLE_AMR
      const Elem * parent = elem->parent();
      if (parent)
        libmesh_assert_not_equal_to (parent, remote_elem);

      // We can only be strict about active elements' subactive
      // children
      if (elem->active() && elem->has_children())
        for (auto & child : elem->child_ref_range())
          libmesh_assert_not_equal_to (&child, remote_elem);
#endif
    }
}



void libmesh_assert_valid_elem_ids(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_elem_ids()", "MeshTools");

  processor_id_type lastprocid = 0;
  dof_id_type lastelemid = 0;

  for (const auto & elem : mesh.active_element_ptr_range())
    {
      libmesh_assert (elem);
      processor_id_type elemprocid = elem->processor_id();
      dof_id_type elemid = elem->id();

      libmesh_assert_greater_equal (elemid, lastelemid);
      libmesh_assert_greater_equal (elemprocid, lastprocid);

      lastelemid = elemid;
      lastprocid = elemprocid;
    }
}



void libmesh_assert_valid_amr_elem_ids(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_amr_elem_ids()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);

      const Elem * parent = elem->parent();

      if (parent)
        {
          libmesh_assert_greater_equal (elem->id(), parent->id());
          libmesh_assert_greater_equal (elem->processor_id(), parent->processor_id());
        }
    }
}



void libmesh_assert_valid_amr_interior_parents(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_amr_interior_parents()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);

      // We can skip to the next element if we're full-dimension
      // and therefore don't have any interior parents
      if (elem->dim() >= LIBMESH_DIM)
        continue;

      const Elem * ip = elem->interior_parent();

      const Elem * parent = elem->parent();

      if (ip && (ip != remote_elem) && parent)
        {
          libmesh_assert_equal_to (ip->top_parent(),
                                   elem->top_parent()->interior_parent());

          if (ip->level() == elem->level())
            libmesh_assert_equal_to (ip->parent(),
                                     parent->interior_parent());
          else
            {
              libmesh_assert_less (ip->level(), elem->level());
              libmesh_assert_equal_to (ip, parent->interior_parent());
            }
        }
    }
}



void libmesh_assert_contiguous_dof_ids(const MeshBase & mesh, unsigned int sysnum)
{
  LOG_SCOPE("libmesh_assert_contiguous_dof_ids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  dof_id_type min_dof_id = std::numeric_limits<dof_id_type>::max(),
              max_dof_id = std::numeric_limits<dof_id_type>::min();

  // Figure out what our local dof id range is
  for (const auto * node : mesh.local_node_ptr_range())
    {
      for (auto v : make_range(node->n_vars(sysnum)))
        for (auto c : make_range(node->n_comp(sysnum, v)))
          {
            dof_id_type id = node->dof_number(sysnum, v, c);
            min_dof_id = std::min (min_dof_id, id);
            max_dof_id = std::max (max_dof_id, id);
          }
    }

  // Make sure no other processors' ids are inside it
  for (const auto * node : mesh.node_ptr_range())
    {
      if (node->processor_id() == mesh.processor_id())
        continue;
      for (auto v : make_range(node->n_vars(sysnum)))
        for (auto c : make_range(node->n_comp(sysnum, v)))
          {
            dof_id_type id = node->dof_number(sysnum, v, c);
            libmesh_assert (id < min_dof_id ||
                            id > max_dof_id);
          }
    }
}



template <>
void libmesh_assert_topology_consistent_procids<Elem>(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_topology_consistent_procids()", "MeshTools");

  // This parameter is not used when !LIBMESH_ENABLE_AMR
  libmesh_ignore(mesh);

  // If we're adaptively refining, check processor ids for consistency
  // between parents and children.
#ifdef LIBMESH_ENABLE_AMR

  // Ancestor elements we won't worry about, but subactive and active
  // elements ought to have parents with consistent processor ids
  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert(elem);

      if (!elem->active() && !elem->subactive())
        continue;

      const Elem * parent = elem->parent();

      if (parent)
        {
          libmesh_assert(parent->has_children());
          processor_id_type parent_procid = parent->processor_id();
          bool matching_child_id = false;
          // If we've got a remote_elem then we don't know whether
          // it's responsible for the parent's processor id; all
          // we can do is assume it is and let its processor fail
          // an assert if there's something wrong.
          for (auto & child : parent->child_ref_range())
            if (&child == remote_elem ||
                child.processor_id() == parent_procid)
              matching_child_id = true;
          libmesh_assert(matching_child_id);
        }
    }
#endif
}



template <>
void libmesh_assert_topology_consistent_procids<Node>(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_topology_consistent_procids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  // We want this test to be valid even when called after nodes have
  // been added asynchronously but before they're renumbered.
  //
  // Plus, some code (looking at you, stitch_meshes) modifies
  // DofObject ids without keeping max_elem_id()/max_node_id()
  // consistent, but that's done in a safe way for performance
  // reasons, so we'll play along and just figure out new max ids
  // ourselves.
  dof_id_type parallel_max_node_id = 0;
  for (const auto & node : mesh.node_ptr_range())
    parallel_max_node_id = std::max<dof_id_type>(parallel_max_node_id,
                                                 node->id()+1);
  mesh.comm().max(parallel_max_node_id);


  std::vector<bool> node_touched_by_me(parallel_max_node_id, false);

  for (const auto & elem : as_range(mesh.local_elements_begin(),
                                    mesh.local_elements_end()))
    {
      libmesh_assert (elem);

      for (auto & node : elem->node_ref_range())
        {
          dof_id_type nodeid = node.id();
          node_touched_by_me[nodeid] = true;
        }
    }
  std::vector<bool> node_touched_by_anyone(node_touched_by_me);
  mesh.comm().max(node_touched_by_anyone);

  for (const auto & node : mesh.local_node_ptr_range())
    {
      libmesh_assert(node);
      dof_id_type nodeid = node->id();
      libmesh_assert(!node_touched_by_anyone[nodeid] ||
                     node_touched_by_me[nodeid]);
    }
}



void libmesh_assert_canonical_node_procids (const MeshBase & mesh)
{
  for (const auto & elem : mesh.active_element_ptr_range())
    for (auto & node : elem->node_ref_range())
      libmesh_assert_equal_to
        (node.processor_id(),
         node.choose_processor_id(node.processor_id(),
                                  elem->processor_id()));
}



#ifdef LIBMESH_ENABLE_AMR
void libmesh_assert_valid_refinement_tree(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_refinement_tree()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert(elem);
      if (elem->has_children())
        for (auto & child : elem->child_ref_range())
          if (&child != remote_elem)
            libmesh_assert_equal_to (child.parent(), elem);
      if (elem->active())
        {
          libmesh_assert(!elem->ancestor());
          libmesh_assert(!elem->subactive());
        }
      else if (elem->ancestor())
        {
          libmesh_assert(!elem->subactive());
        }
      else
        libmesh_assert(elem->subactive());

      if (elem->p_refinement_flag() == Elem::JUST_REFINED)
        libmesh_assert_greater(elem->p_level(), 0);
    }
}
#else
void libmesh_assert_valid_refinement_tree(const MeshBase &)
{
}
#endif // LIBMESH_ENABLE_AMR

#endif // !NDEBUG



#ifdef DEBUG

void libmesh_assert_no_links_to_elem(const MeshBase & mesh,
                                     const Elem * bad_elem)
{
  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);
      libmesh_assert_not_equal_to (elem->parent(), bad_elem);
      for (auto n : elem->neighbor_ptr_range())
        libmesh_assert_not_equal_to (n, bad_elem);

#ifdef LIBMESH_ENABLE_AMR
      if (elem->has_children())
        for (auto & child : elem->child_ref_range())
          libmesh_assert_not_equal_to (&child, bad_elem);
#endif
    }
}


void libmesh_assert_equal_points (const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_equal_points()", "MeshTools");

  dof_id_type pmax_node_id = mesh.max_node_id();
  mesh.comm().max(pmax_node_id);

  for (dof_id_type i=0; i != pmax_node_id; ++i)
    {
      const Point * p = mesh.query_node_ptr(i);

      libmesh_assert(mesh.comm().semiverify(p));
    }
}


void libmesh_assert_equal_connectivity (const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_equal_connectivity()", "MeshTools");

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      const Elem * e = mesh.query_elem_ptr(i);

      std::vector<dof_id_type> nodes;
      if (e)
        for (auto n : e->node_index_range())
          nodes.push_back(e->node_id(n));

      libmesh_assert(mesh.comm().semiverify(e ? &nodes : nullptr));
    }
}


void libmesh_assert_connected_nodes (const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_connected_nodes()", "MeshTools");

  std::set<const Node *> used_nodes;

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);

      for (auto & n : elem->node_ref_range())
        used_nodes.insert(&n);
    }

  for (const auto & node : mesh.node_ptr_range())
    {
      libmesh_assert(node);
      libmesh_assert(used_nodes.count(node));
    }
}



void libmesh_assert_valid_constraint_rows (const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  const auto & constraint_rows = mesh.get_constraint_rows();

  bool have_constraint_rows = !constraint_rows.empty();
  mesh.comm().max(have_constraint_rows);
  if (!have_constraint_rows)
    return;

  for (auto & row : constraint_rows)
    {
      const Node * node = row.first;
      libmesh_assert(node == mesh.node_ptr(node->id()));

      for (auto & pr : row.second)
        {
          const Elem * spline_elem = pr.first.first;
          libmesh_assert(spline_elem == mesh.elem_ptr(spline_elem->id()));
        }
    }

  dof_id_type pmax_node_id = mesh.max_node_id();
  mesh.comm().max(pmax_node_id);

  for (dof_id_type i=0; i != pmax_node_id; ++i)
    {
      const Node * node = mesh.query_node_ptr(i);

      bool have_constraint = constraint_rows.count(node);

      const std::size_t my_n_constraints = have_constraint ?
        libmesh_map_find(constraint_rows, node).size() : std::size_t(-1);
      const std::size_t * n_constraints = node ?
        &my_n_constraints : nullptr;

      libmesh_assert(mesh.comm().semiverify(n_constraints));
    }
}



void libmesh_assert_valid_boundary_ids(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_boundary_ids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  const BoundaryInfo & boundary_info = mesh.get_boundary_info();

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      const Elem * elem = mesh.query_elem_ptr(i);
      const unsigned int my_n_nodes = elem ? elem->n_nodes() : 0;
      const unsigned int my_n_edges = elem ? elem->n_edges() : 0;
      const unsigned int my_n_sides = elem ? elem->n_sides() : 0;
      unsigned int
        n_nodes = my_n_nodes,
        n_edges = my_n_edges,
        n_sides = my_n_sides;

      mesh.comm().max(n_nodes);
      mesh.comm().max(n_edges);
      mesh.comm().max(n_sides);

      if (elem)
        {
          libmesh_assert_equal_to(my_n_nodes, n_nodes);
          libmesh_assert_equal_to(my_n_edges, n_edges);
          libmesh_assert_equal_to(my_n_sides, n_sides);
        }

      // Let's test all IDs on the element with one communication
      // rather than n_nodes + n_edges + n_sides communications, to
      // cut down on latency in dbg modes.
      std::vector<boundary_id_type> all_bcids;

      for (unsigned int n=0; n != n_nodes; ++n)
        {
          std::vector<boundary_id_type> bcids;
          if (elem)
            {
              boundary_info.boundary_ids(elem->node_ptr(n), bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }
          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));

          all_bcids.insert(all_bcids.end(), bcids.begin(),
                           bcids.end());
          // Separator
          all_bcids.push_back(BoundaryInfo::invalid_id);
        }

      for (unsigned short e=0; e != n_edges; ++e)
        {
          std::vector<boundary_id_type> bcids;

          if (elem)
            {
              boundary_info.edge_boundary_ids(elem, e, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));

          all_bcids.insert(all_bcids.end(), bcids.begin(),
                           bcids.end());
          // Separator
          all_bcids.push_back(BoundaryInfo::invalid_id);

          if (elem)
            {
              boundary_info.raw_edge_boundary_ids(elem, e, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));
        }

      for (unsigned short s=0; s != n_sides; ++s)
        {
          std::vector<boundary_id_type> bcids;

          if (elem)
            {
              boundary_info.boundary_ids(elem, s, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));

          if (elem)
            {
              boundary_info.raw_boundary_ids(elem, s, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));
        }

      for (unsigned short sf=0; sf != 2; ++sf)
        {
          std::vector<boundary_id_type> bcids;

          if (elem)
            {
              boundary_info.shellface_boundary_ids(elem, sf, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));

          if (elem)
            {
              boundary_info.raw_shellface_boundary_ids(elem, sf, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));
        }

      libmesh_assert(mesh.comm().semiverify
                     (elem ? &all_bcids : nullptr));
    }
}


void libmesh_assert_valid_dof_ids(const MeshBase & mesh, unsigned int sysnum)
{
  LOG_SCOPE("libmesh_assert_valid_dof_ids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    assert_semiverify_dofobj(mesh.comm(),
                             mesh.query_elem_ptr(i),
                             sysnum);

  dof_id_type pmax_node_id = mesh.max_node_id();
  mesh.comm().max(pmax_node_id);

  for (dof_id_type i=0; i != pmax_node_id; ++i)
    assert_semiverify_dofobj(mesh.comm(),
                             mesh.query_node_ptr(i),
                             sysnum);
}


#ifdef LIBMESH_ENABLE_UNIQUE_ID
void libmesh_assert_valid_unique_ids(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_unique_ids()", "MeshTools");

  libmesh_parallel_only(mesh.comm());

  // Storage for semi-local DofObject ids.
  std::unordered_set<unique_id_type> semilocal_unique_ids;

  auto gather_elem_ids = [&]()
  {
    for (auto const & elem : mesh.active_element_ptr_range())
      {
        auto [it, inserted] = semilocal_unique_ids.insert(elem->unique_id());
        libmesh_assert(inserted);
        libmesh_ignore(it);
      }
  };

  auto gather_node_ids = [&]()
  {
    for (auto const & node : mesh.node_ptr_range())
      {
        auto [it, inserted] = semilocal_unique_ids.insert(node->unique_id());
        libmesh_assert(inserted);
        libmesh_ignore(it);
      }
  };

  auto verify_elems = [&]()
  {
    dof_id_type pmax_elem_id = mesh.max_elem_id();
    mesh.comm().max(pmax_elem_id);

    for (auto i : make_range(pmax_elem_id))
      {
        const Elem * elem = mesh.query_elem_ptr(i);
        assert_dofobj_unique_id(mesh.comm(), elem, semilocal_unique_ids);
      }
  };

  auto verify_nodes = [&]()
  {
    dof_id_type pmax_node_id = mesh.max_node_id();
    mesh.comm().max(pmax_node_id);

    for (auto i : make_range(pmax_node_id))
      {
        const Node * node = mesh.query_node_ptr(i);
        assert_dofobj_unique_id(mesh.comm(), node, semilocal_unique_ids);
      }
  };

  if (!mesh.allow_node_and_elem_unique_id_overlap())
  {
    // First collect all the unique_ids we can see and make sure there's
    // no duplicates
    gather_elem_ids();
    gather_node_ids();

    // Then make sure elements/nodes are all in sync and remote
    // elements don't duplicate semilocal
    verify_elems();
    verify_nodes();
  }
  else
  {
    // If the mesh allows Node and Elem unique_ids to overlap, then we only
    // check for validity and uniqueness of an Elem (resp. Node) unique id
    // within the set of Elem (resp. Node) unique_ids.
    gather_elem_ids();
    verify_elems();

    // Clear id list before checking Nodes
    semilocal_unique_ids.clear();

    // Finally, check Nodes
    gather_node_ids();
    verify_nodes();
  }
}
#endif

void libmesh_assert_consistent_distributed(const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());

  dof_id_type parallel_max_elem_id = mesh.max_elem_id();
  mesh.comm().max(parallel_max_elem_id);

  for (dof_id_type i=0; i != parallel_max_elem_id; ++i)
    {
      const Elem * elem = mesh.query_elem_ptr(i);
      processor_id_type pid =
        elem ? elem->processor_id() : DofObject::invalid_processor_id;
      mesh.comm().min(pid);
      libmesh_assert(elem || pid != mesh.processor_id());
    }

  dof_id_type parallel_max_node_id = mesh.max_node_id();
  mesh.comm().max(parallel_max_node_id);

  for (dof_id_type i=0; i != parallel_max_node_id; ++i)
    {
      const Node * node = mesh.query_node_ptr(i);
      processor_id_type pid =
        node ? node->processor_id() : DofObject::invalid_processor_id;
      mesh.comm().min(pid);
      libmesh_assert(node || pid != mesh.processor_id());
    }
}


void libmesh_assert_consistent_distributed_nodes(const MeshBase & mesh)
{
  libmesh_parallel_only(mesh.comm());
  auto locator = mesh.sub_point_locator();

  dof_id_type parallel_max_elem_id = mesh.max_elem_id();
  mesh.comm().max(parallel_max_elem_id);

  for (dof_id_type i=0; i != parallel_max_elem_id; ++i)
    {
      const Elem * elem = mesh.query_elem_ptr(i);

      const unsigned int my_n_nodes = elem ? elem->n_nodes() : 0;
      unsigned int n_nodes = my_n_nodes;
      mesh.comm().max(n_nodes);

      if (n_nodes)
        libmesh_assert(mesh.comm().semiverify(elem ? &my_n_nodes : nullptr));

      for (unsigned int n=0; n != n_nodes; ++n)
        {
          const Node * node = elem ? elem->node_ptr(n) : nullptr;
          processor_id_type pid =
            node ? node->processor_id() : DofObject::invalid_processor_id;
          mesh.comm().min(pid);
          libmesh_assert(node || pid != mesh.processor_id());
        }
    }
}



template <>
void libmesh_assert_parallel_consistent_procids<Elem>(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_parallel_consistent_procids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  // Some code (looking at you, stitch_meshes) modifies DofObject ids
  // without keeping max_elem_id()/max_node_id() consistent, but
  // that's done in a safe way for performance reasons, so we'll play
  // along and just figure out new max ids ourselves.
  dof_id_type parallel_max_elem_id = 0;
  for (const auto & elem : mesh.element_ptr_range())
    parallel_max_elem_id = std::max<dof_id_type>(parallel_max_elem_id,
                                                 elem->id()+1);
  mesh.comm().max(parallel_max_elem_id);

  // Check processor ids for consistency between processors

  for (dof_id_type i=0; i != parallel_max_elem_id; ++i)
    {
      const Elem * elem = mesh.query_elem_ptr(i);

      processor_id_type min_id =
        elem ? elem->processor_id() :
        std::numeric_limits<processor_id_type>::max();
      mesh.comm().min(min_id);

      processor_id_type max_id =
        elem ? elem->processor_id() :
        std::numeric_limits<processor_id_type>::min();
      mesh.comm().max(max_id);

      if (elem)
        {
          libmesh_assert_equal_to (min_id, elem->processor_id());
          libmesh_assert_equal_to (max_id, elem->processor_id());
        }

      if (min_id == mesh.processor_id())
        libmesh_assert(elem);
    }
}



void libmesh_assert_parallel_consistent_new_node_procids(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_parallel_consistent_new_node_procids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  // We want this test to hit every node when called even after nodes
  // have been added asynchronously but before everything has been
  // renumbered.
  dof_id_type parallel_max_elem_id = mesh.max_elem_id();
  mesh.comm().max(parallel_max_elem_id);

  std::vector<bool> elem_touched_by_anyone(parallel_max_elem_id, false);

  for (dof_id_type i=0; i != parallel_max_elem_id; ++i)
    {
      const Elem * elem = mesh.query_elem_ptr(i);

      const unsigned int my_n_nodes = elem ? elem->n_nodes() : 0;
      unsigned int n_nodes = my_n_nodes;
      mesh.comm().max(n_nodes);

      if (n_nodes)
        libmesh_assert(mesh.comm().semiverify(elem ? &my_n_nodes : nullptr));

      for (unsigned int n=0; n != n_nodes; ++n)
        {
          const Node * node = elem ? elem->node_ptr(n) : nullptr;
          const processor_id_type pid = node ? node->processor_id() : 0;
          libmesh_assert(mesh.comm().semiverify (node ? &pid : nullptr));
        }
    }
}

template <>
void libmesh_assert_parallel_consistent_procids<Node>(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_parallel_consistent_procids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  // We want this test to be valid even when called even after nodes
  // have been added asynchronously but before they're renumbered
  //
  // Plus, some code (looking at you, stitch_meshes) modifies
  // DofObject ids without keeping max_elem_id()/max_node_id()
  // consistent, but that's done in a safe way for performance
  // reasons, so we'll play along and just figure out new max ids
  // ourselves.
  dof_id_type parallel_max_node_id = 0;
  for (const auto & node : mesh.node_ptr_range())
    parallel_max_node_id = std::max<dof_id_type>(parallel_max_node_id,
                                                 node->id()+1);
  mesh.comm().max(parallel_max_node_id);

  std::vector<bool> node_touched_by_anyone(parallel_max_node_id, false);

  for (const auto & elem : as_range(mesh.local_elements_begin(),
                                    mesh.local_elements_end()))
    {
      libmesh_assert (elem);

      for (auto & node : elem->node_ref_range())
        {
          dof_id_type nodeid = node.id();
          node_touched_by_anyone[nodeid] = true;
        }
    }
  mesh.comm().max(node_touched_by_anyone);

  // Check processor ids for consistency between processors
  // on any node an element touches
  for (dof_id_type i=0; i != parallel_max_node_id; ++i)
    {
      if (!node_touched_by_anyone[i])
        continue;

      const Node * node = mesh.query_node_ptr(i);
      const processor_id_type pid = node ? node->processor_id() : 0;

      libmesh_assert(mesh.comm().semiverify (node ? &pid : nullptr));
    }
}



#ifdef LIBMESH_ENABLE_AMR
void libmesh_assert_valid_refinement_flags(const MeshBase & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_refinement_flags()", "MeshTools");

  libmesh_parallel_only(mesh.comm());
  if (mesh.n_processors() == 1)
    return;

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  std::vector<unsigned char> my_elem_h_state(pmax_elem_id, 255);
  std::vector<unsigned char> my_elem_p_state(pmax_elem_id, 255);

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);
      dof_id_type elemid = elem->id();

      my_elem_h_state[elemid] =
        static_cast<unsigned char>(elem->refinement_flag());

      my_elem_p_state[elemid] =
        static_cast<unsigned char>(elem->p_refinement_flag());
    }
  std::vector<unsigned char> min_elem_h_state(my_elem_h_state);
  mesh.comm().min(min_elem_h_state);

  std::vector<unsigned char> min_elem_p_state(my_elem_p_state);
  mesh.comm().min(min_elem_p_state);

  for (dof_id_type i=0; i!= pmax_elem_id; ++i)
    {
      libmesh_assert(my_elem_h_state[i] == 255 ||
                     my_elem_h_state[i] == min_elem_h_state[i]);
      libmesh_assert(my_elem_p_state[i] == 255 ||
                     my_elem_p_state[i] == min_elem_p_state[i]);
    }
}
#else
void libmesh_assert_valid_refinement_flags(const MeshBase &)
{
}
#endif // LIBMESH_ENABLE_AMR



void libmesh_assert_valid_neighbors(const MeshBase & mesh,
                                    bool assert_valid_remote_elems)
{
  LOG_SCOPE("libmesh_assert_valid_neighbors()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);
      elem->libmesh_assert_valid_neighbors();
    }

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      const Elem * elem = mesh.query_elem_ptr(i);

      const unsigned int my_n_neigh = elem ? elem->n_neighbors() : 0;
      unsigned int n_neigh = my_n_neigh;
      mesh.comm().max(n_neigh);
      if (elem)
        libmesh_assert_equal_to (my_n_neigh, n_neigh);

      for (unsigned int n = 0; n != n_neigh; ++n)
        {
          dof_id_type my_neighbor = DofObject::invalid_id;
          dof_id_type * p_my_neighbor = nullptr;

          // If we have a non-remote_elem neighbor link, then we can
          // verify it.
          if (elem && elem->neighbor_ptr(n) != remote_elem)
            {
              p_my_neighbor = &my_neighbor;
              if (elem->neighbor_ptr(n))
                my_neighbor = elem->neighbor_ptr(n)->id();

              // But wait - if we haven't set remote_elem links yet then
              // some nullptr links on ghost elements might be
              // future-remote_elem links, so we can't verify those.
              if (!assert_valid_remote_elems &&
                  !elem->neighbor_ptr(n) &&
                  elem->processor_id() != mesh.processor_id())
                p_my_neighbor = nullptr;
            }
          libmesh_assert(mesh.comm().semiverify(p_my_neighbor));
        }
    }
}
#endif // DEBUG



// Functors for correct_node_proc_ids
namespace {

typedef std::unordered_map<dof_id_type, processor_id_type> proc_id_map_type;

struct SyncNodeSet
{
  typedef unsigned char datum; // bool but without bit twiddling issues

  SyncNodeSet(std::unordered_set<const Node *> & _set,
              MeshBase & _mesh) :
    node_set(_set), mesh(_mesh) {}

  std::unordered_set<const Node *> & node_set;

  MeshBase & mesh;

  // ------------------------------------------------------------
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & data)
  {
    // Find whether each requested node belongs in the set
    data.resize(ids.size());

    for (auto i : index_range(ids))
      {
        const dof_id_type id = ids[i];

        // We'd better have every node we're asked for
        Node * node = mesh.node_ptr(id);

        // Return if the node is in the set.
        data[i] = node_set.count(node);
      }
  }

  // ------------------------------------------------------------
  bool act_on_data (const std::vector<dof_id_type> & ids,
                    const std::vector<datum> in_set)
  {
    bool data_changed = false;

    // Add nodes we've been informed of to our own set
    for (auto i : index_range(ids))
      {
        if (in_set[i])
          {
            Node * node = mesh.node_ptr(ids[i]);
            if (!node_set.count(node))
              {
                node_set.insert(node);
                data_changed = true;
              }
          }
      }

    return data_changed;
  }
};


struct NodesNotInSet
{
  NodesNotInSet(const std::unordered_set<const Node *> _set)
    : node_set(_set) {}

  bool operator() (const Node * node) const
  {
    if (node_set.count(node))
      return false;
    return true;
  }

  const std::unordered_set<const Node *> node_set;
};


struct SyncProcIdsFromMap
{
  typedef processor_id_type datum;

  SyncProcIdsFromMap(const proc_id_map_type & _map,
                     MeshBase & _mesh) :
    new_proc_ids(_map), mesh(_mesh) {}

  const proc_id_map_type & new_proc_ids;

  MeshBase & mesh;

  // ------------------------------------------------------------
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & data)
  {
    // Find the new processor id of each requested node
    data.resize(ids.size());

    for (auto i : index_range(ids))
      {
        const dof_id_type id = ids[i];

        // Return the node's new processor id if it has one, or its
        // old processor id if not.
        if (const auto it = new_proc_ids.find(id);
            it != new_proc_ids.end())
          data[i] = it->second;
        else
          {
            // We'd better find every node we're asked for
            const Node & node = mesh.node_ref(id);
            data[i] = node.processor_id();
          }
      }
  }

  // ------------------------------------------------------------
  void act_on_data (const std::vector<dof_id_type> & ids,
                    const std::vector<datum> proc_ids)
  {
    // Set the node processor ids we've now been informed of
    for (auto i : index_range(ids))
      {
        Node & node = mesh.node_ref(ids[i]);
        node.processor_id() = proc_ids[i];
      }
  }
};
}



void correct_node_proc_ids (MeshBase & mesh)
{
  LOG_SCOPE("correct_node_proc_ids()","MeshTools");

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  // We require all processors to agree on nodal processor ids before
  // going into this algorithm.
#ifdef DEBUG
  libmesh_assert_parallel_consistent_procids<Node>(mesh);
#endif

  // If we have any unpartitioned elements at this
  // stage there is a problem
  libmesh_assert (n_elem(mesh.unpartitioned_elements_begin(),
                         mesh.unpartitioned_elements_end()) == 0);

  // Fix nodes' processor ids.  Coarsening may have left us with nodes
  // which are no longer touched by any elements of the same processor
  // id, and for DofMap to work we need to fix that.

  // This is harder now that libMesh no longer requires a distributed
  // mesh to ghost all nodal neighbors: it is possible for two active
  // elements on two different processors to share the same node in
  // such a way that neither processor knows the others' element
  // exists!

  // While we're at it, if this mesh is configured to allow
  // repartitioning, we'll repartition *all* the nodes' processor ids
  // using the canonical Node heuristic, to try and improve DoF load
  // balancing.  But if the mesh is disallowing repartitioning, we
  // won't touch processor_id on any node where it's valid, regardless
  // of whether or not it's canonical.
  bool repartition_all_nodes = !mesh.skip_noncritical_partitioning();
  std::unordered_set<const Node *> valid_nodes;

  // If we aren't allowed to repartition, then we're going to leave
  // every node we can at its current processor_id, and *only*
  // repartition the nodes whose current processor id is incompatible
  // with DoFMap (because it doesn't touch an active element, e.g. due
  // to coarsening)
  if (!repartition_all_nodes)
    {
      for (const auto & elem : mesh.active_element_ptr_range())
        for (const auto & node : elem->node_ref_range())
          if (elem->processor_id() == node.processor_id())
            valid_nodes.insert(&node);

      SyncNodeSet syncv(valid_nodes, mesh);

      Parallel::sync_dofobject_data_by_id
        (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), syncv);
    }

  // We build up a set of compatible processor ids for each node
  proc_id_map_type new_proc_ids;

  for (auto & elem : mesh.active_element_ptr_range())
    {
      processor_id_type pid = elem->processor_id();

      for (auto & node : elem->node_ref_range())
        {
          const dof_id_type id = node.id();
          if (auto it = new_proc_ids.find(id);
              it == new_proc_ids.end())
            new_proc_ids.emplace(id, pid);
          else
            it->second = node.choose_processor_id(it->second, pid);
        }
    }

  // Sort the new pids to push to each processor
  std::map<processor_id_type, std::vector<std::pair<dof_id_type, processor_id_type>>>
    ids_to_push;

  for (const auto & node : mesh.node_ptr_range())
    if (const auto it = std::as_const(new_proc_ids).find(node->id());
        it != new_proc_ids.end() && node->processor_id() != DofObject::invalid_processor_id)
      ids_to_push[node->processor_id()].emplace_back(node->id(), /*pid=*/it->second);

  auto action_functor =
    [& mesh, & new_proc_ids]
    (processor_id_type,
     const std::vector<std::pair<dof_id_type, processor_id_type>> & data)
    {
      for (const auto & [id, pid] : data)
        {
          if (const auto it = new_proc_ids.find(id);
              it == new_proc_ids.end())
            new_proc_ids.emplace(id, pid);
          else
            {
              const Node & node = mesh.node_ref(id);
              it->second = node.choose_processor_id(it->second, pid);
            }
        }
    };

  Parallel::push_parallel_vector_data
    (mesh.comm(), ids_to_push, action_functor);

  // Now new_proc_ids is correct for every node we used to own.  Let's
  // ask every other processor about the nodes they used to own.  But
  // first we'll need to keep track of which nodes we used to own,
  // lest we get them confused with nodes we newly own.
  std::unordered_set<Node *> ex_local_nodes;
  for (auto & node : mesh.local_node_ptr_range())
    if (const auto it = new_proc_ids.find(node->id());
        it != new_proc_ids.end() && it->second != mesh.processor_id())
      ex_local_nodes.insert(node);

  SyncProcIdsFromMap sync(new_proc_ids, mesh);
  if (repartition_all_nodes)
    Parallel::sync_dofobject_data_by_id
      (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), sync);
  else
    {
      NodesNotInSet nnis(valid_nodes);

      Parallel::sync_dofobject_data_by_id
        (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), nnis, sync);
    }

  // And finally let's update the nodes we used to own.
  for (const auto & node : ex_local_nodes)
    {
      if (valid_nodes.count(node))
        continue;

      const dof_id_type id = node->id();
      const proc_id_map_type::iterator it = new_proc_ids.find(id);
      libmesh_assert(it != new_proc_ids.end());
      node->processor_id() = it->second;
    }

  // We should still have consistent nodal processor ids coming out of
  // this algorithm, but if we're allowed to repartition the mesh then
  // they should be canonically correct too.
#ifdef DEBUG
  libmesh_assert_valid_procids<Node>(mesh);
  //if (repartition_all_nodes)
  //  libmesh_assert_canonical_node_procids(mesh);
#endif
}



void Private::globally_renumber_nodes_and_elements (MeshBase & mesh)
{
  MeshCommunication().assign_global_indices(mesh);
}

} // namespace MeshTools

} // namespace libMesh
