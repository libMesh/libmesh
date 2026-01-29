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

#ifndef LIBMESH_ELEM_INTERNAL_H
#define LIBMESH_ELEM_INTERNAL_H

// libMesh includes
#include "libmesh/elem.h"

namespace libMesh
{

// Forward declarations
class PointLocatorBase;
class PeriodicBoundaries;

/**
 * The ElemInternal namespace holds helper functions that are used
 * internally by the Elem class. These should not be called directly,
 * call the appropriate member functions on the Elem class instead.
 */
namespace ElemInternal
{

template<class T>
void
family_tree (T elem,
             std::vector<T> & family,
             bool reset = true)
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!elem->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // Add this element to the family tree.
  family.push_back(elem);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!elem->active())
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote())
        ElemInternal::family_tree (&c, family, false);
}



template<class T>
void
total_family_tree(T elem,
                  std::vector<T> & family,
                  bool reset)
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // Add this element to the family tree.
  family.push_back(elem);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (elem->has_children())
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote())
        ElemInternal::total_family_tree (&c, family, false);
}



template<class T>
void
active_family_tree(T elem,
                   std::vector<T> & active_family,
                   bool reset = true)
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!elem->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    active_family.clear();

  // Add this element to the family tree if it is active
  if (elem->active())
    active_family.push_back(elem);

  // Otherwise recurse into the element's children.
  // Do not clear the vector any more.
  else
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote())
        ElemInternal::active_family_tree (&c, active_family, false);
}



template <class T>
void
family_tree_by_side (T elem,
                     std::vector<T> & family,
                     unsigned int s,
                     bool reset)
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!elem->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  libmesh_assert_less (s, elem->n_sides());

  // Add this element to the family tree.
  family.push_back(elem);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!elem->active())
    {
      const unsigned int nc = elem->n_children();
      for (unsigned int c = 0; c != nc; c++)
        if (!elem->child_ptr(c)->is_remote() && elem->is_child_on_side(c, s))
          ElemInternal::family_tree_by_side (elem->child_ptr(c), family, s, false);
    }
}



template <class T>
void
total_family_tree_by_side (T elem,
                           std::vector<T> & family,
                           unsigned int s,
                           bool reset)
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  libmesh_assert_less (s, elem->n_sides());

  // Add this element to the family tree.
  family.push_back(elem);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (elem->has_children())
    {
      const unsigned int nc = elem->n_children();
      for (unsigned int c = 0; c != nc; c++)
        if (!elem->child_ptr(c)->is_remote() && elem->is_child_on_side(c, s))
          ElemInternal::total_family_tree_by_side (elem->child_ptr(c), family, s, false);
    }
}



template <class T>
void
active_family_tree_by_side (T elem,
                            std::vector<T> & family,
                            unsigned int side,
                            bool reset = true)
{
  // The "family tree" doesn't include subactive or remote elements
  libmesh_assert(!elem->subactive());
  libmesh_assert(!elem->is_remote());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  libmesh_assert_less (side, elem->n_sides());

  // Add an active element to the family tree.
  if (elem->active())
    family.push_back(elem);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else
    {
      const unsigned int nc = elem->n_children();
      for (unsigned int c = 0; c != nc; c++)
        if (!elem->child_ptr(c)->is_remote() && elem->is_child_on_side(c, side))
          ElemInternal::active_family_tree_by_side(elem->child_ptr(c), family, side, false);
    }
}


template <class T>
void
family_tree_by_neighbor (T elem,
                         std::vector<T> & family,
                         T neighbor_in,
                         bool reset = true)
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!elem->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  libmesh_assert (elem->has_neighbor(neighbor_in));

  // Add this element to the family tree.
  family.push_back(elem);

  // Recurse into the elements children, if it's not active.
  // Do not clear the vector any more.
  if (!elem->active())
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote() && c.has_neighbor(neighbor_in))
        ElemInternal::family_tree_by_neighbor (&c, family, neighbor_in, false);
}


template <class T>
void
total_family_tree_by_neighbor (T elem,
                               std::vector<T> & family,
                               T neighbor_in,
                               bool reset = true)
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  libmesh_assert (elem->has_neighbor(neighbor_in));

  // Add this element to the family tree.
  family.push_back(elem);

  // Recurse into the elements children, if it has any.
  // Do not clear the vector any more.
  if (elem->has_children())
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote() && c.has_neighbor(neighbor_in))
        ElemInternal::total_family_tree_by_neighbor (&c, family, neighbor_in, false);
}


template <class T>
void
family_tree_by_subneighbor (T elem,
                            std::vector<T> & family,
                            T neighbor_in,
                            T subneighbor,
                            bool reset = true)
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!elem->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // To simplify this function we need an existing neighbor
  libmesh_assert (neighbor_in);
  libmesh_assert (!neighbor_in->is_remote());
  libmesh_assert (elem->has_neighbor(neighbor_in));

  // This only makes sense if subneighbor descends from neighbor
  libmesh_assert (subneighbor);
  libmesh_assert (!subneighbor->is_remote());
  libmesh_assert (neighbor_in->is_ancestor_of(subneighbor));

  // Add this element to the family tree if applicable.
  if (neighbor_in == subneighbor)
    family.push_back(elem);

  // Recurse into the elements children, if it's not active.
  // Do not clear the vector any more.
  if (!elem->active())
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote())
        for (auto child_neigh : c.neighbor_ptr_range())
          if (child_neigh &&
              (child_neigh == neighbor_in || (child_neigh->parent() == neighbor_in &&
                                              child_neigh->is_ancestor_of(subneighbor))))
            ElemInternal::family_tree_by_subneighbor (&c, family, child_neigh, subneighbor, false);
}



template <class T>
void
total_family_tree_by_subneighbor(T elem,
                                 std::vector<T> & family,
                                 T neighbor_in,
                                 T subneighbor,
                                 bool reset = true)
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // To simplify this function we need an existing neighbor
  libmesh_assert (neighbor_in);
  libmesh_assert (!neighbor_in->is_remote());
  libmesh_assert (elem->has_neighbor(neighbor_in));

  // This only makes sense if subneighbor descends from neighbor
  libmesh_assert (subneighbor);
  libmesh_assert (!subneighbor->is_remote());
  libmesh_assert (neighbor_in->is_ancestor_of(subneighbor));

  // Add this element to the family tree if applicable.
  if (neighbor_in == subneighbor)
    family.push_back(elem);

  // Recurse into the elements children, if it has any.
  // Do not clear the vector any more.
  if (elem->has_children())
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote())
        for (auto child_neigh : c.neighbor_ptr_range())
          if (child_neigh &&
              (child_neigh == neighbor_in ||
               (child_neigh->parent() == neighbor_in && child_neigh->is_ancestor_of(subneighbor))))
            c.total_family_tree_by_subneighbor(family, child_neigh, subneighbor, false);
}



template <class T>
void
active_family_tree_by_neighbor(T elem,
                               std::vector<T> & family,
                               T neighbor_in,
                               bool reset = true)
{
  // The "family tree" doesn't include subactive elements or
  // remote_elements
  libmesh_assert(!elem->subactive());
  libmesh_assert(!elem->is_remote());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
#ifndef NDEBUG
  if (elem->level() >= neighbor_in->level())
    libmesh_assert (elem->has_neighbor(neighbor_in));
#endif

  // Add an active element to the family tree.
  if (elem->active())
    family.push_back(elem);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else if (!elem->active())
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote() && c.has_neighbor(neighbor_in))
        ElemInternal::active_family_tree_by_neighbor (&c, family, neighbor_in, false);
}



template <class T>
void
active_family_tree_by_topological_neighbor (T elem,
                                            std::vector<T> & family,
                                            T neighbor_in,
                                            const MeshBase & mesh,
                                            const PointLocatorBase & point_locator,
                                            const PeriodicBoundaries * pb,
                                            bool reset = true)
{
  // The "family tree" doesn't include subactive elements or
  // remote_elements
  libmesh_assert(!elem->subactive());
  libmesh_assert(!elem->is_remote());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a topological neighbor
#ifndef NDEBUG
  if (elem->level() >= neighbor_in->level())
    libmesh_assert (elem->has_topological_neighbor(neighbor_in,
                                                   mesh,
                                                   point_locator,
                                                   pb));
#endif

  // Add an active element to the family tree.
  if (elem->active())
    family.push_back(elem);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else if (!elem->active())
    for (auto & c : elem->child_ref_range())
      if (!c.is_remote() &&
          c.has_topological_neighbor(neighbor_in, mesh, point_locator, pb))
        c.active_family_tree_by_topological_neighbor
          (family, neighbor_in, mesh, point_locator, pb, false);
}



template <class T>
void
find_point_neighbors(T this_elem,
                     std::set<T> & neighbor_set,
                     T start_elem)
{
  libmesh_assert(start_elem);
  libmesh_assert(start_elem->active());
  libmesh_assert(start_elem->contains_vertex_of(this_elem, true) ||
                 this_elem->contains_vertex_of(start_elem, true));

  neighbor_set.clear();
  neighbor_set.insert(start_elem);

  std::set<T> untested_set, next_untested_set;
  untested_set.insert(start_elem);

  while (!untested_set.empty())
    {
      // Loop over all the elements in the patch that haven't already
      // been tested
      for (const auto & elem : untested_set)
          for (auto current_neighbor : elem->neighbor_ptr_range())
            {
              if (current_neighbor &&
                  !current_neighbor->is_remote())    // we have a real neighbor on this side
                {
                  if (current_neighbor->active())                // ... if it is active
                    {
                      if (this_elem->contains_vertex_of(current_neighbor, true) // ... and touches us
                          || current_neighbor->contains_vertex_of(this_elem, true))
                        {
                          // Make sure we'll test it
                          if (!neighbor_set.count(current_neighbor))
                            next_untested_set.insert (current_neighbor);

                          // And add it
                          neighbor_set.insert (current_neighbor);
                        }
                    }
#ifdef LIBMESH_ENABLE_AMR
                  else                                 // ... the neighbor is *not* active,
                    {                                  // ... so add *all* neighboring
                                                       // active children
                      std::vector<T> active_neighbor_children;

                      current_neighbor->active_family_tree_by_neighbor
                        (active_neighbor_children, elem);

                      for (const auto & current_child : active_neighbor_children)
                        {
                          if (this_elem->contains_vertex_of(current_child, true) ||
                              current_child->contains_vertex_of(this_elem, true))
                            {
                              // Make sure we'll test it
                              if (!neighbor_set.count(current_child))
                                next_untested_set.insert (current_child);

                              neighbor_set.insert (current_child);
                            }
                        }
                    }
#endif // #ifdef LIBMESH_ENABLE_AMR
                }
            }
      untested_set.swap(next_untested_set);
      next_untested_set.clear();
    }
}



template <class T>
void
find_interior_neighbors(T this_elem,
                        std::set<T> & neighbor_set)
{
  neighbor_set.clear();

  if ((this_elem->dim() >= LIBMESH_DIM) ||
      !this_elem->interior_parent())
    return;

  T ip = this_elem->interior_parent();
  libmesh_assert (ip->contains_vertex_of(this_elem, true) ||
                  this_elem->contains_vertex_of(ip, true));

  libmesh_assert (!ip->subactive());

#ifdef LIBMESH_ENABLE_AMR
  while (!ip->active()) // only possible with AMR, be careful because
    {                   // ip->child_ptr(c) is only good with AMR.
      for (auto & child : ip->child_ref_range())
        {
          if (child.contains_vertex_of(this_elem, true) ||
              this_elem->contains_vertex_of(&child, true))
            {
              ip = &child;
              break;
            }
        }
    }
#endif

  ElemInternal::find_point_neighbors(this_elem, neighbor_set, ip);

  // Now we have all point neighbors from the interior manifold, but
  // we need to weed out any neighbors that *only* intersect us at one
  // point (or at one edge, if we're a 1-D element in 3D).
  //
  // The refinement hierarchy helps us here: if the interior element
  // has a lower or equal refinement level then we can discard it iff
  // it doesn't contain all our vertices.  If it has a higher
  // refinement level then we can discard it iff we don't contain at
  // least dim()+1 of its vertices
  auto it = neighbor_set.begin();
  const auto end = neighbor_set.end();

  while (it != end)
    {
      auto current = it++;
      T current_elem = *current;

      // This won't invalidate iterator it, because it is already
      // pointing to the next element
      if (current_elem->level() > this_elem->level())
        {
          unsigned int vertices_contained = 0;
          for (auto & n : current_elem->node_ref_range())
            if (this_elem->contains_point(n))
              vertices_contained++;

          if (vertices_contained <= this_elem->dim())
            {
              neighbor_set.erase(current);
              continue;
            }
        }
      else
        {
          for (auto & n : this_elem->node_ref_range())
            {
              if (!current_elem->contains_point(n))
                {
                  neighbor_set.erase(current);
                  break;
                }
            }
        }
    }
}

} // namespace ElemInternal
} // namespace libMesh

#endif
