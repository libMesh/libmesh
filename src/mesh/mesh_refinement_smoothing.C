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



// Local includes
#include "libmesh/libmesh_config.h"

// only compile these functions if the user requests AMR support
#ifdef LIBMESH_ENABLE_AMR

#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/parallel.h"
#include "libmesh/remote_elem.h"

namespace libMesh
{



//-----------------------------------------------------------------
// Mesh refinement methods
bool MeshRefinement::limit_level_mismatch_at_node (const unsigned int max_mismatch)
{
  // This function must be run on all processors at once
  parallel_object_only();

  bool flags_changed = false;


  // Vector holding the maximum element level that touches a node.
  std::vector<unsigned char> max_level_at_node (_mesh.n_nodes(), 0);
  std::vector<unsigned char> max_p_level_at_node (_mesh.n_nodes(), 0);

  // Loop over all the active elements & fill the vector
  for (auto & elem : _mesh.active_element_ptr_range())
    {
      const unsigned char elem_level =
        cast_int<unsigned char>(elem->level() +
                                ((elem->refinement_flag() == Elem::REFINE) ? 1 : 0));
      const unsigned char elem_p_level =
        cast_int<unsigned char>(elem->p_level() +
                                ((elem->p_refinement_flag() == Elem::REFINE) ? 1 : 0));

      // Set the max_level at each node
      for (unsigned int n=0; n<elem->n_nodes(); n++)
        {
          const dof_id_type node_number = elem->node_id(n);

          libmesh_assert_less (node_number, max_level_at_node.size());

          max_level_at_node[node_number] =
            std::max (max_level_at_node[node_number], elem_level);
          max_p_level_at_node[node_number] =
            std::max (max_p_level_at_node[node_number], elem_p_level);
        }
    }


  // Now loop over the active elements and flag the elements
  // who violate the requested level mismatch. Alternatively, if
  // _enforce_mismatch_limit_prior_to_refinement is true, swap refinement flags
  // accordingly.
  for (auto & elem : _mesh.active_element_ptr_range())
    {
      const unsigned int elem_level = elem->level();
      const unsigned int elem_p_level = elem->p_level();

      // Skip the element if it is already fully flagged
      // unless we are enforcing mismatch prior to refinement and may need to
      // remove the refinement flag(s)
      if (elem->refinement_flag() == Elem::REFINE &&
          elem->p_refinement_flag() == Elem::REFINE
          && !_enforce_mismatch_limit_prior_to_refinement)
        continue;

      // Loop over the nodes, check for possible mismatch
      for (unsigned int n=0; n<elem->n_nodes(); n++)
        {
          const dof_id_type node_number = elem->node_id(n);

          // Flag the element for refinement if it violates
          // the requested level mismatch
          if ((elem_level + max_mismatch) < max_level_at_node[node_number]
              && elem->refinement_flag() != Elem::REFINE)
            {
              elem->set_refinement_flag (Elem::REFINE);
              flags_changed = true;
            }
          if ((elem_p_level + max_mismatch) < max_p_level_at_node[node_number]
              && elem->p_refinement_flag() != Elem::REFINE)
            {
              elem->set_p_refinement_flag (Elem::REFINE);
              flags_changed = true;
            }

          // Possibly enforce limit mismatch prior to refinement
          flags_changed |= this->enforce_mismatch_limit_prior_to_refinement(elem, POINT, max_mismatch);
        }
    }

  // If flags changed on any processor then they changed globally
  this->comm().max(flags_changed);

  return flags_changed;
}



bool MeshRefinement::limit_level_mismatch_at_edge (const unsigned int max_mismatch)
{
  // This function must be run on all processors at once
  parallel_object_only();

  bool flags_changed = false;


  // Maps holding the maximum element level that touches an edge
  std::map<std::pair<unsigned int, unsigned int>, unsigned char>
    max_level_at_edge;
  std::map<std::pair<unsigned int, unsigned int>, unsigned char>
    max_p_level_at_edge;

  // Loop over all the active elements & fill the maps
  for (auto & elem : _mesh.active_element_ptr_range())
    {
      const unsigned char elem_level =
        cast_int<unsigned char>(elem->level() +
                                ((elem->refinement_flag() == Elem::REFINE) ? 1 : 0));
      const unsigned char elem_p_level =
        cast_int<unsigned char>(elem->p_level() +
                                ((elem->p_refinement_flag() == Elem::REFINE) ? 1 : 0));

      // Set the max_level at each edge
      for (unsigned int n=0; n<elem->n_edges(); n++)
        {
          std::unique_ptr<const Elem> edge = elem->build_edge_ptr(n);
          dof_id_type childnode0 = edge->node_id(0);
          dof_id_type childnode1 = edge->node_id(1);
          if (childnode1 < childnode0)
            std::swap(childnode0, childnode1);

          for (const Elem * p = elem; p != nullptr; p = p->parent())
            {
              std::unique_ptr<const Elem> pedge = p->build_edge_ptr(n);
              dof_id_type node0 = pedge->node_id(0);
              dof_id_type node1 = pedge->node_id(1);

              if (node1 < node0)
                std::swap(node0, node1);

              // If elem does not share this edge with its ancestor
              // p, refinement levels of elements sharing p's edge
              // are not restricted by refinement levels of elem.
              // Furthermore, elem will not share this edge with any
              // of p's ancestors, so we can safely break out of the
              // for loop early.
              if (node0 != childnode0 && node1 != childnode1)
                break;

              childnode0 = node0;
              childnode1 = node1;

              std::pair<unsigned int, unsigned int> edge_key =
                std::make_pair(node0, node1);

              if (max_level_at_edge.find(edge_key) ==
                  max_level_at_edge.end())
                {
                  max_level_at_edge[edge_key] = elem_level;
                  max_p_level_at_edge[edge_key] = elem_p_level;
                }
              else
                {
                  max_level_at_edge[edge_key] =
                    std::max (max_level_at_edge[edge_key], elem_level);
                  max_p_level_at_edge[edge_key] =
                    std::max (max_p_level_at_edge[edge_key], elem_p_level);
                }
            }
        }
    }


  // Now loop over the active elements and flag the elements
  // who violate the requested level mismatch
  for (auto & elem : _mesh.active_element_ptr_range())
    {
      const unsigned int elem_level = elem->level();
      const unsigned int elem_p_level = elem->p_level();

      // Skip the element if it is already fully flagged
      if (elem->refinement_flag() == Elem::REFINE &&
          elem->p_refinement_flag() == Elem::REFINE
          && !_enforce_mismatch_limit_prior_to_refinement)
        continue;

      // Loop over the nodes, check for possible mismatch
      for (unsigned int n=0; n<elem->n_edges(); n++)
        {
          std::unique_ptr<Elem> edge = elem->build_edge_ptr(n);
          dof_id_type node0 = edge->node_id(0);
          dof_id_type node1 = edge->node_id(1);
          if (node1 < node0)
            std::swap(node0, node1);

          std::pair<dof_id_type, dof_id_type> edge_key =
            std::make_pair(node0, node1);

          // Flag the element for refinement if it violates
          // the requested level mismatch
          if ((elem_level + max_mismatch) < max_level_at_edge[edge_key]
              && elem->refinement_flag() != Elem::REFINE)
            {
              elem->set_refinement_flag (Elem::REFINE);
              flags_changed = true;
            }

          if ((elem_p_level + max_mismatch) < max_p_level_at_edge[edge_key]
              && elem->p_refinement_flag() != Elem::REFINE)
            {
              elem->set_p_refinement_flag (Elem::REFINE);
              flags_changed = true;
            }

          // Possibly enforce limit mismatch prior to refinement
          flags_changed |= this->enforce_mismatch_limit_prior_to_refinement(elem, EDGE, max_mismatch);
        } // loop over edges
    } // loop over active elements

  // If flags changed on any processor then they changed globally
  this->comm().max(flags_changed);

  return flags_changed;
}



bool MeshRefinement::limit_overrefined_boundary(const signed char max_mismatch)
{
  // This function must be run on all processors at once
  parallel_object_only();

  bool flags_changed = false;

  // Loop over all the active elements & look for mismatches to fix.
  for (auto & elem : _mesh.active_element_ptr_range())
    {
      // If we don't have an interior_parent then there's nothing to
      // be mismatched with.
      if ((elem->dim() >= LIBMESH_DIM) ||
          !elem->interior_parent())
        continue;

      const unsigned char elem_level =
        cast_int<unsigned char>(elem->level() +
                                ((elem->refinement_flag() == Elem::REFINE) ? 1 : 0));
      const unsigned char elem_p_level =
        cast_int<unsigned char>(elem->p_level() +
                                ((elem->p_refinement_flag() == Elem::REFINE) ? 1 : 0));

      // get all relevant interior elements
      std::set<Elem *> neighbor_set;
      elem->find_interior_neighbors(neighbor_set);

      for (auto & neighbor : neighbor_set)
        if (max_mismatch >= 0)
          {
            if ((elem_level > neighbor->level() + max_mismatch) &&
                (neighbor->refinement_flag() != Elem::REFINE))
              {
                neighbor->set_refinement_flag(Elem::REFINE);
                flags_changed = true;
              }

            if ((elem_p_level > neighbor->p_level() + max_mismatch) &&
                (neighbor->p_refinement_flag() != Elem::REFINE))
              {
                neighbor->set_p_refinement_flag(Elem::REFINE);
                flags_changed = true;
              }
          }
    }

  // If flags changed on any processor then they changed globally
  this->comm().max(flags_changed);

  return flags_changed;
}



bool MeshRefinement::limit_underrefined_boundary(const signed char max_mismatch)
{
  // This function must be run on all processors at once
  parallel_object_only();

  bool flags_changed = false;

  // Loop over all the active elements & look for mismatches to fix.
  for (auto & elem : _mesh.active_element_ptr_range())
    {
      // If we don't have an interior_parent then there's nothing to
      // be mismatched with.
      if ((elem->dim() >= LIBMESH_DIM) ||
          !elem->interior_parent())
        continue;

      // get all relevant interior elements
      std::set<const Elem *> neighbor_set;
      elem->find_interior_neighbors(neighbor_set);

      std::set<const Elem *>::iterator n_it = neighbor_set.begin();
      for (; n_it != neighbor_set.end(); ++n_it)
        {
          // FIXME - non-const versions of the Elem set methods
          // would be nice
          const Elem * neighbor = *n_it;

          const unsigned char neighbor_level =
            cast_int<unsigned char>(neighbor->level() +
                                    ((neighbor->refinement_flag() == Elem::REFINE) ? 1 : 0));

          const unsigned char neighbor_p_level =
            cast_int<unsigned char>(neighbor->p_level() +
                                    ((neighbor->p_refinement_flag() == Elem::REFINE) ? 1 : 0));

          if (max_mismatch >= 0)
            {
              if ((neighbor_level >
                   elem->level() + max_mismatch) &&
                  (elem->refinement_flag() != Elem::REFINE))
                {
                  elem->set_refinement_flag(Elem::REFINE);
                  flags_changed = true;
                }

              if ((neighbor_p_level >
                   elem->p_level() + max_mismatch) &&
                  (elem->p_refinement_flag() != Elem::REFINE))
                {
                  elem->set_p_refinement_flag(Elem::REFINE);
                  flags_changed = true;
                }
            }
        } // loop over interior neighbors
    }

  // If flags changed on any processor then they changed globally
  this->comm().max(flags_changed);

  return flags_changed;
}



bool MeshRefinement::eliminate_unrefined_patches ()
{
  // This function must be run on all processors at once
  parallel_object_only();

  bool flags_changed = false;

  // Note: we *cannot* use a reference to the real pointer here, since
  // the pointer may be reseated below and we don't want to reseat
  // pointers held by the Mesh.
  for (Elem * elem : _mesh.active_element_ptr_range())
    {
      // First, see if there's any possibility we might have to flag
      // this element for h and p refinement - do we have any visible
      // neighbors?  Next we'll check to see if any of those neighbors
      // are as coarse or coarser than us.
      bool h_flag_me = false,
           p_flag_me = false;
      for (auto neighbor : elem->neighbor_ptr_range())
        {
          // Quit if the element is not a local boundary
          if (neighbor != nullptr && neighbor != remote_elem)
            {
              h_flag_me = true;
              p_flag_me = true;
              break;
            }
        }

      // Skip the element if it is already fully flagged for refinement
      if (elem->p_refinement_flag() == Elem::REFINE)
        p_flag_me = false;
      if (elem->refinement_flag() == Elem::REFINE)
        {
          h_flag_me = false;
          if (!p_flag_me)
            continue;
        }
      // Test the parent if that is already flagged for coarsening
      else if (elem->refinement_flag() == Elem::COARSEN)
        {
          libmesh_assert(elem->parent());
          elem = elem->parent();
          // FIXME - this doesn't seem right - RHS
          if (elem->refinement_flag() != Elem::COARSEN_INACTIVE)
            continue;
          p_flag_me = false;
        }

      const unsigned int my_level = elem->level();
      int my_p_adjustment = 0;
      if (elem->p_refinement_flag() == Elem::REFINE)
        my_p_adjustment = 1;
      else if (elem->p_refinement_flag() == Elem::COARSEN)
        {
          libmesh_assert_greater (elem->p_level(), 0);
          my_p_adjustment = -1;
        }
      const unsigned int my_new_p_level = elem->p_level() +
        my_p_adjustment;

      // Check all the element neighbors
      for (auto neighbor : elem->neighbor_ptr_range())
        {
          // Quit if the element is on a local boundary
          if (neighbor == nullptr || neighbor == remote_elem)
            {
              h_flag_me = false;
              p_flag_me = false;
              break;
            }
          // if the neighbor will be equally or less refined than
          // we are, then we will not become an unrefined island.
          // So if we are still considering h refinement:
          if (h_flag_me &&
              // If our neighbor is already at a lower level,
              // it can't end up at a higher level even if it
              // is flagged for refinement once
              ((neighbor->level() < my_level) ||
               // If our neighbor is at the same level but isn't
               // flagged for refinement, it won't end up at a
               // higher level
               ((neighbor->active()) &&
                (neighbor->refinement_flag() != Elem::REFINE)) ||
               // If our neighbor is currently more refined but is
               // a parent flagged for coarsening, it will end up
               // at the same level.
               (neighbor->refinement_flag() == Elem::COARSEN_INACTIVE)))
            {
              // We've proven we won't become an unrefined island,
              // so don't h refine to avoid that.
              h_flag_me = false;

              // If we've also proven we don't need to p refine,
              // we don't need to check more neighbors
              if (!p_flag_me)
                break;
            }
          if (p_flag_me)
            {
              // if active neighbors will have a p level
              // equal to or lower than ours, then we do not need to p
              // refine ourselves.
              if (neighbor->active())
                {
                  int p_adjustment = 0;
                  if (neighbor->p_refinement_flag() == Elem::REFINE)
                    p_adjustment = 1;
                  else if (neighbor->p_refinement_flag() == Elem::COARSEN)
                    {
                      libmesh_assert_greater (neighbor->p_level(), 0);
                      p_adjustment = -1;
                    }
                  if (my_new_p_level >= neighbor->p_level() + p_adjustment)
                    {
                      p_flag_me = false;
                      if (!h_flag_me)
                        break;
                    }
                }
              // If we have inactive neighbors, we need to
              // test all their active descendants which neighbor us
              else if (neighbor->ancestor())
                {
                  if (neighbor->min_new_p_level_by_neighbor(elem,
                                                            my_new_p_level + 2) <= my_new_p_level)
                    {
                      p_flag_me = false;
                      if (!h_flag_me)
                        break;
                    }
                }
            }
        }

      if (h_flag_me)
        {
          // Parents that would create islands should no longer
          // coarsen
          if (elem->refinement_flag() == Elem::COARSEN_INACTIVE)
            {
              for (auto & child : elem->child_ref_range())
                {
                  libmesh_assert_equal_to (child.refinement_flag(),
                                           Elem::COARSEN);
                  child.set_refinement_flag(Elem::DO_NOTHING);
                }
              elem->set_refinement_flag(Elem::INACTIVE);
            }
          else
            elem->set_refinement_flag(Elem::REFINE);
          flags_changed = true;
        }
      if (p_flag_me)
        {
          if (elem->p_refinement_flag() == Elem::COARSEN)
            elem->set_p_refinement_flag(Elem::DO_NOTHING);
          else
            elem->set_p_refinement_flag(Elem::REFINE);
          flags_changed = true;
        }
    }

  // If flags changed on any processor then they changed globally
  this->comm().max(flags_changed);

  return flags_changed;
}



bool MeshRefinement::enforce_mismatch_limit_prior_to_refinement(Elem * elem,
                                                                NeighborType nt,
                                                                unsigned max_mismatch)
{
  // Eventual return value
  bool flags_changed = false;

  // If we are enforcing the limit prior to refinement then we
  // need to remove flags from any elements marked for refinement that
  // would cause a mismatch
  if (_enforce_mismatch_limit_prior_to_refinement
      && elem->refinement_flag() == Elem::REFINE)
    {
      // get all the relevant neighbors since we may have to refine
      // elements off edges or corners as well
      std::set<const Elem *> neighbor_set;

      if (nt == POINT)
        elem->find_point_neighbors(neighbor_set);
      else if (nt == EDGE)
        elem->find_edge_neighbors(neighbor_set);
      else
        libmesh_error_msg("Unrecognized NeighborType: " << nt);

      // Loop over the neighbors of element e
      for (const auto & neighbor : neighbor_set)
        {
          if ((elem->level() + 1 - max_mismatch) > neighbor->level())
            {
              elem->set_refinement_flag(Elem::DO_NOTHING);
              flags_changed = true;
            }
          if ((elem->p_level() + 1 - max_mismatch) > neighbor->p_level())
            {
              elem->set_p_refinement_flag(Elem::DO_NOTHING);
              flags_changed = true;
            }
        } // loop over edge/point neighbors
    } // if _enforce_mismatch_limit_prior_to_refinement

  return flags_changed;
}

} // namespace libMesh


#endif
