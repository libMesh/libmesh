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



// C++ includes

// Local includes
#include "libmesh/elem.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/remote_elem.h"

namespace libMesh
{


//--------------------------------------------------------------------
// Elem methods

/**
 * The following functions only apply when
 * AMR is enabled and thus are not present
 * otherwise.
 */
#ifdef LIBMESH_ENABLE_AMR

void Elem::set_p_level(unsigned int p)
{
  // Maintain the parent's p level as the minimum of it's children
  if (this->parent() != nullptr)
    {
      unsigned int parent_p_level = this->parent()->p_level();

      // If our new p level is less than our parents, our parents drops
      if (parent_p_level > p)
        {
          this->parent()->set_p_level(p);

          // And we should keep track of the drop, in case we need to
          // do a projection later.
          this->parent()->set_p_refinement_flag(Elem::JUST_COARSENED);
        }
      // If we are the lowest p level and it increases, so might
      // our parent's, but we have to check every other child to see
      else if (parent_p_level == _p_level && _p_level < p)
        {
          _p_level = cast_int<unsigned char>(p);
          parent_p_level = cast_int<unsigned char>(p);
          for (auto & c : this->parent()->child_ref_range())
            parent_p_level = std::min(parent_p_level,
                                      c.p_level());

          // When its children all have a higher p level, the parent's
          // should rise
          if (parent_p_level > this->parent()->p_level())
            {
              this->parent()->set_p_level(parent_p_level);

              // And we should keep track of the rise, in case we need to
              // do a projection later.
              this->parent()->set_p_refinement_flag(Elem::JUST_REFINED);
            }

          return;
        }
    }

  this->hack_p_level(p);
}



void Elem::refine (MeshRefinement & mesh_refinement)
{
  libmesh_assert_equal_to (this->refinement_flag(), Elem::REFINE);
  libmesh_assert (this->active());

  const unsigned int nc = this->n_children();

  // Create my children if necessary
  if (!_children)
    {
      _children = std::make_unique<Elem *[]>(nc);

      unsigned int parent_p_level = this->p_level();
      const unsigned int nei = this->n_extra_integers();
      for (unsigned int c = 0; c != nc; c++)
        {
          auto current_child = Elem::build(this->type(), this);
          _children[c] = current_child.get();

          current_child->set_refinement_flag(Elem::JUST_REFINED);
          current_child->set_p_level(parent_p_level);
          current_child->set_p_refinement_flag(this->p_refinement_flag());

          for (auto cnode : current_child->node_index_range())
            {
              Node * node =
                mesh_refinement.add_node(*this, c, cnode,
                                         current_child->processor_id());
              node->set_n_systems (this->n_systems());
              current_child->set_node(cnode, node);
            }

          Elem * added_child = mesh_refinement.add_elem (std::move(current_child));
          added_child->set_n_systems(this->n_systems());
          libmesh_assert_equal_to (added_child->n_extra_integers(),
                                   this->n_extra_integers());
          for (unsigned int i=0; i != nei; ++i)
            added_child->set_extra_integer(i, this->get_extra_integer(i));
        }
    }
  else
    {
      unsigned int parent_p_level = this->p_level();
      for (unsigned int c = 0; c != nc; c++)
        {
          Elem * current_child = this->child_ptr(c);
          if (current_child != remote_elem)
            {
              libmesh_assert(current_child->subactive());
              current_child->set_refinement_flag(Elem::JUST_REFINED);
              current_child->set_p_level(parent_p_level);
              current_child->set_p_refinement_flag(this->p_refinement_flag());
            }
        }
    }

  // Get any new remote neighbor links right; even find_neighbors
  // relies on those
  for (unsigned int s : this->side_index_range())
    {
      if (this->neighbor_ptr(s) != remote_elem)
        continue;

      for (unsigned int c = 0; c != nc; c++)
        {
          Elem * current_child = this->child_ptr(c);
          if (current_child != remote_elem &&
              this->is_child_on_side(c, s))
            current_child->set_neighbor
              (s, const_cast<RemoteElem *>(remote_elem));
        }
    }

  // Un-set my refinement flag now
  this->set_refinement_flag(Elem::INACTIVE);

  // Leave the p refinement flag set - we will need that later to get
  // projection operations correct
  // this->set_p_refinement_flag(Elem::INACTIVE);

#ifndef NDEBUG
  for (unsigned int c = 0; c != nc; c++)
    if (this->child_ptr(c) != remote_elem)
      {
        libmesh_assert_equal_to (this->child_ptr(c)->parent(), this);
        libmesh_assert(this->child_ptr(c)->active());
      }
#endif
  libmesh_assert (this->ancestor());
}



void Elem::coarsen()
{
  libmesh_assert_equal_to (this->refinement_flag(), Elem::COARSEN_INACTIVE);
  libmesh_assert (!this->active());

  // We no longer delete children until MeshRefinement::contract()

  unsigned int parent_p_level = 0;

  const unsigned int n_n = this->n_nodes();

  // re-compute hanging node nodal locations
  for (unsigned int c = 0, nc = this->n_children(); c != nc; ++c)
    {
      Elem * mychild = this->child_ptr(c);
      if (mychild == remote_elem)
        continue;
      for (auto cnode : mychild->node_index_range())
        {
          Point new_pos;
          bool calculated_new_pos = false;

          for (unsigned int n=0; n<n_n; n++)
            {
              // The value from the embedding matrix
              const Real em_val = this->embedding_matrix(c,cnode,n);

              // The node location is somewhere between existing vertices
              if ((em_val != 0.) && (em_val != 1.))
                {
                  new_pos.add_scaled (this->point(n), em_val);
                  calculated_new_pos = true;
                }
            }

          if (calculated_new_pos)
            {
              //Move the existing node back into it's original location
              for (unsigned int i=0; i<LIBMESH_DIM; i++)
                {
                  Point & child_node = mychild->point(cnode);
                  child_node(i)=new_pos(i);
                }
            }
        }
    }

  for (auto & mychild : this->child_ref_range())
    {
      if (&mychild == remote_elem)
        continue;
      libmesh_assert_equal_to (mychild.refinement_flag(), Elem::COARSEN);
      mychild.set_refinement_flag(Elem::INACTIVE);
      if (mychild.p_level() > parent_p_level)
        parent_p_level = mychild.p_level();
    }

  this->set_refinement_flag(Elem::JUST_COARSENED);
  this->set_p_level(parent_p_level);

  libmesh_assert (this->active());
}



void Elem::contract()
{
  // Subactive elements get deleted entirely, not contracted
  libmesh_assert (this->active());

  // Active contracted elements no longer can have children
  _children.reset(nullptr);

  if (this->refinement_flag() == Elem::JUST_COARSENED)
    this->set_refinement_flag(Elem::DO_NOTHING);
}

#endif // #ifdef LIBMESH_ENABLE_AMR


} // namespace libMesh
