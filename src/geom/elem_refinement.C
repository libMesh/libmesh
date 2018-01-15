// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

void Elem::refine (MeshRefinement & mesh_refinement)
{
  libmesh_assert_equal_to (this->refinement_flag(), Elem::REFINE);
  libmesh_assert (this->active());

  const unsigned int nc = this->n_children();

  // Create my children if necessary
  if (!_children)
    {
      _children = new Elem *[nc];

      unsigned int parent_p_level = this->p_level();
      for (unsigned int c = 0; c != nc; c++)
        {
          _children[c] = Elem::build(this->type(), this).release();
          Elem * current_child = this->child_ptr(c);

          current_child->set_refinement_flag(Elem::JUST_REFINED);
          current_child->set_p_level(parent_p_level);
          current_child->set_p_refinement_flag(this->p_refinement_flag());

          for (unsigned int nc=0; nc<current_child->n_nodes(); nc++)
            {
              Node * node =
                mesh_refinement.add_node(*this, c, nc,
                                         current_child->processor_id());
              node->set_n_systems (this->n_systems());
              current_child->set_node(nc) = node;
            }

          mesh_refinement.add_elem (current_child);
          current_child->set_n_systems(this->n_systems());
        }
    }
  else
    {
      unsigned int parent_p_level = this->p_level();
      for (unsigned int c = 0; c != nc; c++)
        {
          Elem * current_child = this->child_ptr(c);
          libmesh_assert(current_child->subactive());
          current_child->set_refinement_flag(Elem::JUST_REFINED);
          current_child->set_p_level(parent_p_level);
          current_child->set_p_refinement_flag(this->p_refinement_flag());
        }
    }

  // Un-set my refinement flag now
  this->set_refinement_flag(Elem::INACTIVE);

  // Leave the p refinement flag set - we will need that later to get
  // projection operations correct
  // this->set_p_refinement_flag(Elem::INACTIVE);

  for (unsigned int c = 0; c != nc; c++)
    {
      libmesh_assert_equal_to (this->child_ptr(c)->parent(), this);
      libmesh_assert(this->child_ptr(c)->active());
    }
  libmesh_assert (this->ancestor());
}



void Elem::coarsen()
{
  libmesh_assert_equal_to (this->refinement_flag(), Elem::COARSEN_INACTIVE);
  libmesh_assert (!this->active());

  // We no longer delete children until MeshRefinement::contract()
  // delete [] _children;
  // _children = libmesh_nullptr;

  unsigned int parent_p_level = 0;

  // re-compute hanging node nodal locations
  for (unsigned int c = 0, nc = this->n_children(); c != nc; ++c)
    {
      Elem * mychild = this->child_ptr(c);
      if (mychild == remote_elem)
        continue;
      for (unsigned int nc=0; nc != mychild->n_nodes(); nc++)
        {
          Point new_pos;
          bool calculated_new_pos = false;

          for (unsigned int n=0; n<this->n_nodes(); n++)
            {
              // The value from the embedding matrix
              const float em_val = this->embedding_matrix(c,nc,n);

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
                  Point & child_node = mychild->point(nc);
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
  delete [] _children;
  _children = libmesh_nullptr;

  if (this->refinement_flag() == Elem::JUST_COARSENED)
    this->set_refinement_flag(Elem::DO_NOTHING);
}

#endif // #ifdef LIBMESH_ENABLE_AMR


} // namespace libMesh
