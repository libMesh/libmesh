// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

void Elem::refine (MeshRefinement& mesh_refinement)
{
  libmesh_assert_equal_to (this->refinement_flag(), Elem::REFINE);
  libmesh_assert (this->active());

  // Create my children if necessary
  if (!_children)
    {
      _children = new Elem*[this->n_children()];

      unsigned int parent_p_level = this->p_level();
      for (unsigned int c=0; c<this->n_children(); c++)
        {
          _children[c] = Elem::build(this->type(), this).release();
          _children[c]->set_refinement_flag(Elem::JUST_REFINED);
          _children[c]->set_p_level(parent_p_level);
          _children[c]->set_p_refinement_flag(this->p_refinement_flag());
        }

      // Compute new nodal locations
      // and asssign nodes to children
      // Make these static.  It is unlikely the
      // sizes will change from call to call, so having these
      // static should save on reallocations
      std::vector<std::vector<Point> >        p    (this->n_children());
      std::vector<std::vector<Node*> >        nodes(this->n_children());


      // compute new nodal locations
      for (unsigned int c=0; c<this->n_children(); c++)
        {
          Elem *current_child = this->child(c);
          p[c].resize    (current_child->n_nodes());
          nodes[c].resize(current_child->n_nodes());

          for (unsigned int nc=0; nc<current_child->n_nodes(); nc++)
            {
              // zero entries
              p[c][nc].zero();
              nodes[c][nc] = NULL;

              for (unsigned int n=0; n<this->n_nodes(); n++)
                {
                  // The value from the embedding matrix
                  const float em_val = this->embedding_matrix(c,nc,n);

                  if (em_val != 0.)
                    {
                      p[c][nc].add_scaled (this->point(n), em_val);

                      // We may have found the node, in which case we
                      // won't need to look it up later.
                      if (em_val == 1.)
                        nodes[c][nc] = this->get_node(n);
                    }
                }
            }

          // assign nodes to children & add them to the mesh
          const Real pointtol = this->hmin() * TOLERANCE;
          for (unsigned int nc=0; nc<current_child->n_nodes(); nc++)
            {
              if (nodes[c][nc] != NULL)
                {
                  current_child->set_node(nc) = nodes[c][nc];
                }
              else
                {
                  current_child->set_node(nc) =
                    mesh_refinement.add_point(p[c][nc],
                                              current_child->processor_id(),
                                              pointtol);
                  current_child->get_node(nc)->set_n_systems
                    (this->n_systems());
                }
            }

          mesh_refinement.add_elem (current_child);
          current_child->set_n_systems(this->n_systems());
        }
    }
  else
    {
      unsigned int parent_p_level = this->p_level();
      for (unsigned int c=0; c<this->n_children(); c++)
        {
          Elem *current_child = this->child(c);
          libmesh_assert(current_child->subactive());
          current_child->set_refinement_flag(Elem::JUST_REFINED);
          current_child->set_p_level(parent_p_level);
          current_child->set_p_refinement_flag(this->p_refinement_flag());
        }
    }

  // Un-set my refinement flag now
  this->set_refinement_flag(Elem::INACTIVE);
  this->set_p_refinement_flag(Elem::INACTIVE);

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      libmesh_assert_equal_to (this->child(c)->parent(), this);
      libmesh_assert(this->child(c)->active());
    }
  libmesh_assert (this->ancestor());
}



void Elem::coarsen()
{
  libmesh_assert_equal_to (this->refinement_flag(), Elem::COARSEN_INACTIVE);
  libmesh_assert (!this->active());

  // We no longer delete children until MeshRefinement::contract()
  // delete [] _children;
  // _children = NULL;

  unsigned int parent_p_level = 0;

  // re-compute hanging node nodal locations
  for (unsigned int c=0; c<this->n_children(); c++)
    {
      Elem *mychild = this->child(c);
      if (mychild == remote_elem)
        continue;
      for (unsigned int nc=0; nc<mychild->n_nodes(); nc++)
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

          if(calculated_new_pos)
            {
              //Move the existing node back into it's original location
              for(unsigned int i=0; i<LIBMESH_DIM; i++)
                {
                  Point & child_node = *(mychild->get_node(nc));
                  child_node(i)=new_pos(i);
                }
            }
        }
    }

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      Elem *mychild = this->child(c);
      if (mychild == remote_elem)
        continue;
      libmesh_assert_equal_to (mychild->refinement_flag(), Elem::COARSEN);
      mychild->set_refinement_flag(Elem::INACTIVE);
      if (mychild->p_level() > parent_p_level)
        parent_p_level = mychild->p_level();
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
  _children = NULL;

  if (this->refinement_flag() == Elem::JUST_COARSENED)
    this->set_refinement_flag(Elem::DO_NOTHING);
}

#endif // #ifdef LIBMESH_ENABLE_AMR


} // namespace libMesh
