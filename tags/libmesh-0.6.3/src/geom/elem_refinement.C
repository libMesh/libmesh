// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include "elem.h"
#include "mesh_refinement.h"
#include "remote_elem.h"


//--------------------------------------------------------------------
// Elem methods

/**
 * The following functions only apply when
 * AMR is enabled and thus are not present
 * otherwise.
 */ 
#ifdef ENABLE_AMR

void Elem::refine (MeshRefinement& mesh_refinement)
{
  libmesh_assert (this->refinement_flag() == Elem::REFINE);
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
          Elem *child = this->child(c);
	  p[c].resize    (child->n_nodes());
	  nodes[c].resize(child->n_nodes());

	  for (unsigned int nc=0; nc<child->n_nodes(); nc++)
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
	  for (unsigned int nc=0; nc<child->n_nodes(); nc++)
	    {
	      if (nodes[c][nc] != NULL)
	        {
		  child->set_node(nc) = nodes[c][nc];
	        }
	      else
	        {
		  child->set_node(nc) =
		    mesh_refinement.add_point(p[c][nc],
					      child->processor_id(),
                                              pointtol);
		  child->get_node(nc)->set_n_systems
                    (this->n_systems());
	        }
	    }
      
	  mesh_refinement.add_elem (child);
          child->set_n_systems(this->n_systems());
        }
    }
  else
    {
      unsigned int parent_p_level = this->p_level();
      for (unsigned int c=0; c<this->n_children(); c++)
        {	
          Elem *child = this->child(c);
          libmesh_assert(child->subactive());
          child->set_refinement_flag(Elem::JUST_REFINED);
          child->set_p_level(parent_p_level);
          child->set_p_refinement_flag(this->p_refinement_flag());
        }
    }

  // Un-set my refinement flag now
  this->set_refinement_flag(Elem::INACTIVE);
  this->set_p_refinement_flag(Elem::INACTIVE);

  for (unsigned int c=0; c<this->n_children(); c++)
    {	
      libmesh_assert(this->child(c)->parent() == this);
      libmesh_assert(this->child(c)->active());
    }
  libmesh_assert (this->ancestor());
}



unsigned int Elem::_cast_node_address_to_unsigned_int(const unsigned int n)
{
  // An unsigned int associated with the
  // address of the node n.  We can't use the
  // node number since they can change, so we use the
  // Node's address.  (We also can't use the x,y,z
  // location of the node since that can change too!)

#if SIZEOF_INT == SIZEOF_VOID_P

  // 32-bit machines
  const unsigned int n_id =
    reinterpret_cast<unsigned int>(this->get_node(n));

#elif SIZEOF_LONG_INT == SIZEOF_VOID_P

  // 64-bit machines 
  // Another big prime number less than max_unsigned_int
  // for key creation on 64-bit machines
  const unsigned int bp3 = 4294967291;
  const unsigned int n_id =			
    reinterpret_cast<long unsigned int>(this->get_node(n))%bp3;

#else
  // Huh?
#error			WHAT KIND OF CRAZY MACHINE IS THIS? CANNOT COMPILE

#endif

  return n_id;
}





void Elem::coarsen()
{
  libmesh_assert (this->refinement_flag() == Elem::COARSEN_INACTIVE);
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
	for(unsigned int i=0; i<DIM; i++)
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
      libmesh_assert (mychild->refinement_flag() == Elem::COARSEN);
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
  if (_children)
    {
      delete [] _children;
      _children = NULL;
    }
  if (this->refinement_flag() == Elem::JUST_COARSENED)
    this->set_refinement_flag(Elem::DO_NOTHING);
}

#endif // #ifdef ENABLE_AMR


