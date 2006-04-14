// $Id: elem_refinement.C,v 1.20 2006-04-14 22:31:50 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
  assert (this->refinement_flag() == Elem::REFINE);
  assert (this->active());
  
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
      std::vector<std::vector<unsigned int> > keys (this->n_children());
      std::vector<std::vector<Node*> >        nodes(this->n_children());
    

      // compute new nodal locations
      for (unsigned int c=0; c<this->n_children(); c++)
        {	
	  p[c].resize    (this->child(c)->n_nodes());
	  keys[c].resize (this->child(c)->n_nodes());
	  nodes[c].resize(this->child(c)->n_nodes());

	  for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
	    {
	      // zero entries
	      p[c][nc].zero();
	      keys[c][nc]  = 0;
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
		  
		      // Otherwise build the key to look for the node
		      else
		        {
                          const unsigned int n_id = _cast_node_address_to_unsigned_int(n);
			  
                          // Compute the key for this new node nc.  This will
			  // be used to locate the node if it already exists
			  // in the mesh.
			  keys[c][nc] +=
			  (((static_cast<unsigned int>(em_val*100000.)%_bp1) *
			    (n_id%_bp1))%_bp1)*_bp2;
		        }
		    }
	        }
	    }
      
	// assign nodes to children & add them to the mesh
	  for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
	    {
	      if (nodes[c][nc] != NULL)
	        {
		  this->child(c)->set_node(nc) = nodes[c][nc];
	        }
	      else
	        {
		  this->child(c)->set_node(nc) =
		    mesh_refinement.add_point(p[c][nc],
					      keys[c][nc]);
		  this->child(c)->get_node(nc)->set_n_systems
                    (this->n_systems());
	        }
	    }
      
	  mesh_refinement.add_elem (this->child(c));
          this->child(c)->set_n_systems(this->n_systems());
        }
    }
  else
    {
      unsigned int parent_p_level = this->p_level();
      for (unsigned int c=0; c<this->n_children(); c++)
        {	
          assert(this->child(c)->subactive());
          this->child(c)->set_refinement_flag(Elem::JUST_REFINED);
          this->child(c)->set_p_level(parent_p_level);
          this->child(c)->set_p_refinement_flag(this->p_refinement_flag());
        }
    }

  // Un-set my refinement flag now
  this->set_refinement_flag(Elem::INACTIVE);
  this->set_p_refinement_flag(Elem::INACTIVE);

  for (unsigned int c=0; c<this->n_children(); c++)
    {	
      assert(this->child(c)->parent() == this);
      assert(this->child(c)->active());
    }
  assert (this->ancestor());
}




void Elem::compute_children_node_keys()
{
  // No node keys need to be set if I have no children 
  if (!_children)
    return;

  // Temporary storage for keys. 
  std::vector<std::vector<unsigned int> > keys (this->n_children());


  // compute new nodal locations
  for (unsigned int c=0; c<this->n_children(); c++)
  {	
    keys[c].resize (this->child(c)->n_nodes());

    for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
    {
      for (unsigned int n=0; n<this->n_nodes(); n++)
      {
        // The value from the embedding matrix
        const float em_val = this->embedding_matrix(c,nc,n);
         
        // The node location is somewhere between existing vertices
        if ((em_val != 0.) && (em_val != 1.)) 
        {
          const unsigned int n_id = _cast_node_address_to_unsigned_int(n);

          // Compute the key for this new node nc.  This will
          // be used to locate the node if it already exists
          // in the mesh.
          keys[c][nc] +=
            (((static_cast<unsigned int>(em_val*100000.)%_bp1) *
              (n_id%_bp1))%_bp1)*_bp2;
        }
      }
    }
  }

  // By now, all the keys for the child nodes have been computed.
  // keys[c][nc] will have zero entries for nodes that already belonged
  // to the parent.
  for (unsigned int c=0; c<this->n_children(); c++)
    for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
      if (keys[c][nc] != 0)
        this->child(c)->get_node(nc)->set_key() = keys[c][nc];
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
  assert (this->refinement_flag() == Elem::COARSEN_INACTIVE);
  assert (!this->active());

  // We no longer delete children until MeshRefinement::contract()
  // delete [] _children;
  // _children = NULL;

  unsigned int parent_p_level = 0;
  for (unsigned int c=0; c<this->n_children(); c++)
    {	
      assert (this->child(c)->refinement_flag() == Elem::COARSEN);
      this->child(c)->set_refinement_flag(Elem::INACTIVE);
      if (this->child(c)->p_level() > parent_p_level)
        parent_p_level = this->child(c)->p_level();
    }

  this->set_refinement_flag(Elem::JUST_COARSENED);
  this->set_p_level(parent_p_level);

  assert (this->active());
}



void Elem::contract()
{
  // Subactive elements get deleted entirely, not contracted
  assert (this->active());

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


