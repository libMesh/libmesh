// $Id: elem_refinement.C,v 1.12 2004-10-28 20:06:14 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
  assert (_children == NULL);
  
  // Two big prime numbers less than
  // sqrt(max_unsigned_int) for key creation.
  const unsigned int bp1 = 65449;
  const unsigned int bp2 = 48661;
  
  
  // Create my children
  {
    _children = new Elem*[this->n_children()];

    for (unsigned int c=0; c<this->n_children(); c++)
      {
	_children[c] = Elem::build(this->type(), this).release();
	_children[c]->set_refinement_flag(Elem::JUST_REFINED);
      }
  }


  // Compute new nodal locations
  // and asssign nodes to children
  {
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
			// An unsigned int associated with the
			// address of the node n.  We can't use the
			// node number since they can change.
			
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
			WHAT KIND OF CRAZY MACHINE IS THIS? CANNOT COMPILE
			
#endif
			  // Compute the key for this new node nc.  This will
			  // be used to locate the node if it already exists
			  // in the mesh.
			  keys[c][nc] +=
			  (((static_cast<unsigned int>(em_val*100000.)%bp1) *
			    (n_id%bp1))%bp1)*bp2;
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
	      }
	  }
      
	mesh_refinement.add_elem (this->child(c));
      }
  }

  // Un-set my refinement flag now
  this->set_refinement_flag(Elem::DO_NOTHING);

  assert (!this->active());
}



void Elem::coarsen()
{
  assert (this->refinement_flag() == Elem::COARSEN);
  assert (!this->active());

  // Delete the storage for my children
  delete [] _children;

  _children = NULL;

  this->set_refinement_flag(Elem::DO_NOTHING);

  assert (this->active());
}

#endif // #ifdef ENABLE_AMR


