// $Id: dof_map.C,v 1.47 2003-07-18 18:53:00 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



// C++ Includes -------------------------------------
#include <set>
#include <algorithm>
#include <math.h>


// Local Includes -----------------------------------
#include "dof_map.h"
#include "elem.h"
#include "mesh_base.h"
#include "fe_interface.h"
#include "sparse_matrix.h"
#include "libmesh.h"
#include "mesh_logging.h"




// ------------------------------------------------------------
// DofMap member functions
void DofMap::attach_matrix (SparseMatrix<Number>& matrix)
{
  _matrix = &matrix;
  
  _matrix->attach_dof_map (*this);
}



void DofMap::attach_other_matrix (std::pair<std::string, SparseMatrix<Number>*> _new)
{
  if (_other_matrices.count(_new.first))
    {
      std::cerr << "ERROR: matrix "
		<< _new.first
		<< " has already been added to this system!"
		<< std::endl;
      error();
    }
  else
    {
      _other_matrices.insert(_new);
      _new.second->attach_dof_map(*this);
    }
}


void DofMap::reinit(const MeshBase& mesh)
{
  assert (mesh.is_prepared());
  
  START_LOG("reinit()", "DofMap");
  
  //this->clear();

  const unsigned int n_var = this->n_variables();
  const unsigned int dim   = mesh.mesh_dimension();


#ifdef ENABLE_AMR
  
  //------------------------------------------------------------
  // Clear the old_dof_objects for all the nodes
  // and elements so that we can overwrite them
  {
    const_node_iterator       node_it (mesh.nodes_begin());
    const const_node_iterator node_end(mesh.nodes_end());
    
    for ( ; node_it != node_end; ++node_it)
      {
	(*node_it)->clear_old_dof_object();
	assert ((*node_it)->old_dof_object == NULL);
      }
    
    const_elem_iterator       elem_it (mesh.elements_begin());
    const const_elem_iterator elem_end(mesh.elements_end());
    
    for ( ; elem_it != elem_end; ++elem_it)
      {
	(*elem_it)->clear_old_dof_object();
	assert ((*elem_it)->old_dof_object == NULL);
      }
  }

  
  //------------------------------------------------------------
  // Set the old_dof_objects for the elements that
  // weren't just created
  {
    const_elem_iterator       elem_it (mesh.elements_begin());
    const const_elem_iterator elem_end(mesh.elements_end());
    
    for ( ; elem_it != elem_end; ++elem_it)
      {
	Elem* elem = *elem_it;

        // Skip the elements that were just refined
	if (elem->refinement_flag() == Elem::JUST_REFINED) continue;

	for (unsigned int n=0; n<elem->n_nodes(); n++)
	  {
	    Node* node = elem->get_node(n);

	    if (node->old_dof_object == NULL)
	      node->set_old_dof_object();
	  }

	assert (elem->old_dof_object == NULL);
	elem->set_old_dof_object();
      }
  }

#endif // #ifdef ENABLE_AMR

  
  //------------------------------------------------------------
  // Then set the number of variables for each \p DofObject
  // equal to n_variables() for this system.  This will
  // handle new \p DofObjects that may have just been created
  {
    // All the nodes
    const_node_iterator       node_it  (mesh.nodes_begin());
    const const_node_iterator node_end (mesh.nodes_end());

    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_vars(this->sys_number(),n_var);

    // All the elements
    const_elem_iterator       elem_it (mesh.elements_begin());
    const const_elem_iterator elem_end(mesh.elements_end());

    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_vars(this->sys_number(),n_var);
  }
  

  //------------------------------------------------------------
  // Next allocate space for the DOF indices
  for (unsigned int var=0; var<this->n_variables(); var++)
    {
      const FEType& fe_type = this->variable_type(var);

      // For all the active elements
      const_active_elem_iterator       elem_it (mesh.elements_begin());
      const const_active_elem_iterator elem_end(mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	{
	  Elem*    elem       = *elem_it;
	  const ElemType type = elem->type();
	     
	  // Allocate the nodal DOFs
	  for (unsigned int n=0; n<elem->n_nodes(); n++)
	    {
	      const unsigned int dofs_at_node =
		FEInterface::n_dofs_at_node(dim, fe_type, type, n);
	      
	      elem->get_node(n)->set_n_comp(this->sys_number(),
					    var,
					    dofs_at_node);
	    }
	     
	  // Allocate the element DOFs
	  const unsigned int dofs_per_elem =
	    FEInterface::n_dofs_per_elem(dim, fe_type, type);
	  
	  elem->set_n_comp(this->sys_number(), var, dofs_per_elem);
	}
    }


  //------------------------------------------------------------
  // Finally, clear all the current DOF indices
  {
    // All the nodes
    const_node_iterator       node_it  (mesh.nodes_begin());
    const const_node_iterator node_end (mesh.nodes_end());

    for ( ; node_it != node_end; ++node_it)
      (*node_it)->invalidate_dofs(this->sys_number());
    
    // All the elements
    const_active_elem_iterator       elem_it (mesh.elements_begin());
    const const_active_elem_iterator elem_end(mesh.elements_end());

    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->invalidate_dofs(this->sys_number());
  }

  
  STOP_LOG("reinit()", "DofMap");
}



void DofMap::clear()
{
  // we don't want to clear
  // the coupling matrix!
  // It should not change...
  //_dof_coupling.clear();
  
  _first_df.clear();
  _last_df.clear();
  _send_list.clear();
  _n_nz.clear();
  _n_oz.clear();


#ifdef ENABLE_AMR

  _dof_constraints.clear();

#endif

  _matrix = NULL;

  _other_matrices.clear();

  _n_dfs = 0;
}




void DofMap::distribute_dofs(const MeshBase& mesh)
{
  assert (mesh.is_prepared());
  
  const unsigned int proc_id = libMesh::processor_id();
  const unsigned int n_proc  = libMesh::n_processors();
  const unsigned int sys_num = this->sys_number();
  const unsigned int n_vars  = this->n_variables();
  
  assert (n_vars > 0);
  assert (proc_id < n_proc);
  
  // re-init in case the mesh has changed
  this->reinit(mesh);

  // Log how long it takes to distribute the degrees of freedom
  START_LOG("distribute_dofs()", "DofMap");

  // The DOF counter, will be incremented as we encounter
  // new degrees of freedom
  unsigned int next_free_dof = 0;

  // Resize the _first_df and _last_df arrays, fill with 0.
  _first_df.resize(n_proc); std::fill(_first_df.begin(), _first_df.end(), 0);
  _last_df.resize (n_proc); std::fill(_last_df.begin(),  _last_df.end(),  0);

  _send_list.clear();
  
  //-------------------------------------------------------------------------
  // DOF numbering
  for (unsigned int processor=0; processor<n_proc; processor++)
    {
      _first_df[processor] = next_free_dof;
      
      for (unsigned var=0; var<n_vars; var++)
	{
	  const_active_pid_elem_iterator       elem_it (mesh.elements_begin(),
							processor);
	  const const_active_pid_elem_iterator elem_end(mesh.elements_end(),
							processor);

	  for ( ; elem_it != elem_end; ++elem_it)
	    {
	      // Only number dofs connected to active
	      // elements on this processor.
	      
	      Elem* elem                 = *elem_it;
	      const unsigned int n_nodes = elem->n_nodes();
	      
	      // First number the nodal DOFS
	      for (unsigned int n=0; n<n_nodes; n++)
		{
		  Node* node = elem->get_node(n);
		  
		  for (unsigned int index=0; index<node->n_comp(sys_num,var);
		       index++)
		    
		    // only assign a dof number if there isn't one
		    // already there
		    if (node->dof_number(sys_num,var,index) ==
			DofObject::invalid_id)
		      {
			node->set_dof_number(sys_num,
					     var,
					     index,
					     next_free_dof++);

			// If this DOF is on the local processor add it
			// to the _send_list
			if (processor == proc_id)
			  _send_list.push_back(node->dof_number(sys_num,
								var,
								index));
		      }
		}
		  
	      // Now number the element DOFS
	      for (unsigned int index=0; index<elem->n_comp(sys_num,var);
		   index++)
		{		  
		  // No way we could have already numbered this DOF!
		  assert (elem->dof_number(sys_num,var,index) ==
			  DofObject::invalid_id);
		  
		  elem->set_dof_number(sys_num,
				       var,
				       index,
				       next_free_dof++);
		  
		  // If this DOF is on the local processor add it
		  // to the _send_list
		  if (processor == proc_id)
		    _send_list.push_back(elem->dof_number(sys_num,
							  var,
							  index));
		}
	    }
	}
      
      _last_df[processor] = (next_free_dof-1);
    }

  // Set the total number of degrees of freedom
  _n_dfs = next_free_dof;
  //-------------------------------------------------------------------------

  

  //-------------------------------------------------------------------------
  // The _send_list now only contains entries from the local processor.
  // We need to add the DOFs from elements that live on neighboring processors
  // that are neighbors of the elements on the local processor
  {
    const_active_local_elem_iterator       elem_it (mesh.elements_begin());
    const const_active_local_elem_iterator elem_end(mesh.elements_end());

    std::vector<unsigned int> di;

    // Loop over the active local elements
    for ( ; elem_it != elem_end; ++elem_it)
      {
	Elem* elem = *elem_it;

	// Loop over the neighbors of those elements
	for (unsigned int s=0; s<elem->n_neighbors(); s++)
	  if (elem->neighbor(s) != NULL)
	    
	    // If the neighbor lives on a different processoe
	    if (elem->neighbor(s)->processor_id() != libMesh::processor_id())
	      {
		// Get the DOF indices for this neighboring element
		this->dof_indices (elem->neighbor(s), di);

		// Insert the DOF indices into the send list
		_send_list.insert (_send_list.end(),
				   di.begin(), di.end());
	      }
      }
  }
  //-------------------------------------------------------------------------


  
  // Note that in the code above nodes on processor boundaries
  // that are shared by multiple elements are added for each element.
  // Here we need to clean up that data structure
  this->sort_send_list ();

  // All done. Stop logging.
  STOP_LOG("distribute_dofs()", "DofMap");    
}



void DofMap::sort_send_list ()
{
  // First sort the send list.  After this
  // duplicated elements will be adjacent in the
  // vector
  std::sort(_send_list.begin(), _send_list.end());

  // Now use std::unique to remove duplicate entries  
  std::vector<unsigned int>::iterator new_end =
    std::unique (_send_list.begin(), _send_list.end());

  // Remove the end of the send_list.  Use the "swap trick"
  // from Effective STL
  std::vector<unsigned int> (_send_list.begin(), new_end).swap (_send_list);
}



void DofMap::compute_sparsity(const MeshBase& mesh)
{
  assert (mesh.is_prepared());
  assert (this->n_variables());

  START_LOG("compute_sparsity()", "DofMap");


  /**
   * Compute the sparsity structure of the global matrix.  This can be
   * fed into a PetscMatrix to allocate exacly the number of nonzeros
   * necessary to store the matrix.  This algorithm should be linear
   * in the (# of elements)*(# nodes per element)
   */
  const unsigned int proc_id           = mesh.processor_id();
  const unsigned int n_dofs_on_proc    = this->n_dofs_on_processor(proc_id);
  const unsigned int first_dof_on_proc = this->first_dof(proc_id);
  const unsigned int last_dof_on_proc  = this->last_dof(proc_id);
  
  std::vector<std::vector<unsigned int> > sparsity_pattern (n_dofs_on_proc);


  /**
   * If the user did not explicitly specify the DOF coupling
   * then all the DOFS are coupled to each other.  Furthermore,
   * we can take a shortcut and do this more quickly here.  So
   * we use an if-test.
   */  
  if (_dof_coupling.empty())
    {
      std::vector<unsigned int> element_dofs;

      const_active_elem_iterator       elem_it (mesh.elements_begin());
      const const_active_elem_iterator elem_end(mesh.elements_end());

      for ( ; elem_it != elem_end; ++elem_it)
	{
	  Elem* elem = *elem_it;
	  
	  this->dof_indices (elem, element_dofs);
	  this->find_connected_dofs (element_dofs);
	    
	  const unsigned int n_dofs_on_element = element_dofs.size();
	    
	  for (unsigned int i=0; i<n_dofs_on_element; i++)
	    {
	      const unsigned int ig = element_dofs[i];
	      
	      // Only bother if this matrix row will be stored
	      // on this processor.
	      if ((ig >= first_dof_on_proc) &&
		  (ig <= last_dof_on_proc))
		{
		  // This is what I mean
		  // assert ((ig - first_dof_on_proc) >= 0);
		  // but do the test like this because ig and
		  // first_dof_on_proc are unsigned ints
		  assert (ig >= first_dof_on_proc);
		  assert ((ig - first_dof_on_proc) < sparsity_pattern.size());
		  
		  std::vector<unsigned int>& row =
		    sparsity_pattern[ig - first_dof_on_proc];
		  
		  for (unsigned int j=0; j<n_dofs_on_element; j++)
		    {
		      const unsigned int jg = element_dofs[j];

		      // See if jg is in the sorted range
		      std::pair<std::vector<unsigned int>::iterator,
			        std::vector<unsigned int>::iterator>
			pos = std::equal_range (row.begin(), row.end(), jg);

		      // Insert jg if it wasn't found
		      if (pos.first == pos.second)
			row.insert (pos.first, jg);
		    }
		}
	    }
	}      
    } 


  
  /**
   * This is what we do in the case that the user has specified
   * explicit DOF coupling.
   */  
  else
    {
      assert (_dof_coupling.size() ==
	      this->n_variables());
      
      const unsigned int n_var = this->n_variables();
      
      std::vector<unsigned int> element_dofs_i;
      std::vector<unsigned int> element_dofs_j;
      
      
      const_active_elem_iterator       elem_it (mesh.elements_begin());
      const const_active_elem_iterator elem_end(mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	for (unsigned int vi=0; vi<n_var; vi++)
	  {
	    Elem* elem = *elem_it;
	    
	    // Find element dofs for variable vi
	    this->dof_indices (elem, element_dofs_i, vi);
	    this->find_connected_dofs (element_dofs_i);
	    
	    const unsigned int n_dofs_on_element_i = element_dofs_i.size();
	    
	    for (unsigned int vj=0; vj<n_var; vj++)
	      if (_dof_coupling(vi,vj)) // If vi couples to vj
		{
		  // Find element dofs for variable vj, note that
		  // if vi==vj we already have the dofs.
		  if (vi != vj)
		    {
		      this->dof_indices (elem, element_dofs_j, vj);
		      this->find_connected_dofs (element_dofs_j);
		    }
		  else
		    {
		      element_dofs_j = element_dofs_i;
		    }
		  
		  const unsigned int n_dofs_on_element_j =
		    element_dofs_j.size();
		  
		  for (unsigned int i=0; i<n_dofs_on_element_i; i++)
		    {
		      const unsigned int ig = element_dofs_i[i];
		      
		      // Only bother if this matrix row will be stored
		      // on this processor.
		      if ((ig >= first_dof_on_proc) &&
			  (ig <= last_dof_on_proc))
			{
			  // This is what I mean
			  // assert ((ig - first_dof_on_proc) >= 0);
			  // but do the test like this because ig and
			  // first_dof_on_proc are unsigned ints
			  assert (ig >= first_dof_on_proc);
			  assert (ig < (sparsity_pattern.size() +
					first_dof_on_proc));
			  
			  std::vector<unsigned int>& row =
			    sparsity_pattern[ig - first_dof_on_proc];
			  
			  for (unsigned int j=0; j<n_dofs_on_element_j; j++)
			    {
			      const unsigned int jg = element_dofs_j[j];

			      // See if jg is in the sorted range
			      std::pair<std::vector<unsigned int>::iterator,
			                std::vector<unsigned int>::iterator>
				pos = std::equal_range (row.begin(),
							row.end(),
							jg);

			      // Insert jg if it wasn't found
			      if (pos.first == pos.second)
				row.insert (pos.first, jg);
			    }
			}
		    }
		}
	  }
    }


  
  // Now the full sparsity structure is built for all of the
  // DOFs connected to our rows of the matrix.
  _n_nz.resize (n_dofs_on_proc);
  _n_oz.resize (n_dofs_on_proc);

  // First zero the counters.
  std::fill(_n_nz.begin(), _n_nz.end(), 0);
  std::fill(_n_oz.begin(), _n_oz.end(), 0);
  
  for (unsigned int i=0; i<n_dofs_on_proc; i++)
    {
      // Get the row of the sparsity pattern
      const std::vector<unsigned int>& row = sparsity_pattern[i];

      for (unsigned int j=0; j<row.size(); j++)
	if ((row[j] >= first_dof_on_proc) &&
	    (row[j] <= last_dof_on_proc))
	  _n_nz[i]++;
	else
	  _n_oz[i]++;
    }

  
  STOP_LOG("compute_sparsity()", "DofMap");

  
  // We are done with the sparsity_pattern.  However, quite a
  // lot has gone into computing it.  It is possible that some
  // \p SparseMatrix implementations want to see it.  Let them
  // see it before we throw it away.
  //
  // When additional matrices are handled, let the major matrix
  // be the first to get to know the sparsity pattern, so that
  // it has the best chance to be allocated successfully...
  if (_matrix != NULL)
    _matrix->update_sparsity_pattern (sparsity_pattern);
      
  std::map<std::string, SparseMatrix<Number>* >::iterator
    pos = _other_matrices.begin(),
    end = _other_matrices.end();
      
  for (; pos != end; ++pos)
    pos->second->update_sparsity_pattern (sparsity_pattern);     
}
 


void DofMap::dof_indices (const Elem* elem,
			  std::vector<unsigned int>& di,
			  const unsigned int vn) const
{
  START_LOG("dof_indices()", "DofMap");
  
  assert (elem != NULL);
  
  const unsigned int n_nodes = elem->n_nodes();
  const ElemType type        = elem->type();
  const unsigned int sys_num = this->sys_number();
  const unsigned int n_vars  = this->n_variables();
  
  // Clear the DOF indices vector
  di.clear();

  // Check that sizes match in DEBUG mode
#ifdef DEBUG
  unsigned int tot_size = 0;
#endif
  
  // Get the dof numbers
  for (unsigned int v=0; v<n_vars; v++)
    if ((v == vn) || (vn == static_cast<unsigned int>(-1)))
      { // Do this for all the variables if one was not specified
	// or just for the specified variable

#ifdef DEBUG
	const FEType& fe_type = this->variable_type(v);
	
	tot_size += FEInterface::n_dofs(elem->dim(),
					fe_type,
					type);
#endif
	
	// Get the node-based DOF numbers
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    const Node* node      = elem->get_node(n);
	    const unsigned int nc = node->n_comp(sys_num,v);
	    
	    for (unsigned int i=0; i<nc; i++)
	      {
		assert (node->dof_number(sys_num,v,i) !=
			DofObject::invalid_id);
		
		di.push_back(node->dof_number(sys_num,v,i));
	      }
	  }
	
	// Get the element-based DOF numbers
	{
	  const unsigned int nc = elem->n_comp(sys_num,v);
	  
	  for (unsigned int i=0; i<nc; i++)
	    {
	      assert (elem->dof_number(sys_num,v,i) !=
		      DofObject::invalid_id);
	      
	      di.push_back(elem->dof_number(sys_num,v,i));
	    }
	}
      }

#ifdef DEBUG
  assert (tot_size == di.size());
#endif
  
  STOP_LOG("dof_indices()", "DofMap");  
}
 


#ifdef ENABLE_AMR

void DofMap::old_dof_indices (const Elem* elem,
			      std::vector<unsigned int>& di,
			      const unsigned int vn) const
{
  START_LOG("old_dof_indices()", "DofMap");
  
  assert (elem != NULL);
  assert (elem->old_dof_object != NULL);
	    

  const unsigned int n_nodes = elem->n_nodes();
  const ElemType type        = elem->type();
  const unsigned int sys_num = this->sys_number();
  const unsigned int n_vars  = this->n_variables();
  
  // Clear the DOF indices vector.
  di.clear();

  // Check that sizes match in DEBUG mode
#ifdef DEBUG
  unsigned int tot_size = 0;
#endif
  
  // Get the dof numbers
  for (unsigned int v=0; v<n_vars; v++)
    if ((v == vn) || (vn == static_cast<unsigned int>(-1)))
      {
	// Do this for all the variables if one was not specified
	// or just for the specified variable
	
#ifdef DEBUG
	const FEType& fe_type = this->variable_type(v);
	tot_size += FEInterface::n_dofs(elem->dim(),
					fe_type,
					type);
#endif
	
	// Get the node-based DOF numbers
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    const Node* node      = elem->get_node(n);
	    assert (node->old_dof_object != NULL);
	    const unsigned int nc = node->old_dof_object->n_comp(sys_num,v);
	    
	    for (unsigned int i=0; i<nc; i++)
	      {
		assert (node->old_dof_object->dof_number(sys_num,v,i) !=
			DofObject::invalid_id);
		
		di.push_back(node->old_dof_object->dof_number(sys_num,v,i));
	      }
	  }
	
	// Get the element-based DOF numbers
	{
	  assert (elem->old_dof_object != NULL);
	  const unsigned int nc = elem->old_dof_object->n_comp(sys_num,v);
	  
	  for (unsigned int i=0; i<nc; i++)
	    {
	      assert (elem->old_dof_object->dof_number(sys_num,v,i) !=
		      DofObject::invalid_id);
	      
	      di.push_back(elem->old_dof_object->dof_number(sys_num,v,i));
	    }
	}
      }

#ifdef DEBUG
  assert (tot_size == di.size());
#endif
  
  STOP_LOG("old_dof_indices()", "DofMap");  
}



void DofMap::augment_send_list_for_projection(const MeshBase& mesh)
{
  // Loop over the active local elements in the mesh.
  // If the element was just refined and its parent lives
  // on a different processor then we need to augment the
  // _send_list with the parent's DOF indices so that we
  // can access the parent's data for computing solution
  // projections, etc...

  // The DOF indices for the parent
  std::vector<unsigned int> di;

  // Flag telling us if we need to re-sort the send_list.
  // Note we won't need to re-sort unless a child with a
  // parent belonging to a different processor is found
  bool needs_sorting = false;
  
  
  const_active_local_elem_iterator       elem_it (mesh.elements_begin());
  const const_active_local_elem_iterator elem_end(mesh.elements_end());

  for ( ; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem   = *elem_it;

      // We only need to consider the children that
      // were just refined
      if (elem->refinement_flag() != Elem::JUST_REFINED) continue;
      
      const Elem* parent = elem->parent();

      // If the parent lives on another processor
      // than the child
      if (parent != NULL)
	if (parent->processor_id() != elem->processor_id())
	  {
	    // Get the DOF indices for the parent
	    this->dof_indices (parent, di);

	    // Insert the DOF indices into the send list
	    _send_list.insert (_send_list.end(),
			       di.begin(), di.end());

	    // We will need to re-sort the send list
	    needs_sorting = true;	    
	  }
    }

  // The send-list might need to be sorted again.
  if (needs_sorting) this->sort_send_list ();
}

#endif // #ifdef ENABLE_AMR










void DofMap::find_connected_dofs (std::vector<unsigned int>& elem_dofs) const
{
  
#ifdef ENABLE_AMR

  typedef std::set<unsigned int> RCSet;

  // First insert the DOFS we already depend on into the set.  
  RCSet dof_set (elem_dofs.begin(), elem_dofs.end());

  bool done = true;
    
  // Next insert any dofs those might be constrained in terms
  // of.  Note that in this case we may not be done:  Those may
  // in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  for (unsigned int i=0; i<elem_dofs.size(); i++)
    if (this->is_constrained_dof(elem_dofs[i]))
      {
	// If the DOF is constrained
	DofConstraints::const_iterator
	  pos = _dof_constraints.find(elem_dofs[i]);
	
	assert (pos != _dof_constraints.end());
	
	const DofConstraintRow& constraint_row = pos->second;
	
	assert (!constraint_row.empty());
	
	DofConstraintRow::const_iterator it     = constraint_row.begin();
	DofConstraintRow::const_iterator it_end = constraint_row.end();
	
	
	// Add the DOFs this dof is constrained in terms of.
	// note that these dofs might also be constrained, so
	// we will need to call this function recursively.
	for ( ; it != it_end; ++it)
	  if (!dof_set.count (it->first))
	    {
	      dof_set.insert (it->first);
	      done = false;
	    }
	}
  
  
  // If not done then we need to do more work
  // (obviously :-) )!
  if (!done)
    {
      // Fill the vector with the contents of the set
      elem_dofs.clear();      
      elem_dofs.insert (elem_dofs.end(),
			dof_set.begin(), dof_set.end());
      

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      this->find_connected_dofs (elem_dofs);
      
    } // end if (!done)


#endif // #ifdef ENABLE_AMR


}
