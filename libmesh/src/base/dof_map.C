// $Id: dof_map.C,v 1.28 2003-03-04 15:31:10 benkirk Exp $

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
#include "dense_matrix.h"
#include "dense_vector.h"
#include "libmesh.h"





// ------------------------------------------------------------
// DofMap member functions
void DofMap::attach_matrix (SparseMatrix<Number>& matrix)
{
  _matrix = &matrix;
  
  _matrix->attach_dof_map (*this);
}



void DofMap::reinit(MeshBase& mesh)
{
  assert (mesh.is_prepared());
  
  START_LOG("reinit()", "DofMap");
  
  clear();

  _n_nodes = mesh.n_nodes();
  _n_elem  = mesh.n_elem();

  const unsigned int n_var = this->n_variables();
  const unsigned int dim   = mesh.mesh_dimension();

  // First set the number of variables for each \p DofObject
  // equal to n_variables() for this system.
  {
    // All the nodes
    node_iterator       node_it  (mesh.nodes_begin());
    const node_iterator node_end (mesh.nodes_end());

    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_vars(this->sys_number(),n_var);
    
    // All the elements
    elem_iterator       elem_it (mesh.elements_begin());
    const elem_iterator elem_end(mesh.elements_end());

    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_vars(this->sys_number(),n_var);
  }
  

  
  for (unsigned int var=0; var<this->n_variables(); var++)
    {
      const FEType& fe_type = variable_type(var);

      // For all the active elements
      active_elem_iterator       elem_it (mesh.elements_begin());
      const active_elem_iterator elem_end(mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	{
	  Elem*    elem       = *elem_it;
	  const ElemType type = elem->type();
	     
	  // Allocate the nodal DOFs
	  for (unsigned int n=0; n<elem->n_nodes(); n++)
	    {
	      const unsigned int dofs_at_node =
		FEInterface::n_dofs_at_node(dim, fe_type, type, n);
	      
	      if (dofs_at_node > elem->get_node(n)->n_comp(this->sys_number(),var))
		elem->get_node(n)->set_n_comp(this->sys_number(), var, dofs_at_node);
	    }
	     
	  // Allocate the element DOFs
	  const unsigned int dofs_per_elem =
	    FEInterface::n_dofs_per_elem(dim, fe_type, type);
	  
	  if (dofs_per_elem > elem->n_comp(this->sys_number(), var))
	    elem->set_n_comp(this->sys_number(), var, dofs_per_elem);
	}
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
  
  _n_dfs =
    _n_nodes =
    _n_elem =
    0;
}




void DofMap::distribute_dofs(MeshBase& mesh)
{
  assert (mesh.is_prepared());
  
  const unsigned int proc_id = mesh.processor_id();
  const unsigned int n_proc  = mesh.n_processors();
  
  assert (this->n_variables() != 0);
  assert (proc_id < mesh.n_processors());
  
  // re-init in case the mesh has changed
  this->reinit(mesh);

  assert (_n_nodes);
  assert (_n_elem);

  // Log how long it takes to distribute the degrees of freedom
  START_LOG("distribute_dofs()", "DofMap");

  unsigned int next_free_dof=0;

  // Resize the _first_df and _last_df arrays, fill with 0.
  _first_df.resize(n_proc); std::fill(_first_df.begin(), _first_df.end(), 0);
  _last_df.resize (n_proc); std::fill(_last_df.begin(),  _last_df.end(),  0);

  _send_list.clear();
  
  //------------------------------------------------------------
  // DOF numbering
  for (unsigned int processor=0; processor<n_proc; processor++)
    {
      _first_df[processor] = next_free_dof;
      
      for (unsigned var=0; var<this->n_variables(); var++)
	{
	  active_pid_elem_iterator       elem_it (mesh.elements_begin(), processor);
	  const active_pid_elem_iterator elem_end(mesh.elements_end(),   processor);

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
		  
		  for (unsigned int index=0; index<node->n_comp(this->sys_number(),var); index++)
		    {
		      // only assign a dof number if there isn't one
		      // already there
		      if (node->dof_number(this->sys_number(),var,index) == DofObject::invalid_id)
			{
			  node->set_dof_number(this->sys_number(),var,index,next_free_dof++);
			  
			  if (processor == proc_id)
			    _send_list.push_back(node->dof_number(this->sys_number(),var,index));
			}
		      // if there is an entry there and it isn't on this
		      // processor it also needs to be added to the send list
		      else if (node->dof_number(this->sys_number(),var,index) <
			       _first_df[proc_id])
			{
			  if (processor == proc_id)
			    _send_list.push_back(node->dof_number(this->sys_number(),var,index));
			}
		    }
		}
		  
	      // Now number the element DOFS
	      for (unsigned int index=0; index<elem->n_comp(this->sys_number(),var); index++)
		{		  
		  // No way we could have already numbered this DOF!
		  assert (elem->dof_number(this->sys_number(),var,index) ==
			  DofObject::invalid_id);
		  
		  elem->set_dof_number(this->sys_number(),var,index,next_free_dof++);
		  
		  if (processor == proc_id)
		    _send_list.push_back(elem->dof_number(this->sys_number(),var,index));
		}
	    }
	}
      
      _last_df[processor] = (next_free_dof-1);
    }
  
  _n_dfs = next_free_dof;

  
  // Note that in the code above nodes on processor boundaries
  // that are shared by multiple elements are added for each element.
  // Here we need to clean up that data structure
  {    
    std::sort(_send_list.begin(), _send_list.end());
    std::vector<unsigned int> new_send_list;

    // The new send list will be <= the size of the
    // current _send_list.  Let's reserve memory
    // assuming it will be the same size.  This will
    // allow us to efficiently use the push_back()
    // member (and probably will actually _save_
    // space since vector reallocations generally
    // occur in powers of 2).
    new_send_list.reserve (_send_list.size());
    
    new_send_list.push_back(_send_list[0]);
    
    for (unsigned int i=1; i<_send_list.size(); i++)
      if (_send_list[i] != new_send_list.back())
	new_send_list.push_back(_send_list[i]);

    _send_list = new_send_list;
  }

  // All done. Stop logging.
  STOP_LOG("distribute_dofs()", "DofMap");    
}



void DofMap::compute_sparsity(MeshBase& mesh)
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
  
  std::vector<std::set<unsigned int> > sparsity_pattern (n_dofs_on_proc);


  /**
   * If the user did not explicitly specify the DOF coupling
   * then all the DOFS are coupled to each other.  Furthermore,
   * we can take a shortcut and do this more quickly here.  So
   * we use an if-test.
   */  
  if (_dof_coupling.empty())
    {
      std::vector<unsigned int> element_dofs;
      
      for (unsigned int e=0; e<_n_elem; e++)
	if (mesh.elem(e)->active())
	  {
	    this->dof_indices (mesh.elem(e), element_dofs);
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
		    
		    std::set<unsigned int>& row =
		      sparsity_pattern[ig - first_dof_on_proc];
		    
		    for (unsigned int j=0; j<n_dofs_on_element; j++)
		      {
			const unsigned int jg = element_dofs[j];
			
			row.insert(jg);
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
      
      
      for (unsigned int e=0; e<_n_elem; e++)
	if (mesh.elem(e)->active())
	  for (unsigned int vi=0; vi<n_var; vi++)
	    {
	      // Find element dofs for variable vi
	      this->dof_indices (mesh.elem(e), element_dofs_i, vi);
	      this->find_connected_dofs (element_dofs_i);
	      
	      const unsigned int n_dofs_on_element_i = element_dofs_i.size();

	      for (unsigned int vj=0; vj<n_var; vj++)
		if (_dof_coupling(vi,vj)) // If vi couples to vj
		  {
		    // Find element dofs for variable vj
		    this->dof_indices (mesh.elem(e), element_dofs_j, vj);
		    this->find_connected_dofs (element_dofs_j);	    
		    
		    const unsigned int n_dofs_on_element_j = element_dofs_j.size();
	    
		    for (unsigned int i=0; i<n_dofs_on_element_i; i++)
		      {
			const unsigned int ig = element_dofs_i[i];

			// Only bother if this matrix row will be stored
			// on this processor.
			if ((ig >= first_dof_on_proc) &&
			    (ig <= last_dof_on_proc))
			  {
			    // This is what I mean
			    //assert ((ig - first_dof_on_proc) >= 0);
			    // but do the test like this because ig and
			    // first_dof_on_proc are unsigned ints
			    assert (ig >= first_dof_on_proc);
			    assert (ig < (sparsity_pattern.size() + first_dof_on_proc));
		    
			    std::set<unsigned int>& row =
			      sparsity_pattern[ig - first_dof_on_proc];
		    
			    for (unsigned int j=0; j<n_dofs_on_element_j; j++)
			      {
				const unsigned int jg = element_dofs_j[j];
			
				row.insert(jg);
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
      const std::set<unsigned int>& row = sparsity_pattern[i];
      
      // Now loop over the row and increment the appropriate
      // counters.
      std::set<unsigned int>::const_iterator it_end = row.end();
      
      for (std::set<unsigned int>::const_iterator it=row.begin();
	   it != it_end; ++it)
	if ((*it >= first_dof_on_proc) &&
	    (*it <= last_dof_on_proc))
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
  // NOTE:  The \p SparseMatrix::update_sparsity_pattern() is NOT
  // const, so the \p SparseMatrix may freely trash the
  // sparsity_pattern.  DO NOT expect to use it any more after
  // this call.
  {
    if (_matrix != NULL)
      _matrix->update_sparsity_pattern (sparsity_pattern);
    
    // Explicity clear the sparsity pattern now.
    // it is about to go out of scope, but maybe this
    // will expedite freeing up its memory?
    sparsity_pattern.clear();
  }
}



void DofMap::dof_indices (const Elem* elem,
			  std::vector<unsigned int>& di,
			  const unsigned int vn) const
{
  assert (elem != NULL);


  const unsigned int n_nodes = elem->n_nodes();
  const ElemType type        = elem->type();
  
  // Clear the DOF indices vector.
  di.clear();

  unsigned int tot_size = 0;
  
  // Get the dof numbers
  for (unsigned int v=0; v<this->n_variables(); v++)
    if ((v == vn) || (vn == static_cast<unsigned int>(-1)))
      { // Do this for all the variables if one was not specified
	// or just for the specified variable
	const FEType& fe_type = this->variable_type(v);

	const unsigned int size = FEInterface::n_dofs(elem->dim(), fe_type, type);
	tot_size += size;
	// Reserve space in the di vector
	// so we can use push_back effectively
	//di.reserve (size);
	
	// Get the node-based DOF numbers
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    const Node* node = elem->get_node(n);
	    
	    for (unsigned int i=0; i<node->n_comp(this->sys_number(),v); i++)
	      {
		assert (node->dof_number(this->sys_number(),v,i) !=
			DofObject::invalid_id);
		
		di.push_back(node->dof_number(this->sys_number(),v,i));
	      }
	  }
	
	// Get the element-based DOF numbers	  
	for (unsigned int i=0; i<elem->n_comp(this->sys_number(),v); i++)
	  {
	    assert (elem->dof_number(this->sys_number(),v,i) !=
		    DofObject::invalid_id);
	    
	    di.push_back(elem->dof_number(this->sys_number(),v,i));
	  }
      }

  assert (tot_size == di.size());
}








#ifdef ENABLE_AMR


void DofMap::create_dof_constraints(MeshBase& mesh)
{
  assert (mesh.is_prepared());
  
  START_LOG("create_dof_constraints()", "DofMap");

  const unsigned int dim = mesh.mesh_dimension();

  // Constraints are not necessary in 1D
  if (dim == 1)
    return;
  
  // Here we build the hanging node constraints.  This is done
  // by enforcing the condition u_a = u_b along hanging sides.
  // u_a = u_b is collocated at the nodes of side a, which gives
  // one row of the constraint matrix.
  
  // clear any existing constraints.
  _dof_constraints.clear();
  
  // Look at the element faces.  Check to see if we need to 
  // build constraints.
  for (unsigned int e=0; e<_n_elem; e++)
    if (mesh.elem(e)->active())
      for (unsigned int s=0; s<mesh.elem(e)->n_sides(); s++)
	if (mesh.elem(e)->neighbor(s) != NULL)
	  if (mesh.elem(e)->neighbor(s)->level() < mesh.elem(e)->level()) // constrain dofs shared between
	    {                                                               // this element and ones coarser
	                                                                    // than this element.
	      // Get pointers to the elements of interest and its parent.
	      const Elem* elem   = mesh.elem(e);
	      const Elem* parent = elem->parent();

	      // This can't happen...  Only level-0 elements have NULL
	      // parents, and no level-0 elements can be at a higher
	      // level than their neighbors!
	      assert (parent != NULL);
	      
	      const AutoPtr<Elem> my_side     (elem->build_side(s));
	      const AutoPtr<Elem> parent_side (parent->build_side(s));

	      // Look at all the variables in the system
	      for (unsigned int variable=0; variable<this->n_variables();
		   ++variable)
		{
		  const FEType& fe_type = this->variable_type(variable);
	      
		  for (unsigned int my_dof=0;
		       my_dof<FEInterface::n_dofs(dim-1, fe_type, my_side->type());
		       my_dof++)
		    {
		      assert (my_dof < my_side->n_nodes());

		      // My global node and dof indices.
		      const Node* my_node          = my_side->get_node(my_dof);
		      const unsigned int my_dof_g  = my_node->dof_number(this->sys_number(), variable, 0);

		      // The support point of the DOF
		      const Point& support_point = my_side->point(my_dof);

		      // Figure out where my node lies on their reference element.
		      const Point mapped_point = FEInterface::inverse_map(dim-1, fe_type,
									  parent_side.get(),
									  support_point);

		      // Compute the parent's side shape function values.
		      for (unsigned int their_dof=0;
			   their_dof<FEInterface::n_dofs(dim-1, fe_type, parent_side->type());
			   their_dof++)
			{
			  assert (their_dof < parent_side->n_nodes());
			  
			  // Their global node and dof indices.
			  const Node* their_node          = parent_side->get_node(their_dof);
			  const unsigned int their_dof_g  = their_node->dof_number(this->sys_number(), variable, 0);

			  const Real their_dof_value = FEInterface::shape(dim-1,
									  fe_type,
									  parent_side->type(),
									  their_dof,
									  mapped_point);

			  // Only add non-zero and non-identity values
			  // for Lagrange basis functions.
			  if ((fabs(their_dof_value) > 1.e-5) &&
			      (fabs(their_dof_value) < .999)) 
			    {
			      // A reference to the constraint row.
			      DofConstraintRow& constraint_row = _dof_constraints[my_dof_g];
			      			      
			      const std::pair<unsigned int, float> p (their_dof_g,
								      static_cast<float>(their_dof_value));

			      constraint_row.insert(p);
			    }
			}		      
		    }
		}
	    }
  
  STOP_LOG("create_dof_constraints()", "DofMap");
}



void DofMap::add_constraint_row (const unsigned int dof_number,
				 const DofConstraintRow& constraint_row)
{

  if (!_dof_constraints.count(dof_number))
    std::cerr << "WARNING: DOF " << dof_number << " was already constrained!"
	      << std::endl;

  std::pair<unsigned int, DofConstraintRow> kv(dof_number, constraint_row);

  _dof_constraints.insert(kv);
}



void DofMap::print_dof_constraints() const
{
  std::cout << "DOF CONSTRAINTS OUTPUT:"
	    << std::endl;
  
  for (DofConstraints::const_iterator it=_dof_constraints.begin();
       it != _dof_constraints.end(); ++it)
    {
      const unsigned int i = it->first;
      const DofConstraintRow& row = it->second;

      std::cout << "Constraints for DOF " << i
		<< ": \t";

      for (DofConstraintRow::const_iterator pos=row.begin();
	   pos != row.end(); ++pos)
	std::cout << " (" << pos->first << ","
		  << pos->second << ")\t";

      std::cout << std::endl;
    }
}



void DofMap::constrain_element_matrix (DenseMatrix<Number>& matrix,
				       std::vector<unsigned int>& elem_dofs) const
{
  assert (elem_dofs.size() == matrix.m());
  assert (elem_dofs.size() == matrix.n());
  
  // The constrained matrix is built up as C^T K C.    
  DenseMatrix<Number> C;

  
  this->build_constraint_matrix (C, elem_dofs);

  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      //here();
    
      matrix.left_multiply  (C, true);
      matrix.right_multiply (C, false);
      
      
      assert (matrix.m() == matrix.n());
      assert (matrix.m() == elem_dofs.size());
      assert (matrix.n() == elem_dofs.size());
      
      
      for (unsigned int i=0; i<elem_dofs.size(); i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;
	    
	    // If the DOF is constrained
	    matrix(i,i) = 1.;
	    
	    DofConstraints::const_iterator
	      pos = _dof_constraints.find(elem_dofs[i]);
	    
	    assert (pos != _dof_constraints.end());
	    
	    const DofConstraintRow& constraint_row = pos->second;
	    
	    assert (!constraint_row.empty());
	    
	    for (DofConstraintRow::const_iterator
		   it=constraint_row.begin(); it != constraint_row.end();
		 ++it)
	      for (unsigned int j=0; j<elem_dofs.size(); j++)
		if (elem_dofs[j] == it->first)
		  matrix(i,j) = -it->second;	
	  }
    } // end if is constrained...
}



void DofMap::constrain_element_matrix_and_vector (DenseMatrix<Number>& matrix,
						  std::vector<Number>& rhs,
						  std::vector<unsigned int>& elem_dofs) const
{
  assert (elem_dofs.size() == matrix.m());
  assert (elem_dofs.size() == matrix.n());
  assert (elem_dofs.size() == rhs.size());
  
  // The constrained matrix is built up as C^T K C.
  // The constrained RHS is built up as C^T F
  DenseMatrix<Number> C;
  
  this->build_constraint_matrix (C, elem_dofs);

  
  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // Compute the matrix-matrix-matrix product C^T K C
      matrix.left_multiply  (C, true);
      matrix.right_multiply (C, false);
      
      
      assert (matrix.m() == matrix.n());
      assert (matrix.m() == elem_dofs.size());
      assert (matrix.n() == elem_dofs.size());
      
      
      for (unsigned int i=0; i<elem_dofs.size(); i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;
	    
	    // If the DOF is constrained
	    matrix(i,i) = 1.;
	    
	    DofConstraints::const_iterator
	      pos = _dof_constraints.find(elem_dofs[i]);
	    
	    assert (pos != _dof_constraints.end());
	    
	    const DofConstraintRow& constraint_row = pos->second;
	    
	    assert (!constraint_row.empty());
	    
	    for (DofConstraintRow::const_iterator
		   it=constraint_row.begin(); it != constraint_row.end();
		 ++it)
	      for (unsigned int j=0; j<elem_dofs.size(); j++)
		if (elem_dofs[j] == it->first)
		  matrix(i,j) = -it->second;	
	  }

      
      // Compute the matrix-vector product C^T F
      std::vector<Number> old_rhs(rhs);
      
      rhs.resize(elem_dofs.size());
      
      
      for (unsigned int i=0; i<elem_dofs.size(); i++)
	{
	  // zero before summation
	  rhs[i] = 0.;
	  
	  for (unsigned int j=0; j<old_rhs.size(); j++)
	    rhs[i] += C.transpose(i,j)*old_rhs[j];
	}


      assert (elem_dofs.size() == rhs.size());

      for (unsigned int i=0; i<elem_dofs.size(); i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {	
	    // If the DOF is constrained
	    rhs[i] = 0.;
	  }
    } // end if is constrained...
}



void DofMap::constrain_element_matrix (DenseMatrix<Number>& matrix,
				       std::vector<unsigned int>& row_dofs,
				       std::vector<unsigned int>& col_dofs) const
{
  assert (row_dofs.size() == matrix.m());
  assert (col_dofs.size() == matrix.n());
  
  // The constrained matrix is built up as R^T K C.


  
  
  DenseMatrix<Number> R;
  DenseMatrix<Number> C;

  // Safeguard against the user passing us the same
  // object for row_dofs and col_dofs.  If that is done
  // the calls to build_matrix would fail
  std::vector<unsigned int> orig_row_dofs(row_dofs);
  std::vector<unsigned int> orig_col_dofs(col_dofs);
  
  this->build_constraint_matrix (R, orig_row_dofs);
  this->build_constraint_matrix (C, orig_col_dofs);

  row_dofs = orig_row_dofs;
  col_dofs = orig_col_dofs;
  
    
  // It is possible that the matrix is not constrained at all.
  if ((R.m() == matrix.m()) &&
      (R.n() == row_dofs.size()) &&
      (C.m() == matrix.n()) &&
      (C.n() == col_dofs.size())) // It the matrix is constrained
    {
      //here();
    
      matrix.left_multiply  (R, true);
      matrix.right_multiply (C, false);
      
      
      assert (matrix.m() == row_dofs.size());
      assert (matrix.n() == col_dofs.size());
      
      
      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (this->is_constrained_dof(row_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;
	    
	    // If the DOF is constrained
	    matrix(i,i) = 1.;
	    
	    DofConstraints::const_iterator
	      pos = _dof_constraints.find(row_dofs[i]);
	    
	    assert (pos != _dof_constraints.end());
	    
	    const DofConstraintRow& constraint_row = pos->second;
	    
	    assert (!constraint_row.empty());
	    
	    for (DofConstraintRow::const_iterator
		   it=constraint_row.begin(); it != constraint_row.end();
		 ++it)
	      for (unsigned int j=0; j<col_dofs.size(); j++)
		if (col_dofs[j] == it->first)
		  matrix(i,j) = -it->second;	
	  }
    } // end if is constrained...
}



void DofMap::constrain_element_vector (std::vector<Number>&       rhs,
				       std::vector<unsigned int>& row_dofs) const
{
  assert (rhs.size() == row_dofs.size());

  
  // The constrained RHS is built up as R^T F.  
  DenseMatrix<Number> R;

  this->build_constraint_matrix (R, row_dofs);

  
  // It is possible that the vector is not constrained at all.
  if ((R.m() == rhs.size()) &&
      (R.n() == row_dofs.size())) // if the RHS is constrained
    {
      // Compute the matrix-vector product
      std::vector<Number> old_rhs(rhs);
      
      rhs.resize(row_dofs.size());
      
      
      for (unsigned int i=0; i<row_dofs.size(); i++)
	{
	  // zero before summation
	  rhs[i] = 0.;
	  
	  for (unsigned int j=0; j<old_rhs.size(); j++)
	    rhs[i] += R.transpose(i,j)*old_rhs[j];
	}


      assert (row_dofs.size() == rhs.size());

      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (this->is_constrained_dof(row_dofs[i]))
	  {	
	    // If the DOF is constrained
	    rhs[i] = 0.;
	  }
    } // end if the RHS is constrained.
}



void DofMap::build_constraint_matrix (DenseMatrix<Number>& C,
				      std::vector<unsigned int>& elem_dofs) const
{
  typedef std::set<unsigned int> RCSet;
  RCSet dof_set;

  bool done = true;

  // First insert the DOFS we already depend on into the sets.
  {
    for (unsigned int i=0; i<elem_dofs.size(); i++)
      dof_set.insert(elem_dofs[i]);
  }
  

  // Next insert any dofs those might be constrained in terms
  // of.  Note that in this case we may not be done:  Those may
  // in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  {
    const unsigned int orig_dof_set_size = dof_set.size();
    
    for (unsigned int i=0; i<elem_dofs.size(); i++)
      if (this->is_constrained_dof(elem_dofs[i]))
	{
	  // If the DOF is constrained
	  DofConstraints::const_iterator
	    pos = _dof_constraints.find(elem_dofs[i]);
	  
	  assert (pos != _dof_constraints.end());
	  
	  const DofConstraintRow& constraint_row = pos->second;
	  
	  assert (!constraint_row.empty());
	  
	  for (DofConstraintRow::const_iterator
		 it=constraint_row.begin(); it != constraint_row.end();
	       ++it)
	    dof_set.insert (it->first);
	}

    // If we added any DOFS then we need to do this recursively.
    // It is possible that we just added a DOF that is also
    // constrained!
    if (dof_set.size() != orig_dof_set_size)
      done = false;
  }


  // If not done then we need to do more work
  // (obviously :-) )!
  if (!done)
    {
      std::vector<unsigned int> new_elem_dofs(dof_set.size());
      
      {
	RCSet::const_iterator end_it = dof_set.end(); 
	
	unsigned int i=0;
	
	for (RCSet::const_iterator it=dof_set.begin();
	     it != end_it; ++it)
	  new_elem_dofs[i++] = *it;
      }
      
      // Now we can build the constraint matrix.
      // Note that resize also zeros for a DenseMatrix<Number>.
      C.resize (elem_dofs.size(), new_elem_dofs.size());
      
      // Create the C constraint matrix.
      {
	for (unsigned int i=0; i<elem_dofs.size(); i++)
	  if (this->is_constrained_dof(elem_dofs[i]))
	    {
	      // If the DOF is constrained
	      DofConstraints::const_iterator
		pos = _dof_constraints.find(elem_dofs[i]);
	      
	      assert (pos != _dof_constraints.end());
	      
	      const DofConstraintRow& constraint_row = pos->second;
	      
	      assert (!constraint_row.empty());
	      
	      for (DofConstraintRow::const_iterator
		     it=constraint_row.begin(); it != constraint_row.end();
		   ++it)
		for (unsigned int j=0; j<new_elem_dofs.size(); j++)
		  if (new_elem_dofs[j] == it->first)
		    C(i,j) = it->second;
	    }
	  else
	    {
	      for (unsigned int j=0; j<new_elem_dofs.size(); j++)
		if (new_elem_dofs[j] == elem_dofs[i])
		  C(i,j) =  1.;
	    }	
	//C.print();
      }

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      elem_dofs = new_elem_dofs;
      
      DenseMatrix<Number> Cnew;
      
      this->build_constraint_matrix (Cnew, elem_dofs);

      if ((C.n() == Cnew.m()) &&
	  (Cnew.n() == elem_dofs.size())) // If the constraint matrix
	{                                 // is constrained...
	  //here();
	  C.right_multiply(Cnew, false);
	}
      
      assert (C.n() == elem_dofs.size());
    } // end if (!done)
}

#endif // #ifdef ENABLE_AMR











void DofMap::find_connected_dofs (std::vector<unsigned int>& elem_dofs) const
{
  
#ifdef ENABLE_AMR

  typedef std::set<unsigned int> RCSet;
  RCSet dof_set;

  bool done = true;

  // First insert the DOFS we already depend on into the sets.
  for (unsigned int i=0; i<elem_dofs.size(); i++)
    dof_set.insert(elem_dofs[i]);
  
  
  // Next insert any dofs those might be constrained in terms
  // of.  Note that in this case we may not be done:  Those may
  // in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  {
    const unsigned int orig_dof_set_size = dof_set.size();
    
    for (unsigned int i=0; i<elem_dofs.size(); i++)
      if (this->is_constrained_dof(elem_dofs[i]))
	{
	  // If the DOF is constrained
	  DofConstraints::const_iterator
	    pos = _dof_constraints.find(elem_dofs[i]);
	  
	  assert (pos != _dof_constraints.end());
	  
	  const DofConstraintRow& constraint_row = pos->second;
	  
	  assert (!constraint_row.empty());
	  
	  for (DofConstraintRow::const_iterator
		 it=constraint_row.begin(); it != constraint_row.end();
	       ++it)
	    dof_set.insert (it->first);
	}

    // If we added any DOFS then we need to do this recursively.
    // It is possible that we just added a DOF that is also
    // constrained!
    if (dof_set.size() != orig_dof_set_size)
      done = false;
  }


  // If not done then we need to do more work
  // (obviously :-) )!
  if (!done)
    {
      elem_dofs.resize(dof_set.size());
      
      {
	RCSet::const_iterator end_it = dof_set.end(); 
	
	unsigned int i=0;
	
	for (RCSet::const_iterator it=dof_set.begin();
	     it != end_it; ++it)
	  elem_dofs[i++] = *it;
      }
      

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      this->find_connected_dofs (elem_dofs);
      
    } // end if (!done)


#endif // #ifdef ENABLE_AMR


}
