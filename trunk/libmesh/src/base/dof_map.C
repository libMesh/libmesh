// $Id: dof_map.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <math.h>


// Local Includes -----------------------------------
#include "dof_map.h"
#include "elem.h"
#include "mesh_base.h"
#include "quadrature.h"
#include "fe_interface.h"





// ------------------------------------------------------------
// DofMap static member initializations
const unsigned int DofMap::invalid_number = static_cast<unsigned int>(-1);



// ------------------------------------------------------------
// DofMap member functions
void DofMap::reinit()
{
  perf_log.start_event("reinit()");
  
  clear();

  n_nodes = mesh.n_nodes();
  n_elem  = mesh.n_elem();

  std::vector<unsigned int> max_dofs_per_node(n_components(), 0);
  std::vector<unsigned int> max_dofs_per_elem(n_components(), 0);
  
  // First compute the maximum number of DOFs per node
  // and element for each component.  We may not need to allocate
  // memory at all...
   for (unsigned int c=0; c<n_components(); c++)
     {
       const FEType fe_type = component_type(c);
       
       for (unsigned int e=0; e<n_elem; e++)
	 if (mesh.elem(e)->active())
	   {
	     const Elem*    elem = mesh.elem(e);
	     const ElemType type = elem->type();
	     
	     // Count the nodal DOFs
	     for (unsigned int node=0; node<elem->n_nodes(); node++)
	       {
		 const unsigned int dofs_at_node =
		   FEInterface::n_dofs_at_node(dim, fe_type, type, node);
		 
		 max_dofs_per_node[c] = std::max(max_dofs_per_node[c],
						 dofs_at_node);
	       };
	     
	     // Count the element DOFs
	     const unsigned int dofs_per_elem =
	       FEInterface::n_dofs_per_elem(dim, fe_type, type);
	     
	     max_dofs_per_elem[c] = std::max(max_dofs_per_elem[c],
					     dofs_per_elem);
	   };
     };

  
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  // initialize the node_dofs and elem_dofs maps  
  node_dofs.resize(n_components());
  elem_dofs.resize(n_components());

  
   for (unsigned int c=0; c<n_components(); c++)
    {
      // No need to waste space...
      if (max_dofs_per_node[c] > 0)
	node_dofs[c].resize(n_nodes);

      if (max_dofs_per_elem[c] > 0)
	elem_dofs[c].resize(n_elem);

      const FEType fe_type = component_type(c);
      
      for (unsigned int e=0; e<n_elem; e++)
	if (mesh.elem(e)->active())
	  {
	    const Elem* elem    = mesh.elem(e);
	    const ElemType type = elem->type();
	    
	    // Initialize the nodal dofs
	    if (max_dofs_per_node[c] > 0)
	      for (unsigned int node=0; node<elem->n_nodes(); node++)
		{
		  const unsigned int global_node  = elem->node(node);
		  const unsigned int current_size =
		    node_dofs[c][global_node].size();
		  
		  const unsigned int dofs_per_node =
		    FEInterface::n_dofs_at_node(dim, fe_type, type, node);
		  
		  if (dofs_per_node > current_size)
		    {
		      node_dofs[c][global_node].resize(dofs_per_node);
		      
		      for (unsigned int i=0; i<dofs_per_node; i++)
			node_dofs[c][global_node][i] = invalid_number;
		    };
		};
	    
	    // Initialize the element dofs
	    if (max_dofs_per_elem[c] > 0)
	      {
		const unsigned int current_size =
		  elem_dofs[c][e].size();
		
		const unsigned int dofs_per_elem =
		  FEInterface::n_dofs_per_elem(dim, fe_type, type);
		
		if (dofs_per_elem > current_size)
		  {
		    elem_dofs[c][e].resize(dofs_per_elem);
		    
		    for (unsigned int i=0; i<dofs_per_elem; i++)
		      elem_dofs[c][e][i] = invalid_number;
		  };
	      };
	  };
    };
   
   
#else
  
  // initialize the node_id_map, each entry should get assigned
  // invalid_number
  node_id_map.resize(n_nodes*n_components());

  // Original Code
  // for (unsigned int i=0; i<node_id_map.size(); i++)
  //     node_id_map[i] = invalid_number;


  // New code 
  std::fill(node_id_map.begin(),
	    node_id_map.end(),
	    invalid_number);




  // initialize the elem_id_map, each entry should get assigned
  // invalid_number
  elem_id_map.resize(n_elem*n_components());

  // Original Code
  // for (unsigned int i=0; i<elem_id_map.size(); i++)
  //     elem_id_map[i] = invalid_number;


  // New code 
  std::fill(elem_id_map.begin(),
	    elem_id_map.end(),
	    invalid_number);

#endif
  
  perf_log.stop_event("reinit()");
};



void DofMap::clear()
{
  // we don't want to clear
  // the coupling matrix!
  // It should not change...
  //dof_coupling.clear();
  
#ifdef ENABLE_EXPENSIVE_DATA_STRUCTURES

  node_dofs.clear();
  elem_dofs.clear();
  
#else
  
  node_id_map.clear();
  elem_id_map.clear();

#endif
  
  first_df.clear();
  last_df.clear();
  send_list.clear();
  n_nz.clear();
  n_oz.clear();


#ifdef ENABLE_AMR

  dof_constraints.clear();

#endif
  
  n_dfs =
    n_nodes =
    n_elem =
    0;
};




void DofMap::distribute_dofs(const unsigned int proc_id)
{
  assert (n_components());
  assert (proc_id < mesh.n_processors());
  
  // re-init in case the underlying mesh has changed
  reinit();

  assert (n_nodes);
  assert (n_elem);

  perf_log.start_event("distribute_dofs()");

  unsigned int next_free_dof=0;

  first_df.resize(mesh.n_processors());
  last_df.resize(mesh.n_processors());

  for (unsigned int s=0; s<mesh.n_processors(); ++s)
    first_df[s] = last_df[s] = 0;

  send_list.clear();
  
  //------------------------------------------------------------
  // DOF numbering
  for (unsigned int processor=0; processor<mesh.n_processors();
       ++processor)
    {
      first_df[processor] = next_free_dof;
      
      for (unsigned comp=0; comp<n_components(); comp++)
	{
	  const FEType fe_type = component_type(comp);
	  
	  for (unsigned int e=0; e<n_elem; ++e)
	    if (mesh.elem(e)->active())  
	      if (mesh.elem(e)->processor_id() == processor)  
		{   // Only number dofs connected to active
		    // elements on this processor.
		  
		  const Elem* elem           = mesh.elem(e);
		  const unsigned int n_nodes = elem->n_nodes();
		  const ElemType type        = elem->type();
		  
		  // First number the nodal DOFS
		  for (unsigned int node=0; node<n_nodes; node++)
		    for (unsigned int index=0;
			 index<FEInterface::n_dofs_at_node(dim, fe_type, type, node);
			 index++)
		      {
			// Global node index
			const unsigned int ig = elem->node(node);
			
			// only assign a dof number if there isn't one
			// already there
			if (node_dof_number(ig,comp,index) == invalid_number)
			  {
			    set_node_dof_number(ig,comp,index) =
			      next_free_dof++;
			  
			    if (processor == proc_id)
			      send_list.push_back(node_dof_number(ig,comp,index));
			  }
			// if there is an entry there and it isn't on this
			// processor it also needs to be added to the send list
			else if (node_dof_number(ig,comp,index) <
				 first_df[proc_id])
			  {
			    if (processor == proc_id)
			      send_list.push_back(node_dof_number(ig,comp,index));
			  };
		      };

		  // Now number the element DOFS
		  const unsigned int ig = e;

		  for (unsigned int index=0;
		       index<FEInterface::n_dofs_per_elem(dim, fe_type, type); index++)
		    {
		  
		      // No way we could have already numbered this DOF!
		      assert (elem_dof_number(ig,comp,index) ==
			      invalid_number);
		  
		      set_elem_dof_number(ig,comp,index) =
			next_free_dof++;
		  
		      if (processor == proc_id)
			send_list.push_back(elem_dof_number(ig,comp,index));
		    };
		};
	};
      
      last_df[processor] = (next_free_dof-1);
    };
  
  n_dfs = next_free_dof;

  
  // Note that in the code above nodes on processor boundaries
  // that are shared by multiple elements are added for each element.
  // Here we need to clean up that data structure
  {    
    std::sort(send_list.begin(), send_list.end());
    std::vector<unsigned int> new_send_list;

    new_send_list.push_back(send_list[0]);
    
    for (unsigned int i=1; i<send_list.size(); i++)
      if (send_list[i] != new_send_list.back())
	new_send_list.push_back(send_list[i]);

    send_list = new_send_list;
  };
  
  perf_log.stop_event("distribute_dofs()");
    
  return;
};



void DofMap::compute_sparsity(const unsigned int proc_id)
{
  assert (n_components());

  perf_log.start_event("compute_sparsity()");
  
  /**
   * Compute the sparsity structure of the global matrix.  This can be
   * fed into a PetscMatrix to allocate exacly the number of nonzeros
   * necessary to store the matrix.  This algorithm should be linear
   * in the (# of elements)*(# nodes per element)
   */
  const unsigned int n_dofs_on_proc    = n_dofs_on_processor(proc_id);
  const unsigned int first_dof_on_proc = first_dof(proc_id);
  const unsigned int last_dof_on_proc  = last_dof(proc_id);
  
  std::vector<std::set<unsigned int> > sparsity_pattern (n_dofs_on_proc);


  /**
   * If the user did not explicitly specify the DOF coupling
   * then all the DOFS are coupled to each other.  Furthermore,
   * we can take a shortcut and to this more quickly here.  So
   * we use an if-test.
   */  
  if (dof_coupling.empty())
    {
      std::vector<unsigned int> element_dofs;
      
      for (unsigned int e=0; e<n_elem; e++)
	if (mesh.elem(e)->active())
	  {
	    dof_indices (e, element_dofs);
	    find_connected_dofs (element_dofs);
	    
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
		    //assert ((ig - first_dof_on_proc) >= 0);
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
		      };
		  };
	      };
	  };      
    } 


  
  /**
   * This is what we do in the case that the user has specified
   * explicit DOF coupling.
   */  
  else
    {
      assert (dof_coupling.size() ==
	      n_components());
      
      const unsigned int n_comp = n_components();
      
      std::vector<unsigned int> element_dofs_i;
      std::vector<unsigned int> element_dofs_j;
      
      
      for (unsigned int e=0; e<n_elem; e++)
	if (mesh.elem(e)->active())
	  for (unsigned int ci=0; ci<n_comp; ci++)
	    {
	      // Find element dofs for component ci
	      dof_indices (e, element_dofs_i, ci);
	      find_connected_dofs (element_dofs_i);
	      
	      const unsigned int n_dofs_on_element_i = element_dofs_i.size();

	      for (unsigned int cj=0; cj<n_comp; cj++)
		if (dof_coupling(ci,cj)) // If ci couples to cj
		  {
		    // Find element dofs for component cj
		    dof_indices (e, element_dofs_j, cj);
		    find_connected_dofs (element_dofs_j);	    
		    
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
			    assert ((ig - first_dof_on_proc) < sparsity_pattern.size());
		    
			    std::set<unsigned int>& row =
			      sparsity_pattern[ig - first_dof_on_proc];
		    
			    for (unsigned int j=0; j<n_dofs_on_element_j; j++)
			      {
				const unsigned int jg = element_dofs_j[j];
			
				row.insert(jg);
			      };
			  };
		      };
		  };
	    };
    };      



  // Now the full sparsity structure is built for all of the
  // DOFs connected to our rows of the matrix.
  n_nz.resize (n_dofs_on_proc);
  n_oz.resize (n_dofs_on_proc);
  
  for (unsigned int i=0; i<n_dofs_on_proc; i++)
    {
      // First zero the counters.
      n_nz[i] = n_oz[i] = 0;
      
      // Then get the row of the sparsity pattern
      const std::set<unsigned int>& row = sparsity_pattern[i];
      
      // Now loop over the row and increment the appropriate
      // counters.
      std::set<unsigned int>::const_iterator it_end = row.end();
      
      for (std::set<unsigned int>::const_iterator it=row.begin();
	   it != it_end; ++it)
	if ((*it >= first_dof_on_proc) &&
	    (*it <= last_dof_on_proc))
	  n_nz[i]++;
	else
	  n_oz[i]++;
    };


  // Explicity clear the sparsity pattern now.
  // it is about to go out of scope, but maybe this
  // will expedite freeing up its memory?
  sparsity_pattern.clear();
  
  perf_log.stop_event("compute_sparsity()");
};



void DofMap::dof_indices (const unsigned int e,
			  std::vector<unsigned int>& di,
			  const unsigned int cn) const
{
  assert (e < mesh.n_elem());

  const Elem* elem = mesh.elem(e);
  
  assert (elem != NULL);


  const unsigned int n_nodes = elem->n_nodes();
  const ElemType type        = elem->type();
  
  // Clear the DOF indices vector.
  di.clear();
  
  // Get the dof numbers
  for (unsigned int c=0; c<n_components(); c++)
    if ((c == cn) || (cn == invalid_number))
      { // Do this for all the components if one was not specified
	// or just for the specified component
	const FEType fe_type = component_type(c);
	
	// Get the node-based DOF numbers
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    const unsigned int global_node = elem->node(n);
	    
	    for (unsigned int i=0;
		 i<FEInterface::n_dofs_at_node(dim, fe_type, type, n); i++)
	      {
		assert (node_dof_number(global_node, c, i) !=
			invalid_number);
		
		di.push_back(node_dof_number(global_node, c, i));
	      };
	  };
	
	// Get the element-based DOF numbers	  
	for (unsigned int i=0; i<FEInterface::n_dofs_per_elem(dim, fe_type, type); i++)
	  {
	    assert (elem_dof_number(e, c, i) !=
		    invalid_number);
	    
	    di.push_back(elem_dof_number(e, c, i));
	  };
      };
};








#ifdef ENABLE_AMR


void DofMap::create_dof_constraints()
{
  perf_log.start_event("create_dof_constraints()");


  
  // Here we build the hanging node constraints.  This is done
  // by enforcing the condition u_a = u_b along hanging sides.
  // u_a = u_b is collocated at the nodes of side a, which gives
  // one row of the constraint matrix.
  
  // clear any existing constraints.
  dof_constraints.clear();
  
  // Look at the element faces.  Check to see if we need to 
  // build constraints.
  for (unsigned int e=0; e<n_elem; e++)
    if (mesh.elem(e)->active())
      for (unsigned int s=0; s<mesh.elem(e)->n_sides(); s++)
	if (mesh.elem(e)->neighbor(s) != NULL)
	  if (mesh.elem(e)->neighbor(s)->level() < mesh.elem(e)->level())  // constrain dofs shared between
	    {                                                              // this element and ones coarser
	                                                                   // than this element.
	      // Get pointers to the two elements of interest.  
#ifndef __IBMCPP__
	      std::auto_ptr<Elem> my_side = mesh.elem(e)->build_side(s); 
#else
	      std::auto_ptr<Elem> my_side(mesh.elem(e)->build_side(s)); 
#endif
	      const Elem* neighbor(mesh.elem(e)->neighbor(s));

	      // Look at all the components in the system
	      for (unsigned int component=0; component<n_components();
		   ++component)
		{
		  const FEType fe_type = component_type(component);
	      
		  for (unsigned int my_dof=0;
		       my_dof<FEInterface::n_dofs(dim-1, fe_type, my_side->type());
		       my_dof++)
		    {
		      assert (my_dof < my_side->n_nodes());

		      // My global node and dof indices.
		      const unsigned int my_node_g = my_side->node(my_dof);
		      const unsigned int my_dof_g  = node_dof_number(my_node_g, component);

		      // Figure out where my node lies on their reference element.
		      const Point mapped_point = FEInterface::inverse_map(dim, fe_type, mesh,
									  neighbor,
									  mesh.vertex(my_node_g));

		      // Compute the neighbor shape function values.
		      for (unsigned int their_dof=0;
			   their_dof<FEInterface::n_dofs(dim, fe_type, neighbor->type());
			   their_dof++)
			{
			  assert (their_dof < neighbor->n_nodes());
			  
			  // Their global node and dof indices.
			  const unsigned int their_node_g = neighbor->node(their_dof);
			  const unsigned int their_dof_g  = node_dof_number(their_node_g, component);

			  const real their_dof_value = FEInterface::shape(dim,
									  fe_type,
									  neighbor->type(),
									  their_dof,
									  mapped_point);

			  // Only add non-zero and non-identity values
			  // for Lagrange basis functions.
			  if ((fabs(their_dof_value) > 1.e-5) &&
			      (fabs(their_dof_value) < .999))
			    {
			      // A reference to the constraint row.
			      DofConstraintRow& constraint_row = dof_constraints[my_dof_g];
			      			      
			      const std::pair<unsigned int, float> p (their_dof_g,
								      static_cast<float>(their_dof_value));

			      constraint_row.insert(p);
			    };
			};		      
		    };
		};
	    };
  
  perf_log.stop_event("create_dof_constraints()");
};



void DofMap::add_constraint_row (const unsigned int dof_number,
				 const DofConstraintRow& constraint_row)
{

  if (!dof_constraints.count(dof_number))
    std::cerr << "WARNING: DOF " << dof_number << " was already constrained!"
	      << std::endl;

  std::pair<unsigned int, DofConstraintRow> kv(dof_number, constraint_row);

  dof_constraints.insert(kv);
};



void DofMap::print_dof_constraints() const
{
  std::cout << "DOF CONSTRAINTS OUTPUT:"
	    << std::endl;
  
  for (DofConstraints::const_iterator it=dof_constraints.begin();
       it != dof_constraints.end(); ++it)
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
    };
};



void DofMap::constrain_element_matrix (DenseMatrix& matrix,
				       std::vector<unsigned int>& elem_dofs) const
{
  assert (elem_dofs.size() == matrix.m());
  assert (elem_dofs.size() == matrix.n());
  
  // The constrained matrix is built up as C^T K C.    
  DenseMatrix C;

  
  build_constraint_matrix (C, elem_dofs);

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
	if (is_constrained_dof(elem_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;
	    
	    // If the DOF is constrained
	    matrix(i,i) = 1.;
	    
	    DofConstraints::const_iterator
	      pos = dof_constraints.find(elem_dofs[i]);
	    
	    assert (pos != dof_constraints.end());
	    
	    const DofConstraintRow& constraint_row = pos->second;
	    
	    assert (!constraint_row.empty());
	    
	    for (DofConstraintRow::const_iterator
		   it=constraint_row.begin(); it != constraint_row.end();
		 ++it)
	      for (unsigned int j=0; j<elem_dofs.size(); j++)
		if (elem_dofs[j] == it->first)
		  matrix(i,j) = -it->second;	
	  };
    }; // end if is constrained...
};



void DofMap::constrain_element_matrix_and_vector (DenseMatrix& matrix,
						  std::vector<real>& rhs,
						  std::vector<unsigned int>& elem_dofs) const
{
  assert (elem_dofs.size() == matrix.m());
  assert (elem_dofs.size() == matrix.n());
  assert (elem_dofs.size() == rhs.size());
  
  // The constrained matrix is built up as C^T K C.
  // The constrained RHS is built up as C^T F
  DenseMatrix C;
  
  build_constraint_matrix (C, elem_dofs);

  
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
	if (is_constrained_dof(elem_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;
	    
	    // If the DOF is constrained
	    matrix(i,i) = 1.;
	    
	    DofConstraints::const_iterator
	      pos = dof_constraints.find(elem_dofs[i]);
	    
	    assert (pos != dof_constraints.end());
	    
	    const DofConstraintRow& constraint_row = pos->second;
	    
	    assert (!constraint_row.empty());
	    
	    for (DofConstraintRow::const_iterator
		   it=constraint_row.begin(); it != constraint_row.end();
		 ++it)
	      for (unsigned int j=0; j<elem_dofs.size(); j++)
		if (elem_dofs[j] == it->first)
		  matrix(i,j) = -it->second;	
	  };

      
      // Compute the matrix-vector product C^T F
      std::vector<real> old_rhs(rhs);
      
      rhs.resize(elem_dofs.size());
      
      
      for (unsigned int i=0; i<elem_dofs.size(); i++)
	{
	  // zero before summation
	  rhs[i] = 0.;
	  
	  for (unsigned int j=0; j<old_rhs.size(); j++)
	    rhs[i] += C.transpose(i,j)*old_rhs[j];
	};


      assert (elem_dofs.size() == rhs.size());

      for (unsigned int i=0; i<elem_dofs.size(); i++)
	if (is_constrained_dof(elem_dofs[i]))
	  {	
	    // If the DOF is constrained
	    rhs[i] = 0.;
	  };
    }; // end if is constrained...
};



void DofMap::constrain_element_matrix (DenseMatrix& matrix,
				       std::vector<unsigned int>& row_dofs,
				       std::vector<unsigned int>& col_dofs) const
{
  assert (row_dofs.size() == matrix.m());
  assert (col_dofs.size() == matrix.n());
  
  // The constrained matrix is built up as R^T K C.


  
  
  DenseMatrix R;
  DenseMatrix C;

  // Safeguard against the user passing us the same
  // object for row_dofs and col_dofs.  If that is done
  // the calls to build_matrix would fail
  std::vector<unsigned int> orig_row_dofs(row_dofs);
  std::vector<unsigned int> orig_col_dofs(col_dofs);
  
  build_constraint_matrix (R, orig_row_dofs);
  build_constraint_matrix (C, orig_col_dofs);

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
	if (is_constrained_dof(row_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;
	    
	    // If the DOF is constrained
	    matrix(i,i) = 1.;
	    
	    DofConstraints::const_iterator
	      pos = dof_constraints.find(row_dofs[i]);
	    
	    assert (pos != dof_constraints.end());
	    
	    const DofConstraintRow& constraint_row = pos->second;
	    
	    assert (!constraint_row.empty());
	    
	    for (DofConstraintRow::const_iterator
		   it=constraint_row.begin(); it != constraint_row.end();
		 ++it)
	      for (unsigned int j=0; j<col_dofs.size(); j++)
		if (col_dofs[j] == it->first)
		  matrix(i,j) = -it->second;	
	  };
    }; // end if is constrained...
};



void DofMap::constrain_element_vector (std::vector<real>&         rhs,
				       std::vector<unsigned int>& row_dofs) const
{
  assert (rhs.size() == row_dofs.size());

  
  // The constrained RHS is built up as R^T F.  
  DenseMatrix R;

  build_constraint_matrix (R, row_dofs);

  
  // It is possible that the vector is not constrained at all.
  if ((R.m() == rhs.size()) &&
      (R.n() == row_dofs.size())) // if the RHS is constrained
    {
      // Compute the matrix-vector product
      std::vector<real> old_rhs(rhs);
      
      rhs.resize(row_dofs.size());
      
      
      for (unsigned int i=0; i<row_dofs.size(); i++)
	{
	  // zero before summation
	  rhs[i] = 0.;
	  
	  for (unsigned int j=0; j<old_rhs.size(); j++)
	    rhs[i] += R.transpose(i,j)*old_rhs[j];
	};


      assert (row_dofs.size() == rhs.size());

      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (is_constrained_dof(row_dofs[i]))
	  {	
	    // If the DOF is constrained
	    rhs[i] = 0.;
	  };
    }; // end if the RHS is constrained.
};



void DofMap::build_constraint_matrix (DenseMatrix& C,
				      std::vector<unsigned int>& elem_dofs) const
{
  typedef std::set<unsigned int> RCSet;
  RCSet dof_set;

  bool done = true;

  // First insert the DOFS we already depend on into the sets.
  {
    for (unsigned int i=0; i<elem_dofs.size(); i++)
      dof_set.insert(elem_dofs[i]);
  };
  

  // Next insert any dofs those might be constrained in terms
  // of.  Note that in this case we may not be done:  Those may
  // in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  {
    const unsigned int orig_dof_set_size = dof_set.size();
    
    for (unsigned int i=0; i<elem_dofs.size(); i++)
      if (is_constrained_dof(elem_dofs[i]))
	{
	  // If the DOF is constrained
	  DofConstraints::const_iterator
	    pos = dof_constraints.find(elem_dofs[i]);
	  
	  assert (pos != dof_constraints.end());
	  
	  const DofConstraintRow& constraint_row = pos->second;
	  
	  assert (!constraint_row.empty());
	  
	  for (DofConstraintRow::const_iterator
		 it=constraint_row.begin(); it != constraint_row.end();
	       ++it)
	    dof_set.insert (it->first);
	};

    // If we added any DOFS then we need to do this recursively.
    // It is possible that we just added a DOF that is also
    // constrained!
    if (dof_set.size() != orig_dof_set_size)
      done = false;
  };


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
      };
      
      // Now we can build the constraint matrix.
      // Note that resize also zeros for a DenseMatrix.
      C.resize (elem_dofs.size(), new_elem_dofs.size());
      
      // Create the C constraint matrix.
      {
	for (unsigned int i=0; i<elem_dofs.size(); i++)
	  if (is_constrained_dof(elem_dofs[i]))
	    {
	      // If the DOF is constrained
	      DofConstraints::const_iterator
		pos = dof_constraints.find(elem_dofs[i]);
	      
	      assert (pos != dof_constraints.end());
	      
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
	    };	
	//C.print();
      };

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      elem_dofs = new_elem_dofs;
      
      DenseMatrix Cnew;
      
      build_constraint_matrix (Cnew, elem_dofs);

      if ((C.n() == Cnew.m()) &&
	  (Cnew.n() == elem_dofs.size())) // If the constraint matrix
	{                                 // is constrained...
	  //here();
	  C.right_multiply(Cnew, false);
	};
      
      assert (C.n() == elem_dofs.size());
    }; // end if (!done)
};


// endif of ENABLE_AMR
#endif










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
      if (is_constrained_dof(elem_dofs[i]))
	{
	  // If the DOF is constrained
	  DofConstraints::const_iterator
	    pos = dof_constraints.find(elem_dofs[i]);
	  
	  assert (pos != dof_constraints.end());
	  
	  const DofConstraintRow& constraint_row = pos->second;
	  
	  assert (!constraint_row.empty());
	  
	  for (DofConstraintRow::const_iterator
		 it=constraint_row.begin(); it != constraint_row.end();
	       ++it)
	    dof_set.insert (it->first);
	};

    // If we added any DOFS then we need to do this recursively.
    // It is possible that we just added a DOF that is also
    // constrained!
    if (dof_set.size() != orig_dof_set_size)
      done = false;
  };


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
      };
      

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      find_connected_dofs (elem_dofs);
      
    }; // end if (!done)


#endif

};
