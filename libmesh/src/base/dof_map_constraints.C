// $Id: dof_map_constraints.C,v 1.8 2004-11-03 21:30:19 benkirk Exp $

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



// C++ Includes -------------------------------------
#include <set>
#include <algorithm> // for std::count, std::fill

// Local Includes -----------------------------------
#include "dof_map.h"
#include "elem.h"
#include "mesh_base.h"
#include "fe_interface.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "libmesh.h"
#include "libmesh_logging.h"




// ------------------------------------------------------------
// DofMap member functions

#ifdef ENABLE_AMR

void DofMap::create_dof_constraints(const MeshBase& mesh)
{
  START_LOG("create_dof_constraints()", "DofMap");

  assert (mesh.is_prepared());
  
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
  
  // Look at all the variables in the system
  for (unsigned int variable_number=0; variable_number<this->n_variables();
       ++variable_number)
    {
      const FEType& fe_type = this->variable_type(variable_number);
      
      const_elem_iterator       elem_it (mesh.elements_begin());
      const const_elem_iterator elem_end(mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	FEInterface::compute_constraints (_dof_constraints,
					  this->sys_number(),
					  variable_number,
					  fe_type,
					  *elem_it);
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

  START_LOG("constrain_elem_matrix()", "DofMap");
  
  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // Compute the matrix-matrix-matrix product C^T K C
      matrix.left_multiply_transpose  (C);
      matrix.right_multiply (C);
      
      
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
  
  STOP_LOG("constrain_elem_matrix()", "DofMap");  
}



void DofMap::constrain_element_matrix_and_vector (DenseMatrix<Number>& matrix,
						  DenseVector<Number>& rhs,
						  std::vector<unsigned int>& elem_dofs) const
{
  assert (elem_dofs.size() == matrix.m());
  assert (elem_dofs.size() == matrix.n());
  assert (elem_dofs.size() == rhs.size());
  
  // The constrained matrix is built up as C^T K C.
  // The constrained RHS is built up as C^T F
  DenseMatrix<Number> C;
  
  this->build_constraint_matrix (C, elem_dofs);
  
  START_LOG("cnstrn_elem_mat_vec()", "DofMap");
  
  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // Compute the matrix-matrix-matrix product C^T K C
      matrix.left_multiply_transpose  (C);
      matrix.right_multiply (C);
      
      
      assert (matrix.m() == matrix.n());
      assert (matrix.m() == elem_dofs.size());
      assert (matrix.n() == elem_dofs.size());
      

      // This will put a nonsymmetric entry in the constraint
      // row to ensure that the linear system produces the
      // correct value for the constrained DOF.
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
      DenseVector<Number> old_rhs(rhs);

      // resize the RHS vector & 0 before summation
      rhs.resize(elem_dofs.size());
      rhs.zero();
      
      // compute matrix/vector product
      for (unsigned int i=0; i<elem_dofs.size(); i++)
	for (unsigned int j=0; j<old_rhs.size(); j++)
	  rhs(i) += C.transpose(i,j)*old_rhs(j);
	

      assert (elem_dofs.size() == rhs.size());

      for (unsigned int i=0; i<elem_dofs.size(); i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {	
	    // If the DOF is constrained
	    rhs(i) = 0.;
	  }
    } // end if is constrained...
  
  STOP_LOG("cnstrn_elem_mat_vec()", "DofMap");  
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

  START_LOG("constrain_elem_matrix()", "DofMap");
  
  row_dofs = orig_row_dofs;
  col_dofs = orig_col_dofs;
  
    
  // It is possible that the matrix is not constrained at all.
  if ((R.m() == matrix.m()) &&
      (R.n() == row_dofs.size()) &&
      (C.m() == matrix.n()) &&
      (C.n() == col_dofs.size())) // If the matrix is constrained
    {
      // K_constrained = R^T K C
      matrix.left_multiply_transpose  (R);
      matrix.right_multiply (C);
      
      
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
  
  STOP_LOG("constrain_elem_matrix()", "DofMap");  
}



void DofMap::constrain_element_vector (DenseVector<Number>&       rhs,
				       std::vector<unsigned int>& row_dofs) const
{
  assert (rhs.size() == row_dofs.size());
  
  // The constrained RHS is built up as R^T F.  
  DenseMatrix<Number> R;

  this->build_constraint_matrix (R, row_dofs);

  START_LOG("constrain_elem_vector()", "DofMap");
    
  // It is possible that the vector is not constrained at all.
  if ((R.m() == rhs.size()) &&
      (R.n() == row_dofs.size())) // if the RHS is constrained
    {
      // Compute the matrix-vector product
      DenseVector<Number> old_rhs(rhs);

      // resize RHS & zero before summation
      rhs.resize(row_dofs.size());
      rhs.zero();

      // compute matrix/vector product
      for (unsigned int i=0; i<row_dofs.size(); i++)
	for (unsigned int j=0; j<old_rhs.size(); j++)
	  rhs(i) += R.transpose(i,j)*old_rhs(j);
      
      assert (row_dofs.size() == rhs.size());

      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (this->is_constrained_dof(row_dofs[i]))
	  {	
	    // If the DOF is constrained
	    rhs(i) = 0.;
	  }
    } // end if the RHS is constrained.
  
  STOP_LOG("constrain_elem_vector()", "DofMap");  
}



void DofMap::build_constraint_matrix (DenseMatrix<Number>& C,
				      std::vector<unsigned int>& elem_dofs) const
{
  START_LOG("build_constraint_matrix()", "DofMap");

  // Create a set containing the DOFs we already depend on
  typedef std::set<unsigned int> RCSet;
  RCSet dof_set (elem_dofs.begin(),
		 elem_dofs.end());

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
	
	for (DofConstraintRow::const_iterator
	       it=constraint_row.begin(); it != constraint_row.end();
	     ++it)
	  dof_set.insert (it->first);
      }

  // If we added any DOFS then we need to do this recursively.
  // It is possible that we just added a DOF that is also
  // constrained!
  if (dof_set.size() != elem_dofs.size())
    {
      // Create a new list of element DOFs containing the
      // contents of the current dof_set.
      std::vector<unsigned int> new_elem_dofs (dof_set.begin(),
					       dof_set.end());
      
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
      
      STOP_LOG("build_constraint_matrix()", "DofMap");
      
      this->build_constraint_matrix (Cnew, elem_dofs);
      
      START_LOG("build_constraint_matrix()", "DofMap");  

      if ((C.n() == Cnew.m()) &&
	  (Cnew.n() == elem_dofs.size())) // If the constraint matrix
	{                                 // is constrained...
	  //here();
	  //C.right_multiply(Cnew, false);
	  C.right_multiply(Cnew);
	}
      
      assert (C.n() == elem_dofs.size());
    } // end if (!done)
  
  STOP_LOG("build_constraint_matrix()", "DofMap");  
}

#endif // #ifdef ENABLE_AMR
