// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "fe_type.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "libmesh_logging.h"
#include "system.h" // needed by enforce_constraints_exactly()
#include "mesh.h"   // as is this
#include "numeric_vector.h" // likewise
#include "parallel.h"
#include "point_locator_base.h"
#include "elem_range.h"
#include "threads_allocators.h"



// Anonymous namespace to hold helper classes
namespace {

  class ComputeConstraints
  {
  public:
    ComputeConstraints (DofConstraints &constraints,
			DofMap &dof_map,
#ifdef LIBMESH_ENABLE_PERIODIC
			PeriodicBoundaries &periodic_boundaries,
#endif
			const MeshBase &mesh,
			const unsigned int variable_number) :
      _constraints(constraints),
      _dof_map(dof_map),
#ifdef LIBMESH_ENABLE_PERIODIC
      _periodic_boundaries(periodic_boundaries),
#endif
      _mesh(mesh),
      _variable_number(variable_number)
    {}

    void operator()(const ConstElemRange &range) const
    {
      const System::Variable &var_description = _dof_map.variable(_variable_number);

      for (ConstElemRange::const_iterator it = range.begin(); it!=range.end(); ++it)
	if (var_description.active_on_subdomain((*it)->subdomain_id()))
	  {
#ifdef LIBMESH_ENABLE_AMR
	    FEInterface::compute_constraints (_constraints,
					      _dof_map,
					      _variable_number,
					      *it);
#endif
#ifdef LIBMESH_ENABLE_PERIODIC
	    // FIXME: periodic constraints won't work on a non-serial
	    // mesh unless it's kept ghost elements from opposing
	    // boundaries!
	    FEInterface::compute_periodic_constraints (_constraints,
						       _dof_map,
						       _periodic_boundaries,
						       _mesh,
						       _variable_number,
						       *it);
#endif
	  }
    }
    
  private:
    DofConstraints &_constraints;
    DofMap &_dof_map;
#ifdef LIBMESH_ENABLE_PERIODIC
    PeriodicBoundaries &_periodic_boundaries;
#endif
    const MeshBase &_mesh;
    const unsigned int _variable_number;    
  };
}



// ------------------------------------------------------------
// DofMap member functions

#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)

void DofMap::create_dof_constraints(const MeshBase& mesh)
{
  START_LOG("create_dof_constraints()", "DofMap");

  libmesh_assert (mesh.is_prepared());
  
  const unsigned int dim = mesh.mesh_dimension();

  // Constraints are not necessary in 1D
  if (dim == 1)
  {
    // make sure we stop logging though
    STOP_LOG("create_dof_constraints()", "DofMap");
    return;
  }
  
  // Here we build the hanging node constraints.  This is done
  // by enforcing the condition u_a = u_b along hanging sides.
  // u_a = u_b is collocated at the nodes of side a, which gives
  // one row of the constraint matrix.
  
  // clear any existing constraints.
  _dof_constraints.clear();

  // define the range of elements of interest
  ConstElemRange range;
  {
    // With SerialMesh or a serial ParallelMesh, every processor
    // computes every constraint
    MeshBase::const_element_iterator
      elem_begin = mesh.elements_begin(),
      elem_end   = mesh.elements_end();
    
    // With a parallel ParallelMesh, processors compute only
    // their local constraints
    if (!mesh.is_serial())
      {
	elem_begin = mesh.local_elements_begin();
	elem_end   = mesh.local_elements_end(); 
      }

    // set the range to contain the specified elements
    range.reset (elem_begin, elem_end);
  }
  
  // Look at all the variables in the system.  Reset the element
  // range at each iteration -- there is no need to reconstruct it.
  for (unsigned int variable_number=0; variable_number<this->n_variables();
       ++variable_number, range.reset())
    Threads::parallel_for (range,
			   ComputeConstraints (_dof_constraints,
					       *this,
#ifdef LIBMESH_ENABLE_PERIODIC
					       _periodic_boundaries,
#endif
					       mesh,
					       variable_number));

  // With a parallelized Mesh, we've computed our local constraints,
  // but they may depend on non-local constraints that we'll need to
  // take into account.
  if (!mesh.is_serial())
    this->allgather_recursive_constraints();
  
  STOP_LOG("create_dof_constraints()", "DofMap");
}



void DofMap::add_constraint_row (const unsigned int dof_number,
				 const DofConstraintRow& constraint_row,
				 const bool forbid_constraint_overwrite)
{
  // Optionally allow the user to overwrite constraints.  Defaults to false.
  if (forbid_constraint_overwrite)
    if (this->is_constrained_dof(dof_number))
      {
	std::cerr << "ERROR: DOF " << dof_number << " was already constrained!"
		  << std::endl;
	libmesh_error();
      }
  
  std::pair<unsigned int, DofConstraintRow> kv(dof_number, constraint_row);

  _dof_constraints.insert(kv);
}



void DofMap::print_dof_constraints(std::ostream& os) const
{
  os << "DOF CONSTRAINTS OUTPUT:"
     << std::endl;
  
  for (DofConstraints::const_iterator it=_dof_constraints.begin();
       it != _dof_constraints.end(); ++it)
    {
      const unsigned int i = it->first;
      const DofConstraintRow& row = it->second;

      os << "Constraints for DOF " << i
	 << ": \t";

      for (DofConstraintRow::const_iterator pos=row.begin();
	   pos != row.end(); ++pos)
	os << " (" << pos->first << ","
	   << pos->second << ")\t";

      os << std::endl;
    }
}



void DofMap::constrain_element_matrix (DenseMatrix<Number>& matrix,
				       std::vector<unsigned int>& elem_dofs,
				       bool asymmetric_constraint_rows) const
{
  libmesh_assert (elem_dofs.size() == matrix.m());
  libmesh_assert (elem_dofs.size() == matrix.n());

  // check for easy return
  if (this->n_constrained_dofs() == 0)
    return;
    
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
      
      
      libmesh_assert (matrix.m() == matrix.n());
      libmesh_assert (matrix.m() == elem_dofs.size());
      libmesh_assert (matrix.n() == elem_dofs.size());
      
      
      for (unsigned int i=0; i<elem_dofs.size(); i++)
	// If the DOF is constrained
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;
	    
	    matrix(i,i) = 1.;
	    
            if (asymmetric_constraint_rows)
              {
	        DofConstraints::const_iterator
	          pos = _dof_constraints.find(elem_dofs[i]);
	    
	        libmesh_assert (pos != _dof_constraints.end());
	    
	        const DofConstraintRow& constraint_row = pos->second;
	    
	        libmesh_assert (!constraint_row.empty());
	    
	        for (DofConstraintRow::const_iterator
		       it=constraint_row.begin(); it != constraint_row.end();
		     ++it)
	          for (unsigned int j=0; j<elem_dofs.size(); j++)
		    if (elem_dofs[j] == it->first)
		      matrix(i,j) = -it->second;	
	      }
	  }
    } // end if is constrained...
  
  STOP_LOG("constrain_elem_matrix()", "DofMap");  
}



void DofMap::constrain_element_matrix_and_vector (DenseMatrix<Number>& matrix,
						  DenseVector<Number>& rhs,
						  std::vector<unsigned int>& elem_dofs,
						  bool asymmetric_constraint_rows) const
{
  libmesh_assert (elem_dofs.size() == matrix.m());
  libmesh_assert (elem_dofs.size() == matrix.n());
  libmesh_assert (elem_dofs.size() == rhs.size());

  // check for easy return
  if (this->n_constrained_dofs() == 0)
    return;
  
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
      
      
      libmesh_assert (matrix.m() == matrix.n());
      libmesh_assert (matrix.m() == elem_dofs.size());
      libmesh_assert (matrix.n() == elem_dofs.size());
      

      for (unsigned int i=0; i<elem_dofs.size(); i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;
	    
	    // If the DOF is constrained
	    matrix(i,i) = 1.;
	    
            // This will put a nonsymmetric entry in the constraint
            // row to ensure that the linear system produces the
            // correct value for the constrained DOF.
            if (asymmetric_constraint_rows)
              {
	        DofConstraints::const_iterator
	          pos = _dof_constraints.find(elem_dofs[i]);
	    
	        libmesh_assert (pos != _dof_constraints.end());
	    
	        const DofConstraintRow& constraint_row = pos->second;
	    
// p refinement creates empty constraint rows
//	    libmesh_assert (!constraint_row.empty());
	    
	        for (DofConstraintRow::const_iterator
		       it=constraint_row.begin(); it != constraint_row.end();
		     ++it)
	          for (unsigned int j=0; j<elem_dofs.size(); j++)
		    if (elem_dofs[j] == it->first)
		      matrix(i,j) = -it->second;	
              }
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
	

      libmesh_assert (elem_dofs.size() == rhs.size());

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
				       std::vector<unsigned int>& col_dofs,
				       bool asymmetric_constraint_rows) const
{
  libmesh_assert (row_dofs.size() == matrix.m());
  libmesh_assert (col_dofs.size() == matrix.n());

  // check for easy return
  if (this->n_constrained_dofs() == 0)
    return;
  
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
      
      
      libmesh_assert (matrix.m() == row_dofs.size());
      libmesh_assert (matrix.n() == col_dofs.size());
      
      
      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (this->is_constrained_dof(row_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;
	  
	    // If the DOF is constrained
	    matrix(i,i) = 1.;
	    
            if (asymmetric_constraint_rows)
              {
	        DofConstraints::const_iterator
	          pos = _dof_constraints.find(row_dofs[i]);
	    
	        libmesh_assert (pos != _dof_constraints.end());
	    
	        const DofConstraintRow& constraint_row = pos->second;
	    
	        libmesh_assert (!constraint_row.empty());
	    
	        for (DofConstraintRow::const_iterator
		       it=constraint_row.begin(); it != constraint_row.end();
		     ++it)
	          for (unsigned int j=0; j<col_dofs.size(); j++)
		    if (col_dofs[j] == it->first)
		      matrix(i,j) = -it->second;	
              }
	  }
    } // end if is constrained...
  
  STOP_LOG("constrain_elem_matrix()", "DofMap");  
}



void DofMap::constrain_element_vector (DenseVector<Number>&       rhs,
				       std::vector<unsigned int>& row_dofs,
				       bool) const
{
  libmesh_assert (rhs.size() == row_dofs.size());

  // check for easy return
  if (this->n_constrained_dofs() == 0)
    return;
  
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
      
      libmesh_assert (row_dofs.size() == rhs.size());

      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (this->is_constrained_dof(row_dofs[i]))
	  {	
	    // If the DOF is constrained
	    rhs(i) = 0.;
	  }
    } // end if the RHS is constrained.
  
  STOP_LOG("constrain_elem_vector()", "DofMap");  
}



void DofMap::constrain_element_dyad_matrix (DenseVector<Number>& v,
					    DenseVector<Number>& w,
					    std::vector<unsigned int>& row_dofs,
					    bool) const
{
  libmesh_assert (v.size() == row_dofs.size());
  libmesh_assert (w.size() == row_dofs.size());

  // check for easy return
  if (this->n_constrained_dofs() == 0)
    return;
  
  // The constrained RHS is built up as R^T F.  
  DenseMatrix<Number> R;

  this->build_constraint_matrix (R, row_dofs);

  START_LOG("cnstrn_elem_dyad_mat()", "DofMap");
    
  // It is possible that the vector is not constrained at all.
  if ((R.m() == v.size()) &&
      (R.n() == row_dofs.size())) // if the RHS is constrained
    {
      // Compute the matrix-vector products
      DenseVector<Number> old_v(v);
      DenseVector<Number> old_w(w);

      // resize RHS & zero before summation
      v.resize(row_dofs.size());
      w.resize(row_dofs.size());

      // compute matrix/vector product
      for (unsigned int i=0; i<row_dofs.size(); i++)
	for (unsigned int j=0; j<old_v.size(); j++)
	  {
	    v(i) += R.transpose(i,j)*old_v(j);
	    w(i) += R.transpose(i,j)*old_w(j);
	  }
      
      libmesh_assert (row_dofs.size() == v.size());
      libmesh_assert (row_dofs.size() == w.size());

      /* Constrain only v, not w.  */

      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (this->is_constrained_dof(row_dofs[i]))
	  {	
	    // If the DOF is constrained
	    v(i) = 0.;
	  }
    } // end if the RHS is constrained.
  
  STOP_LOG("cnstrn_elem_dyad_mat()", "DofMap");
}



void DofMap::constrain_nothing (std::vector<unsigned int>& dofs) const
{
  // check for easy return
  if (this->n_constrained_dofs() == 0)
    return;

  // All the work is done by \p build_constraint_matrix.  We just need
  // a dummy matrix.
  DenseMatrix<Number> R;
  this->build_constraint_matrix (R, dofs);
}



void DofMap::enforce_constraints_exactly (const System &system,
                                          NumericVector<Number> *v) const
{
  parallel_only();

  unsigned int local_constraints = this->n_constrained_dofs();
  Parallel::max(local_constraints);
  if (!local_constraints)
    return;
  
  START_LOG("enforce_constraints_exactly()","DofMap");

  if (!v)
    v = system.solution.get();

  NumericVector<Number> *v_local  = NULL; // will be initialized below
  NumericVector<Number> *v_global = NULL; // will be initialized below
  AutoPtr<NumericVector<Number> > v_built;
  if (v->type() == SERIAL)
    {
      v_built = NumericVector<Number>::build();
      v_built->init(this->n_dofs(), this->n_local_dofs(), true, PARALLEL);
      v_built->close();

      for (unsigned int i=v_built->first_local_index();
           i<v_built->last_local_index(); i++)
        v_built->set(i, (*v)(i));
      v_built->close();
      v_global = v_built.get();

      v_local = v;
      libmesh_assert (v_local->closed());
    }
  else if (v->type() == PARALLEL)
    {
      v_built = NumericVector<Number>::build();
      v_built->init (v->size(), v->size(), true, SERIAL);
      v->localize(*v_built);
      v_built->close();
      v_local = v_built.get();

      v_global = v;
    }
  else if (v->type() == GHOSTED)
    {
      v_local = v;
      v_global = v;
    }
  else // unknown v->type()
    {
      std::cerr << "ERROR: Unknown v->type() == " << v->type() 
		<< std::endl;
      libmesh_error();
    }

  // We should never hit these asserts because we should error-out in
  // else clause above.  Just to be sure we don't try to use v_local
  // and v_global uninitialized...
  libmesh_assert(v_local);
  libmesh_assert(v_global);
  libmesh_assert (this == &(system.get_dof_map()));

  DofConstraints::const_iterator c_it = _dof_constraints.begin();
  const DofConstraints::const_iterator c_end = _dof_constraints.end();

  for ( ; c_it != c_end; ++c_it)
    {
      unsigned int constrained_dof = c_it->first;
      if (constrained_dof < this->first_dof() ||
          constrained_dof >= this->end_dof())
        continue;

      const DofConstraintRow constraint_row = c_it->second;

      Number exact_value = 0;
      for (DofConstraintRow::const_iterator
	   j=constraint_row.begin(); j != constraint_row.end();
	   ++j)
        exact_value += j->second * (*v_local)(j->first);

      v_global->set(constrained_dof, exact_value);
    }

  // If the old vector was serial, we probably need to send our values
  // to other processors
  if (v->type() == SERIAL)
    {
      v_global->localize (*v);
    }
  v->close();

  STOP_LOG("enforce_constraints_exactly()","DofMap");
}



std::pair<Real, Real>
DofMap::max_constraint_error (const System &system,
                              NumericVector<Number> *v) const
{
  if (!v)
    v = system.solution.get();
  NumericVector<Number> &vec = *v;

  // We'll assume the vector is closed
  libmesh_assert (vec.closed());

  Real max_absolute_error = 0., max_relative_error = 0.;

  const MeshBase &mesh = system.get_mesh();

  libmesh_assert (this == &(system.get_dof_map()));

  // indices on each element
  std::vector<unsigned int> local_dof_indices;

  MeshBase::const_element_iterator       elem_it  =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator elem_end =
    mesh.active_local_elements_end(); 
      
  for ( ; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;

      this->dof_indices(elem, local_dof_indices);
      std::vector<unsigned int> raw_dof_indices = local_dof_indices;

      // Constraint matrix for each element
      DenseMatrix<Number> C;

      this->build_constraint_matrix (C, local_dof_indices);

      // Continue if the element is unconstrained
      if (!C.m())
        continue;

      libmesh_assert(C.m() == raw_dof_indices.size());
      libmesh_assert(C.n() == local_dof_indices.size());

      for (unsigned int i=0; i!=C.m(); ++i)
        {
          // Recalculate any constrained dof owned by this processor
          unsigned int global_dof = raw_dof_indices[i];
          if (this->is_constrained_dof(global_dof) &&
              global_dof >= vec.first_local_index() &&
              global_dof < vec.last_local_index())
          {
            Number exact_value = 0;
            for (unsigned int j=0; j!=C.n(); ++j)
              {
                if (local_dof_indices[j] != global_dof)
                  exact_value += C(i,j) * 
                    vec(local_dof_indices[j]);
              }

            max_absolute_error = std::max(max_absolute_error,
              std::abs(vec(global_dof) - exact_value));
            max_relative_error = std::max(max_relative_error,
              std::abs(vec(global_dof) - exact_value)
              / std::abs(exact_value));
          }
        }
    }

  return std::pair<Real, Real>(max_absolute_error, max_relative_error);
}



void DofMap::build_constraint_matrix (DenseMatrix<Number>& C,
				      std::vector<unsigned int>& elem_dofs,
				      const bool called_recursively) const
{
  if (!called_recursively) START_LOG("build_constraint_matrix()", "DofMap");

  // Create a set containing the DOFs we already depend on
  typedef std::set<unsigned int> RCSet;
  RCSet dof_set;

  bool we_have_constraints = false;

  // Next insert any other dofs the current dofs might be constrained
  // in terms of.  Note that in this case we may not be done: Those
  // may in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  for (unsigned int i=0; i<elem_dofs.size(); i++)
    if (this->is_constrained_dof(elem_dofs[i]))
      {
        we_have_constraints = true;

	// If the DOF is constrained
	DofConstraints::const_iterator
	  pos = _dof_constraints.find(elem_dofs[i]);
	
	libmesh_assert (pos != _dof_constraints.end());
	
	const DofConstraintRow& constraint_row = pos->second;
	
// Constraint rows in p refinement may be empty
//	libmesh_assert (!constraint_row.empty());
	
	for (DofConstraintRow::const_iterator
	       it=constraint_row.begin(); it != constraint_row.end();
	     ++it)
	  dof_set.insert (it->first);
      }

  // May be safe to return at this point
  // (but remember to stop the perflog)
  if (!we_have_constraints)
    {
      STOP_LOG("build_constraint_matrix()", "DofMap");
      return;
    }

  // delay inserting elem_dofs for efficiency in the case of
  // no constraints.  In that case we don't get here!
  dof_set.insert (elem_dofs.begin(),
		  elem_dofs.end());

  // If we added any DOFS then we need to do this recursively.
  // It is possible that we just added a DOF that is also
  // constrained!
  //
  // Also, we need to handle the special case of an element having DOFs
  // constrained in terms of other, local DOFs
  if ((dof_set.size() != elem_dofs.size()) || // case 1: constrained in terms of other DOFs
      !called_recursively)                    // case 2: constrained in terms of our own DOFs
    {
      // Create a new list of element DOFs containing the
      // contents of the current dof_set.
      std::vector<unsigned int> new_elem_dofs (dof_set.begin(),
					       dof_set.end());
      
      // Now we can build the constraint matrix.
      // Note that resize also zeros for a DenseMatrix<Number>.
      C.resize (elem_dofs.size(), new_elem_dofs.size());
      
      // Create the C constraint matrix.
      for (unsigned int i=0; i<elem_dofs.size(); i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    // If the DOF is constrained
	    DofConstraints::const_iterator
	      pos = _dof_constraints.find(elem_dofs[i]);
	    
	    libmesh_assert (pos != _dof_constraints.end());
	    
	    const DofConstraintRow& constraint_row = pos->second;
	    
// p refinement creates empty constraint rows
//	    libmesh_assert (!constraint_row.empty());
	    
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

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      elem_dofs = new_elem_dofs;
      
      DenseMatrix<Number> Cnew;
      
      this->build_constraint_matrix (Cnew, elem_dofs, true);

      if ((C.n() == Cnew.m()) &&
	  (Cnew.n() == elem_dofs.size())) // If the constraint matrix	                                 
	C.right_multiply(Cnew);           // is constrained...
      
      libmesh_assert (C.n() == elem_dofs.size());
    }
  
  if (!called_recursively) STOP_LOG("build_constraint_matrix()", "DofMap");  
}



void DofMap::allgather_recursive_constraints()
{
  // This function must be run on all processors at once
  parallel_only();

  // Return immediately if there's nothing to gather
  if (libMesh::n_processors() == 1)
    return;

  // We might get to return immediately if none of the processors
  // found any constraints
  unsigned int has_constraints = !_dof_constraints.empty();
  Parallel::max(has_constraints);
  if (!has_constraints)
    return;

  // We might have calculated constraints for constrained dofs
  // which live on other processors.
  // Push these out first.
  {
  std::vector<std::vector<unsigned int> > pushed_ids(libMesh::n_processors());
  std::vector<unsigned int> pushed_on_proc(libMesh::n_processors(), 0);

  // Count the constraints to push to each processor
  unsigned int push_proc_id = 0;
  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    {
      unsigned int constrained = i->first;
      while (constrained >= _end_df[push_proc_id])
        push_proc_id++;
      pushed_on_proc[push_proc_id]++;
    }
  for (unsigned int p = 0; p != libMesh::n_processors(); ++p)
    pushed_ids[p].reserve(pushed_on_proc[p]);

  // Collect the constraints to push to each processor
  push_proc_id = 0;
  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    {
      unsigned int constrained = i->first;
      while (constrained >= _end_df[push_proc_id])
        push_proc_id++;
      pushed_ids[push_proc_id].push_back(constrained);
    }

  // Now trade constraint rows
  for (unsigned int p = 0; p != libMesh::n_processors(); ++p)
    {
      // Push to processor procup while receiving from procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();

      // Pack the constraint rows to push to procup
      std::vector<std::vector<unsigned int> > pushed_keys(pushed_ids[procup].size());
      std::vector<std::vector<Real> > pushed_vals(pushed_ids[procup].size());
      for (unsigned int i = 0; i != pushed_ids[procup].size(); ++i) 
        {
          DofConstraintRow &row = _dof_constraints[pushed_ids[procup][i]];
          unsigned int row_size = row.size();
          pushed_keys[i].reserve(row_size);
          pushed_vals[i].reserve(row_size);
          for (DofConstraintRow::iterator j = row.begin();
               j != row.end(); ++j)
            {
              pushed_keys[i].push_back(j->first);
              pushed_vals[i].push_back(j->second);
            }
        }

      // Trade pushed constraint rows
      std::vector<unsigned int> pushed_ids_to_me;
      std::vector<std::vector<unsigned int> > pushed_keys_to_me;
      std::vector<std::vector<Real> > pushed_vals_to_me;
      Parallel::send_receive(procup, pushed_ids[procup],
                             procdown, pushed_ids_to_me);
      Parallel::send_receive(procup, pushed_keys,
                             procdown, pushed_keys_to_me);
      Parallel::send_receive(procup, pushed_vals,
                             procdown, pushed_vals_to_me);
      libmesh_assert (pushed_ids_to_me.size() == pushed_keys_to_me.size());
      libmesh_assert (pushed_ids_to_me.size() == pushed_vals_to_me.size());

      // Add the constraints that I've been sent
      for (unsigned int i = 0; i != pushed_ids_to_me.size(); ++i)
        {
          libmesh_assert (pushed_keys_to_me[i].size() == pushed_vals_to_me[i].size());

          unsigned int constrained = pushed_ids_to_me[i];

          // If we don't already have a constraint for this dof,
          // add the one we were sent
          if (!this->is_constrained_dof(constrained))
            {
              DofConstraintRow &row = _dof_constraints[constrained];
              for (unsigned int j = 0; j != pushed_keys_to_me[i].size(); ++j)
                {
                  row[pushed_keys_to_me[i][j]] = pushed_vals_to_me[i][j];
                }
            }
        }
    }
  }

  // Now start checking for any other constraints we need
  // to know about, requesting them recursively.

  // Create a set containing the DOFs we already depend on
  typedef std::set<unsigned int> RCSet;
  RCSet unexpanded_set;

  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    {
      unexpanded_set.insert(i->first);
    }

  // We have to keep recursing while the unexpanded set is
  // nonempty on *any* processor
  unsigned int unexpanded_set_nonempty = !unexpanded_set.empty();
  Parallel::max(unexpanded_set_nonempty);

  while (unexpanded_set_nonempty)
    {
      // Request set
      RCSet request_set;

      // Request sets to send to each processor
      std::vector<std::vector<unsigned int> >
        requested_ids(libMesh::n_processors());

      // And the sizes of each
      std::vector<unsigned int> ids_on_proc(libMesh::n_processors(), 0);

      // Fill (and thereby sort and uniq!) the main request set
      for (RCSet::iterator i = unexpanded_set.begin();
           i != unexpanded_set.end(); ++i)
        {
          DofConstraintRow &row = _dof_constraints[*i];
          for (DofConstraintRow::iterator j = row.begin();
               j != row.end(); ++j)
            request_set.insert(j->first);
        }

      // Clear the unexpanded constraint set; we're about to expand it
      unexpanded_set.clear();

      // Count requests by processor
      unsigned int proc_id = 0;
      for (RCSet::iterator i = request_set.begin();
           i != request_set.end(); ++i)
        {
          while (*i >= _end_df[proc_id])
            proc_id++;
          ids_on_proc[proc_id]++;
        }
      for (unsigned int p = 0; p != libMesh::n_processors(); ++p)
        requested_ids[p].reserve(ids_on_proc[p]);

      // Prepare each processor's request set
      proc_id = 0;
      for (RCSet::iterator i = request_set.begin();
           i != request_set.end(); ++i)
        {
          while (*i >= _end_df[proc_id])
            proc_id++;
          requested_ids[proc_id].push_back(*i);
        }

      // Now request constraint rows from other processors
      for (unsigned int p=1; p != libMesh::n_processors(); ++p)
        {
          // Trade my requests with processor procup and procdown
          unsigned int procup = (libMesh::processor_id() + p) %
                                 libMesh::n_processors();
          unsigned int procdown = (libMesh::n_processors() +
                                   libMesh::processor_id() - p) %
                                   libMesh::n_processors();
          std::vector<unsigned int> request_to_fill;
          Parallel::send_receive(procup, requested_ids[procup],
                                 procdown, request_to_fill);

          // Fill those requests
          std::vector<std::vector<unsigned int> > row_keys(request_to_fill.size());
          std::vector<std::vector<Real> > row_vals(request_to_fill.size());
          for (unsigned int i=0; i != request_to_fill.size(); ++i)
            {
              unsigned int constrained = request_to_fill[i];
              if (_dof_constraints.count(constrained))
                {
                  DofConstraintRow &row = _dof_constraints[constrained];
                  unsigned int row_size = row.size();
                  row_keys[i].reserve(row_size);
                  row_vals[i].reserve(row_size);
                  for (DofConstraintRow::iterator j = row.begin();
                       j != row.end(); ++j)
                    {
                      row_keys[i].push_back(j->first);
                      row_vals[i].push_back(j->second);
                    }
                }
            }

          // Trade back the results
          std::vector<std::vector<unsigned int> > filled_keys;
          std::vector<std::vector<Real> > filled_vals;
          Parallel::send_receive(procdown, row_keys,
                                 procup, filled_keys);
          Parallel::send_receive(procdown, row_vals,
                                 procup, filled_vals);
          libmesh_assert (filled_keys.size() == requested_ids[procup].size());
          libmesh_assert (filled_vals.size() == requested_ids[procup].size());

          // Add any new constraint rows we've found
          for (unsigned int i=0; i != requested_ids[procup].size(); ++i)
            {
              libmesh_assert (filled_keys[i].size() == filled_vals[i].size());
              if (!filled_keys[i].empty())
                {
                  unsigned int constrained = requested_ids[procup][i];
                  DofConstraintRow &row = _dof_constraints[constrained];
                  for (unsigned int j = 0; j != filled_keys[i].size(); ++j)
                    row[filled_keys[i][j]] = filled_vals[i][j];

                  // And prepare to check for more recursive constraints
                  unexpanded_set.insert(constrained);
                }
            }
        }

      // We have to keep recursing while the unexpanded set is
      // nonempty on *any* processor
      unexpanded_set_nonempty = !unexpanded_set.empty();
      Parallel::max(unexpanded_set_nonempty);
    }
}



void DofMap::process_constraints ()
{
  // Create a set containing the DOFs we already depend on
  typedef std::set<unsigned int> RCSet;
  RCSet unexpanded_set;

  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    unexpanded_set.insert(i->first);

  while (!unexpanded_set.empty())
    for (RCSet::iterator i = unexpanded_set.begin();
	 i != unexpanded_set.end(); /* nothing */)
      {
	// If the DOF is constrained
	DofConstraints::iterator
	  pos = _dof_constraints.find(*i);
	
	libmesh_assert (pos != _dof_constraints.end());
	
	DofConstraintRow& constraint_row = pos->second;

	std::vector<unsigned int> constraints_to_expand;

	for (DofConstraintRow::const_iterator
	       it=constraint_row.begin(); it != constraint_row.end();
	     ++it)
	  if (this->is_constrained_dof(it->first))
            {
              unexpanded_set.insert(it->first);
	      constraints_to_expand.push_back(it->first);
	    }

	for (unsigned int j = 0; j != constraints_to_expand.size();
	     ++j)
	  {
            unsigned int expandable = constraints_to_expand[j];

	    DofConstraints::const_iterator
	      subpos = _dof_constraints.find(expandable);
	
	    libmesh_assert (subpos != _dof_constraints.end());
	
	    const DofConstraintRow& subconstraint_row = subpos->second;
            
	    for (DofConstraintRow::const_iterator
	           it=subconstraint_row.begin();
		   it != subconstraint_row.end(); ++it)
              {
		constraint_row[it->first] += it->second *
				constraint_row[expandable];
	      }
	    constraint_row.erase(expandable);
          }

	if (constraints_to_expand.empty())
	  unexpanded_set.erase(i++);
	else
	  i++;
      }

  // Now that we have our root constraint dependencies sorted out, add
  // them to the send_list
  this->add_constraints_to_send_list();
}



void DofMap::add_constraints_to_send_list()
{
  // This function must be run on all processors at once
  parallel_only();

  // Return immediately if there's nothing to gather
  if (libMesh::n_processors() == 1)
    return;

  // We might get to return immediately if none of the processors
  // found any constraints
  unsigned int has_constraints = !_dof_constraints.empty();
  Parallel::max(has_constraints);
  if (!has_constraints)
    return;

  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    {
      unsigned int constrained_dof = i->first;

      // We only need the dependencies of our own constrained dofs
      if (constrained_dof < this->first_dof() || 
          constrained_dof >= this->end_dof())
        continue;

      DofConstraintRow& constraint_row = i->second;
      for (DofConstraintRow::const_iterator
	   j=constraint_row.begin(); j != constraint_row.end();
	   ++j)
        {
          unsigned int constraint_dependency = j->first;

          // No point in adding one of our own dofs to the send_list
          if (constraint_dependency >= this->first_dof() && 
              constraint_dependency < this->end_dof())
            continue;

          _send_list.push_back(constraint_dependency);
        }
    }
}



#endif // LIBMESH_ENABLE_AMR || LIBMESH_ENABLE_PERIODIC


#ifdef LIBMESH_ENABLE_AMR

void DofMap::constrain_p_dofs (unsigned int var,
                               const Elem *elem,
                               unsigned int s,
                               unsigned int p)
{
  // We're constraining dofs on elem which correspond to p refinement
  // levels above p - this only makes sense if elem's p refinement
  // level is above p.
  libmesh_assert(elem->p_level() > p);
  libmesh_assert(s < elem->n_sides());

  const unsigned int sys_num = this->sys_number();
  const unsigned int dim = elem->dim();
  ElemType type = elem->type();
  FEType low_p_fe_type = this->variable_type(var);
  FEType high_p_fe_type = this->variable_type(var);
  low_p_fe_type.order = static_cast<Order>(low_p_fe_type.order + p);
  high_p_fe_type.order = static_cast<Order>(high_p_fe_type.order + 
                                            elem->p_level());

  const unsigned int n_nodes = elem->n_nodes();
  for (unsigned int n = 0; n != n_nodes; ++n)
    if (elem->is_node_on_side(n, s))
      {
        const Node * const node = elem->get_node(n);
        const unsigned int low_nc =
	  FEInterface::n_dofs_at_node (dim, low_p_fe_type, type, n);
        const unsigned int high_nc =
	  FEInterface::n_dofs_at_node (dim, high_p_fe_type, type, n);
	
	// since we may be running this method concurretly 
	// on multiple threads we need to acquire a lock 
	// before modifying the _dof_constraints object.
	Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
	
        if (elem->is_vertex(n))
          {
	    // Add "this is zero" constraint rows for high p vertex
            // dofs
            for (unsigned int i = low_nc; i != high_nc; ++i)
              _dof_constraints[node->dof_number(sys_num,var,i)].clear();
          }
        else
          {
            const unsigned int total_dofs = node->n_comp(sys_num, var);
            libmesh_assert(total_dofs >= high_nc);
	    // Add "this is zero" constraint rows for high p
            // non-vertex dofs, which are numbered in reverse
            for (unsigned int j = low_nc; j != high_nc; ++j)
              {
                const unsigned int i = total_dofs - j - 1;
                _dof_constraints[node->dof_number(sys_num,var,i)].clear();
              }
          }
      }
}

#endif // LIBMESH_ENABLE_AMR


#ifdef LIBMESH_ENABLE_PERIODIC

void DofMap::add_periodic_boundary (const PeriodicBoundary& periodic_boundary)
{
  untested();

  PeriodicBoundary boundary = periodic_boundary;
  PeriodicBoundary inverse_boundary;
  inverse_boundary.myboundary = boundary.pairedboundary;
  inverse_boundary.pairedboundary = boundary.myboundary;
  inverse_boundary.translation_vector = -boundary.translation_vector;

  std::pair<unsigned int, PeriodicBoundary> bp
    (boundary.myboundary, boundary);
  std::pair<unsigned int, PeriodicBoundary> ibp
    (boundary.pairedboundary, inverse_boundary);

  _periodic_boundaries.insert(bp);
  _periodic_boundaries.insert(ibp);
}


// ------------------------------------------------------------
// PeriodicBoundaries member functions

PeriodicBoundaries::~PeriodicBoundaries()
{
}

const Elem *PeriodicBoundaries::neighbor(unsigned int boundary_id,
					 const MeshBase &mesh,
                                         const Elem *e,
                                         unsigned int side)
{
  // Find a point on that side (and only that side)

  Point p = e->build_side(side)->centroid();

  PeriodicBoundary *b = this->boundary(boundary_id);
  libmesh_assert (b);
  p += b->translation_vector;

  return mesh.point_locator().operator()(p);
}

#endif
