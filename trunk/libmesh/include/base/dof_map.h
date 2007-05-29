// $Id: dof_map.h,v 1.29 2007-05-29 23:30:09 roystgnr Exp $

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



#ifndef __dof_map_h__
#define __dof_map_h__

// C++ Includes   -----------------------------------
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>

// Local Includes -----------------------------------
#include "auto_ptr.h"
#include "libmesh_common.h"
#include "enum_order.h"
#include "reference_counted_object.h"
#include "libmesh.h" // libMesh::invalid_uint
#include "vector_value.h" // RealVectorValue

// Forward Declarations
class DofMap;
class Elem;
class MeshBase;
class Mesh;
class PointLocatorBase;
class FEType;
class CouplingMatrix;
class System;
template <typename T> class DenseVectorBase;
template <typename T> class DenseVector;
template <typename T> class DenseMatrix;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;




/**
 * This class handles the numbering of degrees of freedom on a mesh.
 * For systems of equations the class supports a fixed number of variables.
 * The degrees of freedom are numbered such that sequential, contiguous blocks
 * correspond to distinct subdomains.  This is so that the resulting data
 * structures will work well with parallel linear algebra packages.
 *
 * @author Benjamin S. Kirk, 2002-2003
 */

// ------------------------------------------------------------
// AMR constraint matrix types

#if defined(ENABLE_AMR) || defined(ENABLE_PERIODIC)
/**
 * A row of the Dof constraint matrix.
 */
typedef std::map<unsigned int, Real> DofConstraintRow;

/** 
 * The constraint matrix storage format. 
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.  Is there some issue with
 * deriving from standard containers, i.e. don't do it because they
 * don't have virtual destructors?
 */
class DofConstraints : public std::map<unsigned int, DofConstraintRow>
{
};
#endif // ENABLE_AMR || ENABLE_PERIODIC

  
// ------------------------------------------------------------
// Periodic boundary conditions information

#ifdef ENABLE_PERIODIC
/**
 * A row of the Dof constraint matrix.
 */
class PeriodicBoundary
{
public:
  // The boundary ID of this boundary and it's counterpart
  unsigned int myboundary,
	       pairedboundary;

  // One of these days we'll support rotated boundaries
  // RealTensor rotation_matrix;

  // The vector which is added to points in myboundary
  // to produce corresponding points in pairedboundary
  RealVectorValue translation_vector;
};


/** 
 * The constraint matrix storage format. 
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.  Is there some issue with
 * deriving from standard containers, i.e. don't do it because they
 * don't have virtual destructors?
 */
class PeriodicBoundaries : public std::map<unsigned int, PeriodicBoundary>
{
public:
  PeriodicBoundary *boundary(unsigned int id);

  // The periodic neighbor of \p e in direction \p side, if it
  // exists.  NULL otherwise
  const Elem *neighbor(unsigned int boundary_id, const Elem *e, unsigned int side);

  /**
   * Reinitialize the underlying data strucures conformal to the current mesh.
   */
  void reinit (MeshBase& mesh);

private:
  AutoPtr<PointLocatorBase> _point_locator;
};
#endif // ENABLE_PERIODIC

  
// ------------------------------------------------------------
// Dof Map class definition

class DofMap : public ReferenceCountedObject<DofMap>
{
public:

  /**
   * Constructor.  Requires the number of the system for which we
   * will be numbering degrees of freedom.
   */
  DofMap(const unsigned int sys_number);

  /**
   * Destructor.
   */
  ~DofMap();

  /**
   * Additional matrices may be handled with this \p DofMap.
   * They are initialized to the same sparsity structure as
   * the major matrix.
   */
  void attach_matrix (SparseMatrix<Number>& matrix);

  /**
   * Distrubute dofs on the current mesh.  Also builds the send list for
   * processor \p proc_id, which defaults to 0 for ease of use in serial
   * applications. 
   */
  void distribute_dofs (MeshBase&);

  /**
   * Computes the sparsity pattern for the matrix corresponding
   * to \p proc_id.  Produces data that can be fed to Petsc for
   * preallocation of sparse matrices.
   */
  void compute_sparsity (const MeshBase&);

  /**
   * Returns a constant reference to the \p _send_list for this processor.  The
   * \p _send_list contains the global indices of all the variables in the
   * global solution vector that influence the current processor.  This
   * information can be used for gathers at each solution step to retrieve
   * solution values needed for computation.
   */
  const std::vector<unsigned int>& get_send_list() const { return _send_list; }
  
  /**
   * Returns a constant reference to the \p _n_nz list for this processor.
   * The vector contains the bandwidth of the on-processor coupling for each
   * row of the global matrix that the current processor owns.  This
   * information can be used to preallocate space for a parallel sparse matrix.
   */
  const std::vector<unsigned int>& get_n_nz() const { return _n_nz; }
  
  /**
   * Returns a constant reference to the \p _n_oz list for this processor.
   * The vector contains the bandwidth of the off-processor coupling for each
   * row of the global matrix that the current processor owns.  This
   * information can be used to preallocate space for a parallel sparse matrix.
   */
  const std::vector<unsigned int>& get_n_oz() const { return _n_oz; }

  /**
   * Add an unknown of order \p order and finite element type
   * \p type to the system of equations.
   */
  void add_variable (const FEType& type);
  
  /**
   * @returns the approximation order for variable \p c.
   */
  Order variable_order (const unsigned int c) const;
  
  /**
   * @returns the finite element type for variable \p c.
   */
  const FEType& variable_type (const unsigned int c) const
  { return *_variable_types[c]; }
  
  /**
   * Returns the number of variables in the global solution vector. Defaults
   * to 1, should be 1 for a scalar equation, 3 for 2D incompressible Navier
   * Stokes (u,v,p), etc...
   */  
  unsigned int n_variables() const
  { return _variable_types.size(); }
  
  /**
   * @returns the total number of degrees of freedom in the problem.
   */
  unsigned int n_dofs() const { return _n_dfs; }

  /**
   * @returns the number of degrees of freedom on this processor.
   */
  unsigned int n_local_dofs () const
  { return this->n_dofs_on_processor (libMesh::processor_id()); }
    
  /**
   * Returns the number of degrees of freedom on subdomain \p proc.
   */
  unsigned int n_dofs_on_processor(const unsigned int proc) const
  { assert(proc < _first_df.size()); return (_end_df[proc] - _first_df[proc]); }

  /**
   * Returns the first dof index that is local to subdomain \p proc.
   */
  unsigned int first_dof(const unsigned int proc = libMesh::processor_id()) const
  { assert(proc < _first_df.size()); return _first_df[proc]; }
  
  /**
   * Returns the last dof index that is local to subdomain \p proc.
   * This function is now deprecated, because it returns nonsense in the rare
   * case where \p proc has no local dof indices.  Use end_dof() instead.
   */
  unsigned int last_dof(const unsigned int proc = libMesh::processor_id()) const
  { assert(proc < _end_df.size()); return (_end_df[proc] - 1); }  

  /**
   * Returns the first dof index that is after all indices local to subdomain \p proc.
   * Analogous to the end() member function of STL containers.
   */
  unsigned int end_dof(const unsigned int proc = libMesh::processor_id()) const
  { assert(proc < _end_df.size()); return _end_df[proc]; }  
  
  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the element. If no variable number is specified then all
   * variables are returned.
   */
  void dof_indices (const Elem* const elem,
		    std::vector<unsigned int>& di,
		    const unsigned int vn = libMesh::invalid_uint) const;

  /**
   * Builds the local element vector \p Ue from the global vector \p Ug,
   * accounting for any constrained degrees of freedom.  For an element
   * without constrained degrees of freedom this is the trivial mapping
   * \f$ Ue[i] = Ug[dof_indices[i]] \f$
   *
   * Note that the user must ensure that the element vector \p Ue is
   * properly sized when calling this method.  This is because there
   * is no \p resize() method in the \p DenseVectorBase<> class. 
   */
  void extract_local_vector (const NumericVector<Number>& Ug,
			     const std::vector<unsigned int>& dof_indices,
			     DenseVectorBase<Number>& Ue) const;


#if defined(ENABLE_AMR) || defined(ENABLE_PERIODIC)

  //--------------------------------------------------------------------
  // Constraint-specific methods
  /**
   * @returns the total number of constrained degrees of freedom
   * in the problem.
   */
  unsigned int n_constrained_dofs() const { return _dof_constraints.size(); }

  /**
   * Rebuilds the raw degree of freedom constraints.
   */ 
  void create_dof_constraints (const MeshBase& mesh);

  /**
   * Postprocesses any constrained degrees of freedom in elem_dofs
   * to be constrained only in terms of unconstrained dofs.
   */
  void process_recursive_constraints ();

  /**
   * Adds a copy of the user-defined row to the constraint matrix.
   * By default, produces an error if the DOF was already constrained.
   */
  void add_constraint_row (const unsigned int dof_number,
			   const DofConstraintRow& constraint_row,
			   const bool forbid_constraint_overwrite = true);
  
  /**
   * @returns true if the degree of freedom dof is constrained,
   * false otherwise.
   */
  bool is_constrained_dof (const unsigned int dof) const;
  
  /**
   * Constrains the element matrix.  This method requires the
   * element matrix to be square, in which case the elem_dofs
   * correspond to the global DOF indices of both the rows and
   * columns of the element matrix.  For this case the rows
   * and columns of the matrix necessarily correspond to variables
   * of the same approximation order.
   *
   * If \p asymmetric_constraint_rows is set to true (as it is by
   * default), constraint row equations will be reinforced in a way
   * which breaks matrix symmetry but makes inexact linear solver
   * solutions more likely to satisfy hanging node constraints.
   */
  void constrain_element_matrix (DenseMatrix<Number>& matrix,
				 std::vector<unsigned int>& elem_dofs,
				 bool asymmetric_constraint_rows = true) const;
  
  /**
   * Constrains the element matrix.  This method allows the
   * element matrix to be non-square, in which case the row_dofs
   * and col_dofs may be of different size and correspond to
   * variables approximated in different spaces.
   */
  void constrain_element_matrix (DenseMatrix<Number>& matrix,
				 std::vector<unsigned int>& row_dofs,
				 std::vector<unsigned int>& col_dofs,
				 bool asymmetric_constraint_rows = true) const;
  
  /**
   * Constrains the element vector.
   */
  void constrain_element_vector (DenseVector<Number>&       rhs,
				 std::vector<unsigned int>& dofs,
				 bool asymmetric_constraint_rows = true) const;
  
  /**
   * Constrains the element matrix and vector.  This method requires
   * the element matrix to be square, in which case the elem_dofs
   * correspond to the global DOF indices of both the rows and
   * columns of the element matrix.  For this case the rows
   * and columns of the matrix necessarily correspond to variables
   * of the same approximation order.
   */
  void constrain_element_matrix_and_vector (DenseMatrix<Number>& matrix,
					    DenseVector<Number>& rhs,
					    std::vector<unsigned int>& elem_dofs,
					    bool asymmetric_constraint_rows = true) const;
  
  /**
   * Constrains the numeric vector \p v, which represents a solution defined on
   * the mesh.  This may need to be used after a linear solve, if your linear
   * solver's solutions do not satisfy your DoF constraints to a tight enough
   * tolerance.
   *
   * If \p v == NULL, the system solution vector is constrained
   */
  void enforce_constraints_exactly (const System &system,
				    NumericVector<Number> *v = NULL) const;
  
  /**
   * Prints the \p _dof_constraints data structure.
   */
  void print_dof_constraints() const;

  /**
   * Tests the constrained degrees of freedom on the numeric vector \p v, which
   * represents a solution defined on the mesh, returning a pair whose first
   * entry is the maximum absolute error on a constrained DoF and whose second
   * entry is the maximum relative error.  Useful for debugging purposes.
   *
   * If \p v == NULL, the system solution vector is tested.
   */
  std::pair<Real, Real> max_constraint_error(const System &system,
				    NumericVector<Number> *v = NULL) const;

#endif // ENABLE_AMR || ENABLE_PERIODIC


#ifdef ENABLE_PERIODIC

  //--------------------------------------------------------------------
  // PeriodicBoundary-specific methods

  /**
   * Adds a copy of the specified periodic boundary to the system.
   */
  void add_periodic_boundary (const PeriodicBoundary& periodic_boundary);

  /**
   * @returns true if the boundary given by \p boundaryid is periodic,
   * false otherwise
   */
  bool is_periodic_boundary (const unsigned int boundaryid) const;

#endif // ENABLE_PERIODIC


#ifdef ENABLE_AMR

  //--------------------------------------------------------------------
  // AMR-specific methods

  /**
   * After a mesh is refined and repartitioned it is possible that the
   * \p _send_list will need to be augmented.  This is the case when an
   * element is refined and its children end up on different processors
   * than the parent.  These children will need values from the parent
   * when projecting the solution onto the refined mesh, hence the parent's
   * DOF indices need to be included in the \p _send_list.
   */
  void augment_send_list_for_projection(const MeshBase &);

  /**
   * Fills the vector di with the global degree of freedom indices
   * for the element using the \p DofMap::old_dof_object.
   * If no variable number is specified then all
   * variables are returned.
   */
  void old_dof_indices (const Elem* const elem,
			std::vector<unsigned int>& di,
			const unsigned int vn = libMesh::invalid_uint) const;
  /**
   * @returns the total number of degrees of freedom on old_dof_objects
   * 
   */
  unsigned int n_old_dofs() const { return _n_old_dfs; }

  /**
   * Constrains degrees of freedom on side \p s of element \p elem which
   * correspond to variable number \p var and to p refinement levels 
   * above \p p.
   */
  void constrain_p_dofs (unsigned int var,
			 const Elem *elem,
			 unsigned int s,
			 unsigned int p);
  
#endif // ENABLE_AMR

  /**
   * Reinitialize the underlying data strucures conformal to the current mesh.
   */
  void reinit (MeshBase& mesh);
  
  /**
   * Free all memory associated with the object, but keep the mesh pointer.
   */
  void clear ();
  
  /**
   * Degree of freedom coupling.  If left empty each DOF
   * couples to all others.  Can be used to reduce memory
   * requirements for sparse matrices.  DOF 0 might only
   * couple to itself, in which case \p dof_coupling(0,0)
   * should be 1 and \p dof_coupling(0,j) = 0 for j not equal
   * to 0.
   *
   * This variable is named as though it were class private,
   * but it is in the public interface.  Also there are no
   * public methods for accessing it...  This typically means
   * you should only use it if you know what you are doing.
   */
  CouplingMatrix* _dof_coupling;

  
private:
  
  /**
   * @returns the number of the system we are responsible for.
   */
  unsigned int sys_number() const;
  
  /**
   * Distributes the global degrees of freedom in a 
   * variable-major ordering.  In this format the local
   * degrees of freedom are in a contiguous block for each
   * variable in the system.
   */  
  void distribute_dofs_var_major (MeshBase&);
  
  /**
   * Distributes the global degrees of freedom in a 
   * node/element-major ordering.  In this format all the 
   * degrees of freedom at a node/element are in contiguous 
   * blocks.  Note in particular that the degrees of freedom
   * for a given variable are not in contiguous blocks, as 
   * in the case of \p distribute_dofs_var_major.
   */  
  void distribute_dofs_node_major (MeshBase&);

  /**
   * Splices the two sorted ranges [begin,middle) and [middle,end)
   * into one sorted range [begin,end).  This method is much like
   * std::inplace_merge except it assumes the intersection
   * of the two sorted ranges is empty and that any element in
   * each range occurs only once in that range.  Additionally,
   * this sort occurs in-place, while std::inplace_merge may
   * use a temporary buffer.
   */ 
  template<typename BidirectionalIterator>
  void sort_sparsity_row (const BidirectionalIterator begin,
			  BidirectionalIterator       middle,
			  const BidirectionalIterator end);
    
  /**
   * Adds entries to the \p _send_list vector corresponding to DoFs
   * on elements neighboring the current processor.
   */
  void add_neighbors_to_send_list(MeshBase& mesh);

  /**
   * Takes the \p _send_list vector (which may have duplicate entries)
   * and sorts it.  The duplicate entries are then removed, resulting in
   * a sorted \p _send_list with unique entries.
   */
  void sort_send_list ();

  
#if defined(ENABLE_AMR) || defined(ENABLE_PERIODIC)

  /**
   * Build the constraint matrix C associated with the element
   * degree of freedom indices elem_dofs. The optional parameter
   * \p called_recursively should be left at the default value
   * \p false.  This is used to handle the special case of
   * an element's degrees of freedom being constrained in terms
   * of other, local degrees of freedom.  The usual case is
   * for an elements DOFs to be constrained by some other,
   * external DOFs.   
   */
  void build_constraint_matrix (DenseMatrix<Number>& C,
				std::vector<unsigned int>& elem_dofs,
				const bool called_recursively=false) const;

  /**
   * Finds all the DOFS associated with the element DOFs elem_dofs.
   * This will account for off-element couplings via hanging nodes.
   */
  void find_connected_dofs (std::vector<unsigned int>& elem_dofs) const;
  
#endif


  /**
   * The finite element type for each variable.
   */
  std::vector<FEType*> _variable_types;

  /**
   * The number of the system we manage DOFs for.
   */
  const unsigned int _sys_number;
  
  /**
   * Additional matrices handled by this object.  These pointers do @e 
   * not handle the memory, instead, \p System, who
   * told \p DofMap about them, owns them.
   */
  std::vector<SparseMatrix<Number>* > _matrices;
  
  /**
   * First DOF index on processor \p p.
   */
  std::vector<unsigned int> _first_df;

  /**
   * Last DOF index (plus 1) on processor \p p.
   */
  std::vector<unsigned int> _end_df;

  /**
   * A list containing all the global DOF indicies that affect the
   * solution on my subdomain.
   */
  std::vector<unsigned int> _send_list;

  /**
   * The number of on-processor nonzeros in my portion of the
   * global matrix.
   */
  std::vector<unsigned int> _n_nz;
  
  /**
   * The number of off-processor nonzeros in my portion of the
   * global matrix.
   */
  std::vector<unsigned int> _n_oz;

  /**
   * Total number of degrees of freedom.
   */
  unsigned int _n_dfs;

#ifdef ENABLE_AMR

  /**
   * Total number of degrees of freedom on old dof objects
   */
  unsigned int _n_old_dfs;

  /**
   * Data structure containing DOF constraints.  The ith
   * entry is the constraint matrix row for DOF i.
   */
  DofConstraints _dof_constraints;

#endif


#ifdef ENABLE_PERIODIC
  /**
   * Data structure containing periodic boundaries.  The ith
   * entry is the constraint matrix row for boundaryid i.
   */
  PeriodicBoundaries _periodic_boundaries;
#endif


#if defined(__GNUC__) && (__GNUC__ < 4) && !defined(__INTEL_COMPILER)
  /**
   * Dummy function that does nothing but can be used to prohibit
   * compiler optimization in some situations where some compilers
   * have optimization bugs.
   */
  static void _dummy_function(void);
#endif
};



// ------------------------------------------------------------
// Dof Map inline member functions
inline
DofMap::DofMap(const unsigned int number) :
  _dof_coupling(NULL),
  _sys_number(number),
//  _matrix(NULL),
  _n_dfs(0) 
#ifdef ENABLE_AMR
  , _n_old_dfs(0)
#endif
{
  _matrices.clear();
}


#if defined(ENABLE_AMR) || defined(ENABLE_PERIODIC)

inline
bool DofMap::is_constrained_dof (const unsigned int dof) const
{
  if (_dof_constraints.count(dof) != 0)
    return true;

  return false;
}

#endif


#ifdef ENABLE_PERIODIC

inline
bool DofMap::is_periodic_boundary (const unsigned int boundaryid) const
{
  if (_periodic_boundaries.count(boundaryid) != 0)
    return true;

  return false;
}

#endif


inline
unsigned int DofMap::sys_number() const
{
  return _sys_number;
}



template<typename BidirectionalIterator>
inline
void DofMap::sort_sparsity_row (const BidirectionalIterator begin,
				BidirectionalIterator       middle,
				const BidirectionalIterator end)
{
  if ((begin == middle) || (middle == end)) return;

  assert (std::distance (begin,  middle) > 0);
  assert (std::distance (middle, end)    > 0);
  assert (std::unique (begin,  middle) == middle);
  assert (std::unique (middle, end)    == end);
  
  while (middle != end)
    {
      BidirectionalIterator
	b = middle,
	a = b-1;

      // Bubble-sort the middle value downward
      while (!(*a < *b)) // *a & *b are less-than comparable, so use <
	{
	  std::swap (*a, *b);

#if defined(__GNUC__) && (__GNUC__ < 4) && !defined(__INTEL_COMPILER)
	  /* Prohibit optimization at this point since gcc 3.3.5 seems
	     to have a bug.  */
	  this->_dummy_function();
#endif

	  if (a == begin) break;
	  
	  b=a;
	  --a;
	}
      
      ++middle;
    }

  // Assure the algorithm worked if we are in DEBUG mode
#ifdef DEBUG
  {
    // SGI STL extension!
    // assert (std::is_sorted(begin,end));
    
    BidirectionalIterator
      prev  = begin,
      first = begin;

    for (++first; first != end; prev=first, ++first)
      if (*first < *prev)
	assert(false);
  }  
#endif

  // Make sure the two ranges did not contain any common elements
  assert (std::unique (begin, end) == end);
}


// ------------------------------------------------------------
// PeriodicBoundary inline member functions

#ifdef ENABLE_PERIODIC

inline
PeriodicBoundary *PeriodicBoundaries::boundary(unsigned int id)
{
  iterator i = this->find(id);
  if (i == this->end())
    return NULL;
  return &i->second;
}

#endif // ENABLE_PERIODIC


#endif // __dof_map_h__
