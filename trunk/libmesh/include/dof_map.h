// $Id: dof_map.h,v 1.28 2003-05-21 13:50:19 benkirk Exp $

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



#ifndef __dof_map_h__
#define __dof_map_h__

// C++ Includes   -----------------------------------
#include <vector>
#include <map>
#include <string>

// Local Includes -----------------------------------
#include "mesh_common.h"
#include "enum_order.h"
#include "fe_type.h"
#include "coupling_matrix.h"
#include "reference_counted_object.h"

// Forward Declarations
class DofMap;
class Elem;
class MeshBase;
template <typename T> class DenseVector;
template <typename T> class DenseMatrix;
template <typename T> class SparseMatrix;




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
   * Attach the matrix that is used with this \p DofMap.
   */
  void attach_matrix (SparseMatrix<Number>& matrix);

  /**
   * Additional matrices may be handled with this \p DofMap.
   * They are initialized to the same sparsity structure as
   * the major matrix.
   */
  void attach_other_matrix (std::pair<std::string, SparseMatrix<Number>*> _new);


#ifdef ENABLE_AMR

  /**
   * A row of the Dof constraint matrix.  Store as a float
   * to save space since the factors are simple rational
   * fractions.
   */
  typedef std::map<unsigned int, float> DofConstraintRow;

  /** 
   * The constraint matrix storage format. 
   */
  typedef std::map<unsigned int, DofConstraintRow> DofConstraints;

#endif

  /**
   * Distrubute dofs on the current mesh.  Also builds the send list for
   * processor \p proc_id, which defaults to 0 for ease of use in serial
   * applications. 
   */
  void distribute_dofs(const MeshBase&);

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
  void add_variable (const FEType& type)
  { _variable_types.push_back (type); }
  
  /**
   * @returns the approximation order for variable \p c.
   */
  Order variable_order (const unsigned int c) const
  { return _variable_types[c].order; }
  
  /**
   * @returns the finite element type for variable \p c.
   */
  const FEType& variable_type (const unsigned int c) const
  { return _variable_types[c]; }
  
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
   * Returns the number of degrees of freedom on subdomain \p proc.
   */
  unsigned int n_dofs_on_processor(const unsigned int proc) const
  { assert(proc < _first_df.size()); return (_last_df[proc] - _first_df[proc]+1); }

  /**
   * Returns the first dof index that is local to subdomain \p proc.
   */
  unsigned int first_dof(unsigned int proc) const
  { assert(proc < _first_df.size()); return _first_df[proc]; }
  
  /**
   * Returns the last dof index that is local to subdomain \p proc.
   */
  unsigned int last_dof(unsigned int proc) const
  { assert(proc < _last_df.size()); return _last_df[proc]; }  


#ifdef ENABLE_AMR

  /**
   * @returns the total number of constrained degrees of freedom
   * in the problem.
   */
  unsigned int n_constrained_dofs() const { return _dof_constraints.size(); }

#endif

  
  /**
   * Fills the vector di with the global degree of freedom indices
   * for the element. If no variable number is specified then all
   * variables are returned.
   */
  void dof_indices (const Elem* elem,
		    std::vector<unsigned int>& di,
		    const unsigned int vn = static_cast<unsigned int>(-1)) const;




#ifdef ENABLE_AMR

  /**
   * Fills the vector di with the global degree of freedom indices
   * for the element using the \p DofMap::old_dof_object.
   * If no variable number is specified then all
   * variables are returned.
   */
  void old_dof_indices (const Elem* elem,
			std::vector<unsigned int>& di,
			const unsigned int vn = static_cast<unsigned int>(-1)) const;

  /**
   * Rebuilds the degree of freedom constraints.
   */ 
  void create_dof_constraints (const MeshBase& mesh);

  /**
   * Adds the user-defined row to the constraint matrix.
   * Issues a warning if the DOF was already constrained.
   */
  void add_constraint_row (const unsigned int dof_number,
			   const DofConstraintRow& constraint_row);
  
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
   */
  void constrain_element_matrix (DenseMatrix<Number>& matrix,
				 std::vector<unsigned int>& elem_dofs) const;
  
  /**
   * Constrains the element matrix.  This method allows the
   * element matrix to be non-square, in which case the row_dofs
   * and col_dofs may be of different size and correspond to
   * variables approximated in different spaces.
   */
  void constrain_element_matrix (DenseMatrix<Number>& matrix,
				 std::vector<unsigned int>& row_dofs,
				 std::vector<unsigned int>& col_dofs) const;
  
  /**
   * Constrains the element vector.
   */
  void constrain_element_vector (DenseVector<Number>&       rhs,
				 std::vector<unsigned int>& dofs) const;
  
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
					    std::vector<unsigned int>& elem_dofs) const;
  
#endif



  /**
   * Reinitialize the underlying data strucures conformal to the current mesh.
   */
  void reinit (const MeshBase& mesh);
  
  /**
   * Free all memory associated with the object, but keep the mesh pointer.
   */
  void clear ();
  


#ifdef ENABLE_AMR

  /**
   * Prints the \p _dof_constraints data structure.
   */
  void print_dof_constraints() const;

#endif
  
  /**
   * Degree of freedom coupling.  If left empty each DOF
   * couples to all others.  Can be used to reduce memory
   * requirements for sparse matrices.  DOF 0 might only
   * couple to itself, in which case \p dof_coupling(0,0)
   * should be 1 and \p dof_coupling(0,j) = 0 for j not equal
   * to 0.
   */
  CouplingMatrix _dof_coupling;

  
private:


  /**
   * @returns the number of the system we are responsible for.
   */
  unsigned int sys_number() const;
    
  
#ifdef ENABLE_AMR

  /**
   * Build the constraint matrix C associated with the element
   * degree of freedom indices elem_dofs.
   */
  void build_constraint_matrix (DenseMatrix<Number>& C,
				std::vector<unsigned int>& elem_dofs) const;

#endif


  /**
   * Finds all the DOFS associated with the element DOFs elem_dofs.
   * This will account for off-element couplings via hanging nodes.
   */
  void find_connected_dofs (std::vector<unsigned int>& elem_dofs) const;
  
  /**
   * The finite element type for each variable.
   */
  std::vector<FEType> _variable_types;

  /**
   * The number of the system we manage DOFs for.
   */
  const unsigned int _sys_number;
  
  /**
   * Pointer to the matrix used with this object.
   */
  SparseMatrix<Number>* _matrix;
  
  /**
   * Additional matrices handled by this object.  These pointers do @e 
   * not handle the memory, instead, \p SystemBase, who
   * told \p DofMap about them, owns them.
   */
  std::map<std::string, SparseMatrix<Number>* > _other_matrices;
  
  /**
   * Total number of degrees of freedom.
   */
  unsigned int _n_dfs;

  /**
   * First DOF index on processor \p p.
   */
  std::vector<unsigned int> _first_df;

  /**
   * Last DOF index (plus 1) on processor \p p.
   */
  std::vector<unsigned int> _last_df;

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


#ifdef ENABLE_AMR

  /**
   * Data structure containing DOF constraints.  The ith
   * entry is the constraint matrix row for DOF i.
   */
  DofConstraints _dof_constraints;

#endif
};



// ------------------------------------------------------------
// Dof Map inline member functions
inline
DofMap::DofMap(const unsigned int number) :
  _sys_number(number),
  _matrix(NULL),
  _n_dfs(0)  
{
  _other_matrices.clear();
}


inline
DofMap::~DofMap()
{
}



#ifdef ENABLE_AMR

inline
bool DofMap::is_constrained_dof (const unsigned int dof) const
{
  if (_dof_constraints.count(dof) != 0)
    return true;

  return false;
}

#endif


inline
unsigned int DofMap::sys_number() const
{
  return _sys_number;
}



#endif
