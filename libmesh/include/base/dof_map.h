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



#ifndef __dof_map_h__
#define __dof_map_h__

// C++ Includes   -----------------------------------
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>

// Local Includes -----------------------------------
#include "libmesh_common.h"
#include "enum_order.h"
#include "reference_counted_object.h"
#include "libmesh.h" // libMesh::invalid_uint
#include "vector_value.h" // RealVectorValue
#include "system.h" // System::Variable
#include "threads.h"
#include "threads_allocators.h"
#include "elem_range.h"


// Forward Declarations
class DofMap;
class DofObject;
class Elem;
class MeshBase;
class Mesh;
class FEType;
class CouplingMatrix;
template <typename T> class DenseVectorBase;
template <typename T> class DenseVector;
template <typename T> class DenseMatrix;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;




// ------------------------------------------------------------
// Sparsity Pattern

/**
 * This defines the sparsity pattern, or graph, of a sparse matrix.
 * The format is quite simple -- the global indices of the nonzero entries
 * in each row are packed into a vector.  The global indices (i,j) of the 
 * nth nonzero entry of row i are given by j = sparsity_pattern[i][n];
 */
namespace SparsityPattern // use a namespace so member classes can be forward-declared.
{
  typedef std::vector<unsigned int, Threads::scalable_allocator<unsigned int> > Row;
  class Graph : public std::vector<Row> {};
  
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
  static void sort_row (const BidirectionalIterator begin,
			BidirectionalIterator       middle,
			const BidirectionalIterator end);

  /**
   * This helper class can be called on multiple threads to compute 
   * the sparsity pattern (or graph) of the sparse matrix resulting
   * from the discretization.  This pattern may be used directly by
   * a particular sparse matrix format (e.g. \p LaspackMatrix)
   * or indirectly (e.g. \p PetscMatrix).  In the latter case the
   * number of nonzeros per row of the matrix is needed for efficient 
   * preallocation.  In this case it suffices to provide estimate
   * (but bounding) values, and in this case the threaded method can
   * take some short-cuts for efficiency.
   */
  class Build
  {
  private:
    const MeshBase &mesh;
    const DofMap &dof_map;
    const CouplingMatrix *dof_coupling;
    const bool implicit_neighbor_dofs;
    const bool need_full_sparsity_pattern;

  public:

    SparsityPattern::Graph sparsity_pattern;
    std::vector<unsigned int> n_nz;
    std::vector<unsigned int> n_oz;
    
    Build (const MeshBase &mesh_in,
	   const DofMap &dof_map_in,
	   const CouplingMatrix *dof_coupling_in,
	   const bool implicit_neighbor_dofs_in,
	   const bool need_full_sparsity_pattern_in) :
      mesh(mesh_in),
      dof_map(dof_map_in),
      dof_coupling(dof_coupling_in),
      implicit_neighbor_dofs(implicit_neighbor_dofs_in),
      need_full_sparsity_pattern(need_full_sparsity_pattern_in)
    {}

    Build (Build &other, Threads::split) :
      mesh(other.mesh),
      dof_map(other.dof_map),
      dof_coupling(other.dof_coupling),
      implicit_neighbor_dofs(other.implicit_neighbor_dofs),
      need_full_sparsity_pattern(other.need_full_sparsity_pattern)
    {}

    void operator()(const ConstElemRange &range);

    void join (const Build &other);
  };

#if defined(__GNUC__) && (__GNUC__ < 4) && !defined(__INTEL_COMPILER)
  /**
   * Dummy function that does nothing but can be used to prohibit
   * compiler optimization in some situations where some compilers
   * have optimization bugs.
   */
  void _dummy_function(void);
#endif
  
}




// ------------------------------------------------------------
// AMR constraint matrix types

#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)
/**
 * A row of the Dof constraint matrix.
 */
typedef std::map<unsigned int, Real, 
                 std::less<unsigned int>, 
                 Threads::scalable_allocator<std::pair<const unsigned int, Real> > > DofConstraintRow;

/** 
 * The constraint matrix storage format. 
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.  Don't delete this from
 * a pointer-to-std::map; the destructor isn't virtual!
 */
class DofConstraints : public std::map<unsigned int, 
                                       DofConstraintRow, 
                                       std::less<unsigned int>, 
                                       Threads::scalable_allocator<std::pair<const unsigned int, DofConstraintRow> > >
{
};
#endif // LIBMESH_ENABLE_AMR || LIBMESH_ENABLE_PERIODIC

  
// ------------------------------------------------------------
// Periodic boundary conditions information

#ifdef LIBMESH_ENABLE_PERIODIC
/**
 * The definition of a periodic boundary.
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

  PeriodicBoundaries() {}

  ~PeriodicBoundaries();

  // The periodic neighbor of \p e in direction \p side, if it
  // exists.  NULL otherwise
  const Elem *neighbor(unsigned int boundary_id, const MeshBase &mesh, const Elem *e, unsigned int side);

private:
};
#endif // LIBMESH_ENABLE_PERIODIC


  
// ------------------------------------------------------------
// DofMap class definition

/**
 * This class handles the numbering of degrees of freedom on a mesh.
 * For systems of equations the class supports a fixed number of variables.
 * The degrees of freedom are numbered such that sequential, contiguous blocks
 * correspond to distinct subdomains.  This is so that the resulting data
 * structures will work well with parallel linear algebra packages.
 *
 * @author Benjamin S. Kirk, 2002-2007
 */
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
   * Takes the \p _send_list vector (which may have duplicate entries)
   * and sorts it.  The duplicate entries are then removed, resulting in
   * a sorted \p _send_list with unique entries.
   */
  void prepare_send_list ();
  
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
  void add_variable (const System::Variable &var);
  
  /**
   * @returns the variable description object for variable \p c.
   */
  const System::Variable& variable (const unsigned int c) const;

  /**
   * @returns the approximation order for variable \p c.
   */
  Order variable_order (const unsigned int c) const;
  
  /**
   * @returns the finite element type for variable \p c.
   */
  const FEType& variable_type (const unsigned int c) const;
  
  /**
   * Returns the number of variables in the global solution vector. Defaults
   * to 1, should be 1 for a scalar equation, 3 for 2D incompressible Navier
   * Stokes (u,v,p), etc...
   */  
  unsigned int n_variables() const
  { return _variables.size(); }
  
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
  { libmesh_assert(proc < _first_df.size()); return (_end_df[proc] - _first_df[proc]); }

  /**
   * Returns the first dof index that is local to subdomain \p proc.
   */
  unsigned int first_dof(const unsigned int proc = libMesh::processor_id()) const
  { libmesh_assert(proc < _first_df.size()); return _first_df[proc]; }
  
  /**
   * Returns the last dof index that is local to subdomain \p proc.
   * This function is now deprecated, because it returns nonsense in the rare
   * case where \p proc has no local dof indices.  Use end_dof() instead.
   */
  unsigned int last_dof(const unsigned int proc = libMesh::processor_id()) const
  { libmesh_assert(proc < _end_df.size()); return (_end_df[proc] - 1); }  

  /**
   * Returns the first dof index that is after all indices local to subdomain \p proc.
   * Analogous to the end() member function of STL containers.
   */
  unsigned int end_dof(const unsigned int proc = libMesh::processor_id()) const
  { libmesh_assert(proc < _end_df.size()); return _end_df[proc]; }  
  
  /**
   * Returns the first local degree of freedom index for variable \p var.
   */
  unsigned int variable_first_local_dof (const unsigned int var) const
  { 
    libmesh_assert ((var+1) < _var_first_local_df.size());
    libmesh_assert (_var_first_local_df[var] != DofObject::invalid_id);
    return _var_first_local_df[var];
  }

  /**
   * Returns (one past) the last local degree of freedom index for variable \p var.
   * Analogous to the end() member function of STL containers.
   */
  unsigned int variable_last_local_dof (const unsigned int var) const
  { 
    libmesh_assert ((var+1) < _var_first_local_df.size());
    libmesh_assert (_var_first_local_df[var+1] != DofObject::invalid_id);
    return _var_first_local_df[var+1];
  }

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the element. If no variable number is specified then all
   * variables are returned.
   */
  void dof_indices (const Elem* const elem,
		    std::vector<unsigned int>& di,
		    const unsigned int vn = libMesh::invalid_uint) const;

  /**
   * Tells other library functions whether or not this problem
   * includes coupling between dofs in neighboring cells, as can
   * currently be specified on the command line or inferred from
   * the use of all discontinuous variables.
   */
  bool use_coupled_neighbor_dofs(const MeshBase& mesh) const;

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


#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)

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
   * Gathers any relevant constraint equations from other processors
   */
  void allgather_recursive_constraints ();

  /**
   * Postprocesses any constrained degrees of freedom in elem_dofs
   * to be constrained only in terms of unconstrained dofs, then adds
   * unconstrained dofs to the send_list and prepares that for use.
   * This should be run after both system (create_dof_constraints) and
   * user constraints have all been added.
   */
  void process_constraints ();

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
   * Prints the \p _dof_constraints data structure.
   */
  void print_dof_constraints(std::ostream& os=std::cout) const;

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

#endif // LIBMESH_ENABLE_AMR || LIBMESH_ENABLE_PERIODIC

  //--------------------------------------------------------------------
  // Constraint-specific methods
  // Some of these methods are enabled (but inlined away to nothing)
  // when constraints are disabled at configure-time.  This is to
  // increase API compatibility of user code with different library
  // builds.

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
   * Constrains a dyadic element matrix B = v w'.  This method
   * requires the element matrix to be square, in which case the
   * elem_dofs correspond to the global DOF indices of both the rows
   * and columns of the element matrix.  For this case the rows and
   * columns of the matrix necessarily correspond to variables of the
   * same approximation order.
   */
  void constrain_element_dyad_matrix (DenseVector<Number>& v,
				      DenseVector<Number>& w,
				      std::vector<unsigned int>& row_dofs,
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
  

#ifdef LIBMESH_ENABLE_PERIODIC

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

#endif // LIBMESH_ENABLE_PERIODIC


#ifdef LIBMESH_ENABLE_AMR

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
  // void augment_send_list_for_projection(const MeshBase &);

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
  
#endif // LIBMESH_ENABLE_AMR

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
   * @invalidates all active DofObject dofs for this system
   */
  void invalidate_dofs(MeshBase& mesh) const;

  /**
   * An adapter function that returns Node pointers by index
   */
  DofObject* node_ptr(MeshBase& mesh, unsigned int i) const;

  /**
   * An adapter function that returns Elem pointers by index
   */
  DofObject* elem_ptr(MeshBase& mesh, unsigned int i) const;

  /**
   * A member function type like node_ptr or elem_ptr
   */
  typedef DofObject* (DofMap::*dofobject_accessor)
    (MeshBase& mesh, unsigned int i) const;
  
  /**
   * Helper function for distributing dofs in parallel
   */
  template<typename iterator_type>
  void set_nonlocal_dof_objects(iterator_type objects_begin,
			        iterator_type objects_end,
			        MeshBase &mesh,
			        dofobject_accessor objects);

  /**
   * Distributes the global degrees of freedom, for dofs on
   * this processor.  In this format the local
   * degrees of freedom are in a contiguous block for each
   * variable in the system.
   * Starts at index next_free_dof, and increments it to
   * the post-final index.
   */  
  void distribute_local_dofs_var_major (unsigned int& next_free_dof,
				        MeshBase& mesh);
  
  /**
   * Distributes the global degrees of freedom, for dofs on
   * this processor.  In this format all the 
   * degrees of freedom at a node/element are in contiguous 
   * blocks.  Note in particular that the degrees of freedom
   * for a given variable are not in contiguous blocks, as 
   * in the case of \p distribute_local_dofs_var_major.
   * Starts at index next_free_dof, and increments it to
   * the post-final index.
   * If build_send_list is true, builds the send list.  If
   * false, clears and reserves the send list
   */  
  void distribute_local_dofs_node_major (unsigned int& next_free_dof,
				         MeshBase& mesh);
    
  /**
   * Adds entries to the \p _send_list vector corresponding to DoFs
   * on elements neighboring the current processor.
   */
  void add_neighbors_to_send_list(MeshBase& mesh);

#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)

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
  void find_connected_dofs (std::vector<unsigned int> &elem_dofs) const;
  
  /**
   * Adds entries to the \p _send_list vector corresponding to DoFs
   * which are dependencies for constraint equations on the current
   * processor.
   */
  void add_constraints_to_send_list();

#endif


  /**
   * The finite element type for each variable.
   */
  std::vector<System::Variable> _variables;

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
   * The first local DOF index for each variable in the \p System.
   */
  std::vector<unsigned int> _var_first_local_df;
  
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

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Total number of degrees of freedom on old dof objects
   */
  unsigned int _n_old_dfs;

#endif

#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)

  /**
   * Data structure containing DOF constraints.  The ith
   * entry is the constraint matrix row for DOF i.
   */
  DofConstraints _dof_constraints;

#endif


#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * Data structure containing periodic boundaries.  The ith
   * entry is the constraint matrix row for boundaryid i.
   */
  PeriodicBoundaries _periodic_boundaries;
#endif

  friend class SparsityPattern::Build;
};



// ------------------------------------------------------------
// Dof Map inline member functions
inline
DofMap::DofMap(const unsigned int number) :
  _dof_coupling(NULL),
  _sys_number(number),
//  _matrix(NULL),
  _n_dfs(0) 
#ifdef LIBMESH_ENABLE_AMR
  , _n_old_dfs(0)
#endif
{
  _matrices.clear();
}


#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)

inline
bool DofMap::is_constrained_dof (const unsigned int dof) const
{
  if (_dof_constraints.count(dof) != 0)
    return true;

  return false;
}

#endif


#ifdef LIBMESH_ENABLE_PERIODIC

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



// ------------------------------------------------------------
// SparsityPattern inline member functions
template<typename BidirectionalIterator>
inline
void SparsityPattern::sort_row (const BidirectionalIterator begin,
				BidirectionalIterator       middle,
				const BidirectionalIterator end)
{
  if ((begin == middle) || (middle == end)) return;

  libmesh_assert (std::distance (begin,  middle) > 0);
  libmesh_assert (std::distance (middle, end)    > 0);
  libmesh_assert (std::unique (begin,  middle) == middle);
  libmesh_assert (std::unique (middle, end)    == end);
  
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
	  SparsityPattern::_dummy_function();
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
    // libmesh_assert (std::is_sorted(begin,end));
    
    BidirectionalIterator
      prev  = begin,
      first = begin;

    for (++first; first != end; prev=first, ++first)
      if (*first < *prev)
	libmesh_assert(false);
  }  
#endif

  // Make sure the two ranges did not contain any common elements
  libmesh_assert (std::unique (begin, end) == end);
}


// ------------------------------------------------------------
// PeriodicBoundary inline member functions
#ifdef LIBMESH_ENABLE_PERIODIC

inline
PeriodicBoundary *PeriodicBoundaries::boundary(unsigned int id)
{
  iterator i = this->find(id);
  if (i == this->end())
    return NULL;
  return &i->second;
}

#endif // LIBMESH_ENABLE_PERIODIC

#if !defined(LIBMESH_ENABLE_AMR) && !defined(LIBMESH_ENABLE_PERIODIC)
  //--------------------------------------------------------------------
  // Constraint-specific methods get inlined into nothing if
  // constraints are disabled, so there's no reason for users not to 
  // use them.

inline void DofMap::constrain_element_matrix (DenseMatrix<Number>&,
				              std::vector<unsigned int>&,
				              bool) const {}
  
inline void DofMap::constrain_element_matrix (DenseMatrix<Number>&,
				              std::vector<unsigned int>&,
				              std::vector<unsigned int>&,
				              bool) const {}
  
inline void DofMap::constrain_element_vector (DenseVector<Number>&,
				              std::vector<unsigned int>&,
				              bool) const {}
  
inline void DofMap::constrain_element_matrix_and_vector (DenseMatrix<Number>&,
					                 DenseVector<Number>&,
					                 std::vector<unsigned int>&,
					                 bool) const {}

inline void DofMap::constrain_element_dyad_matrix (DenseVector<Number>&,
				                   DenseVector<Number>&,
				                   std::vector<unsigned int>&,
				                   bool) const {}

inline void DofMap::enforce_constraints_exactly (const System &,
				                 NumericVector<Number> *) const {}

#endif // !defined(LIBMESH_ENABLE_AMR) && !defined(LIBMESH_ENABLE_PERIODIC)

#endif // __dof_map_h__
