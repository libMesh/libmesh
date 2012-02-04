// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "variable.h"
#include "threads.h"
#include "threads_allocators.h"
#include "elem_range.h"
#include "sparsity_pattern.h"

namespace libMesh
{

// Forward Declarations
class CouplingMatrix;
class DofMap;
class DofObject;
class Elem;
class FEType;
class MeshBase;
class Mesh;
class PeriodicBoundary;
class PeriodicBoundaries;
namespace SparsityPattern { class Build; }
class System;
template <typename T> class DenseVectorBase;
template <typename T> class DenseVector;
template <typename T> class DenseMatrix;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;



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

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
/**
 * A row of the Node constraint mapping.  Currently this just
 * stores the topology of the constrained Nodes, but for forward
 * compatibility we also include coefficients, so we could add
 * Lagrange-positioned-node constraints later.
 */
typedef std::map<const Node *, Real,
                 std::less<const Node *>,
                 Threads::scalable_allocator<std::pair<const Node * const, Real> > > NodeConstraintRow;

/**
 * The Node constraint storage format.
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.  Don't delete this from
 * a pointer-to-std::map; the destructor isn't virtual!
 */
class NodeConstraints : public std::map<const Node *,
                                       NodeConstraintRow,
                                       std::less<const Node *>,
                                       Threads::scalable_allocator<std::pair<const Node * const, NodeConstraintRow> > >
{
};
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#endif // LIBMESH_ENABLE_AMR || LIBMESH_ENABLE_PERIODIC



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
   * Abstract base class to be used to add user-defined implicit
   * degree of freedom couplings.
   */
  class AugmentSparsityPattern
  {
  public:
    virtual ~AugmentSparsityPattern () {};

    /**
     * User-defined function to augment the sparsity pattern.
     */
    virtual void augment_sparsity_pattern (SparsityPattern::Graph & sparsity,
					   std::vector<unsigned int> & n_nz,
					   std::vector<unsigned int> & n_oz) = 0;
  };

  /**
   * Abstract base class to be used to add user-defined parallel
   * degree of freedom couplings.
   */
  class AugmentSendList
  {
  public:
    virtual ~AugmentSendList () {};

    /**
     * User-defined function to augment the send list.
     */
    virtual void augment_send_list (std::vector<unsigned int> & send_list) = 0;
  };
   
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
   * Attach an object to use to populate the
   * sparsity pattern with extra entries.
   *
   * Care must be taken that when adding entries they are sorted into the Rows
   *
   * Further, you _must_ modify n_nz and n_oz properly!
   *
   * This is an advanced function... use at your own peril!
   */
  void attach_extra_sparsity_object (DofMap::AugmentSparsityPattern &asp)
  {
    _augment_sparsity_pattern = &asp;
  }

  /**
   * Attach a function pointer to use as a callback to populate the
   * sparsity pattern with extra entries.
   *
   * Care must be taken that when adding entries they are sorted into the Rows
   *
   * Further, you _must_ modify n_nz and n_oz properly!
   *
   * This is an advanced function... use at your own peril!
   */
  void attach_extra_sparsity_function(void (*func)(SparsityPattern::Graph & sparsity,
                                                   std::vector<unsigned int> & n_nz,
                                                   std::vector<unsigned int> & n_oz,
                                                   void *),
                                      void * context = NULL)
  { _extra_sparsity_function = func; _extra_sparsity_context = context; }

  /**
   * Attach an object to populate the send_list with extra entries.
   * This should only add to the send list, but no checking is done 
   * to enforce this behavior.
   *
   * This is an advanced function... use at your own peril!
   */
  void attach_extra_send_list_object (DofMap::AugmentSendList &asl)
  {
    _augment_send_list = &asl;
  }

  /**
   * Attach a function pointer to use as a callback to populate the
   * send_list with extra entries.
   */
  void attach_extra_send_list_function(void (*func)(std::vector<unsigned int> &, void *), void * context = NULL)
  { _extra_send_list_function = func; _extra_send_list_context = context; }
  
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
  void add_variable (const Variable &var);

  /**
   * @returns the variable description object for variable \p c.
   */
  const Variable& variable (const unsigned int c) const;

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
   * @returns the number of SCALAR dofs.
   */
  unsigned int n_SCALAR_dofs() const { return _n_SCALAR_dofs; }

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

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Returns the first old dof index that is local to subdomain \p proc.
   */
  unsigned int first_old_dof(const unsigned int proc = libMesh::processor_id()) const
  { libmesh_assert(proc < _first_old_df.size()); return _first_old_df[proc]; }
#endif //LIBMESH_ENABLE_AMR

  /**
   * Returns the last dof index that is local to subdomain \p proc.
   * This function is now deprecated, because it returns nonsense in the rare
   * case where \p proc has no local dof indices.  Use end_dof() instead.
   */
  unsigned int last_dof(const unsigned int proc = libMesh::processor_id()) const
  { libmesh_deprecated(); libmesh_assert(proc < _end_df.size()); return (_end_df[proc] - 1); }

  /**
   * Returns the first dof index that is after all indices local to subdomain \p proc.
   * Analogous to the end() member function of STL containers.
   */
  unsigned int end_dof(const unsigned int proc = libMesh::processor_id()) const
  { libmesh_assert(proc < _end_df.size()); return _end_df[proc]; }

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Returns the first old dof index that is after all indices local to subdomain \p proc.
   * Analogous to the end() member function of STL containers.
   */
  unsigned int end_old_dof(const unsigned int proc = libMesh::processor_id()) const
  { libmesh_assert(proc < _end_old_df.size()); return _end_old_df[proc]; }
#endif //LIBMESH_ENABLE_AMR

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
   * Fills the vector \p di with the global degree of freedom indices
   * corresponding to the SCALAR variable vn. If old_dofs=true,
   * the old SCALAR dof indices are returned. Note that we do not
   * need to pass in an element since SCALARs are global variables.
   */
  void SCALAR_dof_indices (std::vector<unsigned int>& di,
		           const unsigned int vn,
		           const bool old_dofs=false) const;

  /**
   * Returns \p true iff all degree of freedom indices in
   * \p dof_indices are either local indices or in the \p send_list.
   * Note that this is an O(logN) operation, not O(1); we don't cache
   * enough information for O(1) right now.
   */
  bool all_semilocal_indices (const std::vector<unsigned int>& dof_indices) const;

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

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * @returns the total number of constrained Nodes
   * in the mesh.
   */
  unsigned int n_constrained_nodes() const { return _node_constraints.size(); }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  /**
   * Rebuilds the raw degree of freedom and DofObject constraints.
   */
  void create_dof_constraints (const MeshBase&);

  /**
   * Gathers any relevant constraint equations from other processors
   */
  void allgather_recursive_constraints (const MeshBase&);

  /**
   * Postprocesses any constrained degrees of freedom
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
   * Returns an iterator pointing to the first DoF constraint row
   */
  DofConstraints::const_iterator constraint_rows_begin() const
    { return _dof_constraints.begin(); }

  /**
   * Returns an iterator pointing just past the last DoF constraint row
   */
  DofConstraints::const_iterator constraint_rows_end() const
    { return _dof_constraints.end(); }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * Returns an iterator pointing to the first Node constraint row
   */
  NodeConstraints::const_iterator node_constraint_rows_begin() const
    { return _node_constraints.begin(); }

  /**
   * Returns an iterator pointing just past the last Node constraint row
   */
  NodeConstraints::const_iterator node_constraint_rows_end() const
    { return _node_constraints.end(); }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  /**
   * @returns true if the degree of freedom dof is constrained,
   * false otherwise.
   */
  bool is_constrained_dof (const unsigned int dof) const;

  /**
   * @returns true if the Node is constrained,
   * false otherwise.
   */
  bool is_constrained_node (const Node* node) const;

  /**
   * Prints the whole \p _dof_constraints and \p
   * _node_constraints data structures.
   */
  void print_dof_constraints(std::ostream& os=libMesh::out) const;

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
   * Does not actually constrain anything, but modifies \p dofs in the
   * same way as any of the constrain functions would do, i.e. adds
   * those dofs in terms of which any of the existing dofs is
   * constrained.
   */
  void constrain_nothing (std::vector<unsigned int>& dofs) const;

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
   * Add a periodic boundary pair
   *
   * @param boundary - primary boundary
   * @param inverse_boundary - inverse boundary
   */
  void add_periodic_boundary (PeriodicBoundary * boundary, PeriodicBoundary * inverse_boundary);

  /**
   * @returns true if the boundary given by \p boundaryid is periodic,
   * false otherwise
   */
  bool is_periodic_boundary (const unsigned int boundaryid) const;

  PeriodicBoundaries * get_periodic_boundaries()
    {
      return _periodic_boundaries;
    }

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
   * Prints summary info about the sparsity bandwidth and constraints.
   */
  void print_info(std::ostream& os=libMesh::out) const;

  /**
   * Gets summary info about the sparsity bandwidth and constraints.
   */
  std::string get_info() const;

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
   * Invalidates all active DofObject dofs for this system
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
   * Finds all the DofObjects associated with the set in \p objs.
   * This will account for off-element couplings via hanging nodes.
   */
  void find_connected_dof_objects (std::vector<const DofObject *> &objs) const;

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
  std::vector<Variable> _variables;

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
   *
   */
  AugmentSparsityPattern *_augment_sparsity_pattern;

  /**
   * A function pointer to a function to call to add extra entries to the sparsity pattern
   */
  void (*_extra_sparsity_function)(SparsityPattern::Graph &,
                                   std::vector<unsigned int> & n_nz,
                                   std::vector<unsigned int> & n_oz,
                                   void *);
  /**
   * A pointer associcated with the extra sparsity that can optionally be passed in
   */
  void * _extra_sparsity_context;

  /**
   * 
   */
  AugmentSendList *_augment_send_list;

  /**
   * A function pointer to a function to call to add extra entries to the send list
   */
  void (*_extra_send_list_function)(std::vector<unsigned int> &, void *);

  /**
   * A pointer associcated with the extra send list that can optionally be passed in
   */
  void * _extra_send_list_context;

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

  /**
   * The total number of SCALAR dofs associated to
   * all SCALAR variables.
   */
  unsigned int _n_SCALAR_dofs;

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Total number of degrees of freedom on old dof objects
   */
  unsigned int _n_old_dfs;

  /**
   * First old DOF index on processor \p p.
   */
  std::vector<unsigned int> _first_old_df;

  /**
   * Last old DOF index (plus 1) on processor \p p.
   */
  std::vector<unsigned int> _end_old_df;

#endif

#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)

  /**
   * Data structure containing DOF constraints.  The ith
   * entry is the constraint matrix row for DOF i.
   */
  DofConstraints _dof_constraints;

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * Data structure containing DofObject constraints.
   */
  NodeConstraints _node_constraints;
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#endif


#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * Data structure containing periodic boundaries.  The ith
   * entry is the constraint matrix row for boundaryid i.
   */
  PeriodicBoundaries *_periodic_boundaries;
#endif

  friend class SparsityPattern::Build;
};


// ------------------------------------------------------------
// Dof Map inline member functions

#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)

inline
bool DofMap::is_constrained_dof (const unsigned int dof) const
{
  if (_dof_constraints.count(dof))
    return true;

  return false;
}

inline
bool DofMap::is_constrained_node (const Node*
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
node
#endif
) const
{
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  if (_node_constraints.count(node))
    return true;
#endif

  return false;
}

#endif


inline
unsigned int DofMap::sys_number() const
{
  return _sys_number;
}


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

} // namespace libMesh

#endif // __dof_map_h__
