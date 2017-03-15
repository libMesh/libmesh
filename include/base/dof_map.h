// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DOF_MAP_H
#define LIBMESH_DOF_MAP_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/enum_order.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/variable.h"
#include "libmesh/threads.h"
#include "libmesh/threads_allocators.h"
#include "libmesh/elem_range.h"
#include "libmesh/ghosting_functor.h"
#include "libmesh/sparsity_pattern.h"
#include "libmesh/parallel_object.h"
#include "libmesh/point.h"

// C++ Includes
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <map>
#include <string>
#include <vector>

namespace libMesh
{

// Forward Declarations
class CouplingMatrix;
class DefaultCoupling;
class DirichletBoundary;
class DirichletBoundaries;
class DofMap;
class DofObject;
class Elem;
class FEType;
class MeshBase;
class PeriodicBoundaryBase;
class PeriodicBoundaries;
class System;
template <typename T> class DenseVectorBase;
template <typename T> class DenseVector;
template <typename T> class DenseMatrix;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;



// ------------------------------------------------------------
// Do we need constraints for anything?

#if defined(LIBMESH_ENABLE_AMR) ||              \
  defined(LIBMESH_ENABLE_PERIODIC) ||           \
  defined(LIBMESH_ENABLE_DIRICHLET)
#  define LIBMESH_ENABLE_CONSTRAINTS 1
#endif

// ------------------------------------------------------------
// AMR constraint matrix types

#ifdef LIBMESH_ENABLE_CONSTRAINTS
/**
 * A row of the Dof constraint matrix.
 */
typedef std::map<dof_id_type, Real,
                 std::less<dof_id_type>,
                 Threads::scalable_allocator<std::pair<const dof_id_type, Real> > > DofConstraintRow;

/**
 * The constraint matrix storage format.
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.  Don't delete this from
 * a pointer-to-std-map; the destructor isn't virtual!
 */
class DofConstraints : public std::map<dof_id_type,
                                       DofConstraintRow,
                                       std::less<dof_id_type>,
                                       Threads::scalable_allocator<std::pair<const dof_id_type, DofConstraintRow> > >
{
};

/**
 * Storage for DofConstraint right hand sides for a particular
 * problem.  Each dof id with a non-zero constraint offset
 * stores it in such a structure.
 */
class DofConstraintValueMap :
    public std::map<dof_id_type, Number,
                    std::less<dof_id_type>,
                    Threads::scalable_allocator<std::pair<const dof_id_type, Number> > >
{
};

/**
 * Storage for DofConstraint right hand sides for all adjoint
 * problems.
 */
class AdjointDofConstraintValues :
    public std::map<unsigned int, DofConstraintValueMap,
                    std::less<unsigned int>,
                    Threads::scalable_allocator
                    <std::pair<const unsigned int, DofConstraintValueMap> > >
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
 * a pointer-to-std-map; the destructor isn't virtual!
 */
class NodeConstraints : public std::map<const Node *,
                                        std::pair<NodeConstraintRow,Point>,
                                        std::less<const Node *>,
                                        Threads::scalable_allocator<std::pair<const Node * const, std::pair<NodeConstraintRow,Point> > > >
{
};
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#endif // LIBMESH_ENABLE_CONSTRAINTS



/**
 * This class handles the numbering of degrees of freedom on a mesh.
 * For systems of equations the class supports a fixed number of variables.
 * The degrees of freedom are numbered such that sequential, contiguous blocks
 * belong to distinct processors.  This is so that the resulting data
 * structures will work well with parallel linear algebra packages.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 * \brief Manages the degrees of freedom (DOFs) in a simulation.
 */
class DofMap : public ReferenceCountedObject<DofMap>,
               public ParallelObject
{
public:

  /**
   * Constructor.  Requires the number of the system for which we
   * will be numbering degrees of freedom & the parent object
   * we are contained in, which defines our communication space.
   */
  explicit
  DofMap(const unsigned int sys_number,
         MeshBase & mesh);

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
    virtual ~AugmentSparsityPattern () {}

    /**
     * User-defined function to augment the sparsity pattern.
     */
    virtual void augment_sparsity_pattern (SparsityPattern::Graph & sparsity,
                                           std::vector<dof_id_type> & n_nz,
                                           std::vector<dof_id_type> & n_oz) = 0;
  };

  /**
   * Abstract base class to be used to add user-defined parallel
   * degree of freedom couplings.
   */
  class AugmentSendList
  {
  public:
    virtual ~AugmentSendList () {}

    /**
     * User-defined function to augment the send list.
     */
    virtual void augment_send_list (std::vector<dof_id_type> & send_list) = 0;
  };

  /**
   * Additional matrices may be handled with this \p DofMap.
   * They are initialized to the same sparsity structure as
   * the major matrix.
   */
  void attach_matrix (SparseMatrix<Number> & matrix);

  /**
   * Matrices should not be attached more than once.  We can test for
   * an already-attached matrix if necessary using \p is_attached
   */
  bool is_attached (SparseMatrix<Number> & matrix);

  /**
   * Distrubute dofs on the current mesh.  Also builds the send list for
   * processor \p proc_id, which defaults to 0 for ease of use in serial
   * applications.
   */
  void distribute_dofs (MeshBase &);

  /**
   * Computes the sparsity pattern for the matrices corresponding to
   * \p proc_id and sends that data to Linear Algebra packages for
   * preallocation of sparse matrices.
   */
  void compute_sparsity (const MeshBase &);

  /**
   * Clears the sparsity pattern
   */
  void clear_sparsity();

  /**
   * Adds a functor which can specify coupling requirements for
   * creation of sparse matrices.
   * Degree of freedom pairs which match the elements and variables
   * returned by these functors will be added to the sparsity pattern,
   * and the degrees of freedom which live on other processors will be
   * added to the send_list for use on ghosted vectors, and the
   * elements which live on other processors will be ghosted on a
   * distributed mesh.
   *
   * GhostingFunctor memory must be managed by the code which calls
   * this function; the GhostingFunctor lifetime is expected to extend
   * until either the functor is removed or the DofMap is destructed.
   */
  void add_coupling_functor(GhostingFunctor & coupling_functor);

  /**
   * Removes a functor which was previously added to the set of
   * coupling functors.
   */
  void remove_coupling_functor(GhostingFunctor & coupling_functor);

  /**
   * Beginning of range of coupling functors
   */
  std::set<GhostingFunctor *>::const_iterator coupling_functors_begin() const
  { return _coupling_functors.begin(); }

  /**
   * End of range of coupling functors
   */
  std::set<GhostingFunctor *>::const_iterator coupling_functors_end() const
  { return _coupling_functors.end(); }

  /**
   * Default coupling functor
   */
  DefaultCoupling & default_coupling() { return *_default_coupling; }

  /**
   * Adds a functor which can specify algebraic ghosting requirements
   * for use with distributed vectors.  Degrees of freedom on other
   * processors which match the elements and variables returned by
   * these functors will be added to the send_list, and the elements
   * on other processors will be ghosted on a distributed mesh.
   *
   * GhostingFunctor memory must be managed by the code which calls
   * this function; the GhostingFunctor lifetime is expected to extend
   * until either the functor is removed or the DofMap is destructed.
   */
  void add_algebraic_ghosting_functor(GhostingFunctor & ghosting_functor);

  /**
   * Removes a functor which was previously added to the set of
   * algebraic ghosting functors.
   */
  void remove_algebraic_ghosting_functor(GhostingFunctor & ghosting_functor);

  /**
   * Beginning of range of algebraic ghosting functors
   */
  std::set<GhostingFunctor *>::const_iterator algebraic_ghosting_functors_begin() const
  { return _algebraic_ghosting_functors.begin(); }

  /**
   * End of range of algebraic ghosting functors
   */
  std::set<GhostingFunctor *>::const_iterator algebraic_ghosting_functors_end() const
  { return _algebraic_ghosting_functors.end(); }

  /**
   * Default algebraic ghosting functor
   */
  DefaultCoupling & default_algebraic_ghosting() { return *_default_evaluating; }

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
  void attach_extra_sparsity_object (DofMap::AugmentSparsityPattern & asp)
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
                                                   std::vector<dof_id_type> & n_nz,
                                                   std::vector<dof_id_type> & n_oz,
                                                   void *),
                                      void * context = libmesh_nullptr)
  { _extra_sparsity_function = func; _extra_sparsity_context = context; }

  /**
   * Attach an object to populate the send_list with extra entries.
   * This should only add to the send list, but no checking is done
   * to enforce this behavior.
   *
   * This is an advanced function... use at your own peril!
   */
  void attach_extra_send_list_object (DofMap::AugmentSendList & asl)
  {
    _augment_send_list = &asl;
  }

  /**
   * Attach a function pointer to use as a callback to populate the
   * send_list with extra entries.
   */
  void attach_extra_send_list_function(void (*func)(std::vector<dof_id_type> &, void *),
                                       void * context = libmesh_nullptr)
  { _extra_send_list_function = func; _extra_send_list_context = context; }

  /**
   * Takes the \p _send_list vector (which may have duplicate entries)
   * and sorts it.  The duplicate entries are then removed, resulting in
   * a sorted \p _send_list with unique entries.  Also calls any user-provided
   * methods for adding to the send list.
   */
  void prepare_send_list ();

  /**
   * Returns a constant reference to the \p _send_list for this processor.  The
   * \p _send_list contains the global indices of all the variables in the
   * global solution vector that influence the current processor.  This
   * information can be used for gathers at each solution step to retrieve
   * solution values needed for computation.
   */
  const std::vector<dof_id_type> & get_send_list() const { return _send_list; }

  /**
   * Returns a constant reference to the \p _n_nz list for this processor.
   * The vector contains the bandwidth of the on-processor coupling for each
   * row of the global matrix that the current processor owns.  This
   * information can be used to preallocate space for a parallel sparse matrix.
   */
  const std::vector<dof_id_type> & get_n_nz() const
  {
    libmesh_assert(_n_nz);
    return *_n_nz;
  }

  /**
   * Returns a constant reference to the \p _n_oz list for this processor.
   * The vector contains the bandwidth of the off-processor coupling for each
   * row of the global matrix that the current processor owns.  This
   * information can be used to preallocate space for a parallel sparse matrix.
   */
  const std::vector<dof_id_type> & get_n_oz() const
  {
    libmesh_assert(_n_oz);
    return *_n_oz;
  }

  // /**
  //  * Add an unknown of order \p order and finite element type
  //  * \p type to the system of equations.
  //  */
  // void add_variable (const Variable & var);

  /**
   * Add a group of unknowns of order \p order and finite element type
   * \p type to the system of equations.
   */
  void add_variable_group (const VariableGroup & var_group);

  /**
   * @returns the \p VariableGroup description object for group \p g.
   */
  const VariableGroup & variable_group (const unsigned int c) const;

  /**
   * @returns the variable description object for variable \p c.
   */
  const Variable & variable (const unsigned int c) const;

  /**
   * @returns the approximation order for variable \p c.
   */
  Order variable_order (const unsigned int c) const;

  /**
   * @returns the approximation order for \p VariableGroup \p vg.
   */
  Order variable_group_order (const unsigned int vg) const;

  /**
   * @returns the finite element type for variable \p c.
   */
  const FEType & variable_type (const unsigned int c) const;

  /**
   * @returns the finite element type for \p VariableGroup \p vg.
   */
  const FEType & variable_group_type (const unsigned int vg) const;

  /**
   * @returns the number of variables in the global solution vector. Defaults
   * to 1, should be 1 for a scalar equation, 3 for 2D incompressible Navier
   * Stokes (u,v,p), etc...
   */
  unsigned int n_variable_groups() const
  { return cast_int<unsigned int>(_variable_groups.size()); }

  /**
   * @returns the number of variables in the global solution vector. Defaults
   * to 1, should be 1 for a scalar equation, 3 for 2D incompressible Navier
   * Stokes (u,v,p), etc...
   */
  unsigned int n_variables() const
  { return cast_int<unsigned int>(_variables.size()); }

  /**
   * @returns true if the variables are capable of being stored in a blocked
   * form.  Presently, this means that there can only be one variable group,
   * and that the group has more than one variable.
   */
  bool has_blocked_representation() const
  {
#ifdef LIBMESH_ENABLE_BLOCKED_STORAGE
    return ((this->n_variable_groups() == 1) && (this->n_variables() > 1));
#else
    return false;
#endif
  }

  /**
   * @returns the block size, if the variables are amenable to block storage.
   * Otherwise 1.
   */
  unsigned int block_size() const
  {
#ifdef LIBMESH_ENABLE_BLOCKED_STORAGE
    return (this->has_blocked_representation() ? this->n_variables() : 1);
#else
    return 1;
#endif
  }

  /**
   * @returns the total number of degrees of freedom in the problem.
   */
  dof_id_type n_dofs() const { return _n_dfs; }

  /**
   * @returns the number of SCALAR dofs.
   */
  dof_id_type n_SCALAR_dofs() const { return _n_SCALAR_dofs; }

  /**
   * @returns the number of degrees of freedom on this processor.
   */
  dof_id_type n_local_dofs () const
  { return this->n_dofs_on_processor (this->processor_id()); }

  /**
   * Returns the number of degrees of freedom on partition \p proc.
   */
  dof_id_type n_dofs_on_processor(const processor_id_type proc) const
  {
    libmesh_assert_less (proc, _first_df.size());
    return cast_int<dof_id_type>(_end_df[proc] - _first_df[proc]);
  }

  /**
   * Returns the first dof index that is local to partition \p proc.
   */
  dof_id_type first_dof(const processor_id_type proc) const
  { libmesh_assert_less (proc, _first_df.size()); return _first_df[proc]; }

  dof_id_type first_dof() const
  { return this->first_dof(this->processor_id()); }

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Returns the first old dof index that is local to partition \p proc.
   */
  dof_id_type first_old_dof(const processor_id_type proc) const
  { libmesh_assert_less (proc, _first_old_df.size()); return _first_old_df[proc]; }

  dof_id_type first_old_dof() const
  { return this->first_old_dof(this->processor_id()); }

#endif //LIBMESH_ENABLE_AMR

  /**
   * Returns the last dof index that is local to processor \p proc.
   * This function is now deprecated, because it returns nonsense in the rare
   * case where \p proc has no local dof indices.  Use end_dof() instead.
   */
  dof_id_type last_dof(const processor_id_type proc) const
  {
    libmesh_deprecated();
    libmesh_assert_less (proc, _end_df.size());
    return cast_int<dof_id_type>(_end_df[proc] - 1);
  }

  dof_id_type last_dof() const
  { return this->last_dof(this->processor_id()); }

  /**
   * Returns the first dof index that is after all indices local to
   * processor \p proc.
   * Analogous to the end() member function of STL containers.
   */
  dof_id_type end_dof(const processor_id_type proc) const
  { libmesh_assert_less (proc, _end_df.size()); return _end_df[proc]; }

  dof_id_type end_dof() const
  { return this->end_dof(this->processor_id()); }

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Returns the first old dof index that is after all indices local
   * to processor \p proc.
   * Analogous to the end() member function of STL containers.
   */
  dof_id_type end_old_dof(const processor_id_type proc) const
  { libmesh_assert_less (proc, _end_old_df.size()); return _end_old_df[proc]; }

  dof_id_type end_old_dof() const
  { return this->end_old_dof(this->processor_id()); }

#endif //LIBMESH_ENABLE_AMR

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the element.
   */
  void dof_indices (const Elem * const elem,
                    std::vector<dof_id_type> & di) const;

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the element.  For one variable, and potentially for a
   * non-default element p refinement level
   */
  void dof_indices (const Elem * const elem,
                    std::vector<dof_id_type> & di,
                    const unsigned int vn,
                    int p_level = -12345) const;

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the node.
   */
  void dof_indices (const Node * const node,
                    std::vector<dof_id_type> & di) const;

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the node.   For one variable \p vn.
   */
  void dof_indices (const Node * const node,
                    std::vector<dof_id_type> & di,
                    const unsigned int vn) const;

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * corresponding to the SCALAR variable vn. If old_dofs=true,
   * the old SCALAR dof indices are returned. Note that we do not
   * need to pass in an element since SCALARs are global variables.
   */
  void SCALAR_dof_indices (std::vector<dof_id_type> & di,
                           const unsigned int vn,
                           const bool old_dofs=false) const;

  /**
   * Returns \p true iff all degree of freedom indices in
   * \p dof_indices are either local indices or in the \p send_list.
   * Note that this is an O(logN) operation, not O(1); we don't cache
   * enough information for O(1) right now.
   */
  bool all_semilocal_indices (const std::vector<dof_id_type> & dof_indices) const;

  /**
   * @returns \p true iff our solutions can be locally evaluated on
   * \p obj (which should be an Elem or a Node) for variable number \p
   * var_num (for all variables, if \p var_num is invalid_uint)
   */
  template <typename DofObjectSubclass>
  bool is_evaluable(const DofObjectSubclass & obj,
                    unsigned int var_num = libMesh::invalid_uint) const;

  /**
   * Allow the implicit_neighbor_dofs flag to be set programmatically.
   * This overrides the --implicit_neighbor_dofs commandline option.
   * We can use this to set the implicit neighbor dofs option differently
   * for different systems, whereas the commandline option is the same
   * for all systems.
   */
  void set_implicit_neighbor_dofs(bool implicit_neighbor_dofs);

  /**
   * Tells other library functions whether or not this problem
   * includes coupling between dofs in neighboring cells, as can
   * currently be specified on the command line or inferred from
   * the use of all discontinuous variables.
   */
  bool use_coupled_neighbor_dofs(const MeshBase & mesh) const;

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
  void extract_local_vector (const NumericVector<Number> & Ug,
                             const std::vector<dof_id_type> & dof_indices,
                             DenseVectorBase<Number> & Ue) const;

  /**
   * Fills an array of those dof indices which belong to the given
   * variable number and live on the current processor.
   */
  void local_variable_indices(std::vector<dof_id_type> & idx,
                              const MeshBase & mesh,
                              unsigned int var_num) const;

#ifdef LIBMESH_ENABLE_CONSTRAINTS

  //--------------------------------------------------------------------
  // Constraint-specific methods
  /**
   * @returns the total number of constrained degrees of freedom
   * in the problem.
   */
  dof_id_type n_constrained_dofs() const;

  /**
   * @returns the number of constrained degrees of freedom
   * on this processor.
   */
  dof_id_type n_local_constrained_dofs() const;

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * @returns the total number of constrained Nodes
   * in the mesh.
   */
  dof_id_type n_constrained_nodes() const
  { return cast_int<dof_id_type>(_node_constraints.size()); }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  /**
   * Rebuilds the raw degree of freedom and DofObject constraints.
   * A time is specified for use in building time-dependent Dirichlet
   * constraints.
   */
  void create_dof_constraints (const MeshBase &, Real time=0);

  /**
   * Gathers constraint equation dependencies from other processors
   */
  void allgather_recursive_constraints (MeshBase &);

  /**
   * Sends constraint equations to constraining processors
   */
  void scatter_constraints (MeshBase &);

  /**
   * Helper function for querying about constraint equations on other
   * processors.  If any id in \p requested_dof_ids is constrained on
   * another processor, its constraint will be added on this processor
   * as well.  If \p look_for_constrainees is true, then constraints
   * will also be returned if the id appears as a constraining value
   * not just if it appears as a constrained value.
   *
   * This function operates recursively: if the constraint for a
   * constrained dof is newly added locally, then any other dofs which
   * constrain it are queried to see if they are in turn constrained,
   * and so on.
   */
  void gather_constraints (MeshBase & mesh,
                           std::set<dof_id_type> & unexpanded_dofs,
                           bool look_for_constrainees);

  /**
   * Postprocesses any constrained degrees of freedom
   * to be constrained only in terms of unconstrained dofs, then adds
   * unconstrained dofs to the send_list and prepares that for use.
   * This should be run after both system (create_dof_constraints) and
   * user constraints have all been added.
   */
  void process_constraints (MeshBase &);

  /**
   * Adds a copy of the user-defined row to the constraint matrix, using
   * an inhomogeneous right-hand-side for the constraint equation.
   */
  void add_constraint_row (const dof_id_type dof_number,
                           const DofConstraintRow & constraint_row,
                           const Number constraint_rhs,
                           const bool forbid_constraint_overwrite);

  /**
   * Adds a copy of the user-defined row to the constraint matrix,
   * using an inhomogeneous right-hand-side for the adjoint constraint
   * equation.
   *
   * \p forbid_constraint_overwrite here only tests for overwriting
   * the rhs.  This method should only be used when an equivalent
   * constraint (with a potentially different rhs) already exists for
   * the primal problem.
   */
  void add_adjoint_constraint_row (const unsigned int qoi_index,
                                   const dof_id_type dof_number,
                                   const DofConstraintRow & constraint_row,
                                   const Number constraint_rhs,
                                   const bool forbid_constraint_overwrite);

  /**
   * Adds a copy of the user-defined row to the constraint matrix, using
   * a homogeneous right-hand-side for the constraint equation.
   * By default, produces an error if the DOF was already constrained.
   */
  void add_constraint_row (const dof_id_type dof_number,
                           const DofConstraintRow & constraint_row,
                           const bool forbid_constraint_overwrite = true)
  { add_constraint_row(dof_number, constraint_row, 0., forbid_constraint_overwrite); }

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

  void stash_dof_constraints()
  {
    libmesh_assert(_stashed_dof_constraints.empty());
    _dof_constraints.swap(_stashed_dof_constraints);
  }

  void unstash_dof_constraints()
  {
    libmesh_assert(_dof_constraints.empty());
    _dof_constraints.swap(_stashed_dof_constraints);
  }

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
  bool is_constrained_dof (const dof_id_type dof) const;

  /**
   * @returns true if the system has any heterogenous constraints for
   * adjoint solution \p qoi_num, false otherwise.
   */
  bool has_heterogenous_adjoint_constraints (const unsigned int qoi_num) const;

  /**
   * @returns the heterogeneous constraint value if the degree of
   * freedom \p dof has a heterogenous constraint for adjoint solution
   * \p qoi_num, zero otherwise.
   */
  Number has_heterogenous_adjoint_constraint (const unsigned int qoi_num,
                                              const dof_id_type dof) const;

  /**
   * @returns a reference to the set of right-hand-side values in
   * primal constraint equations
   */
  DofConstraintValueMap & get_primal_constraint_values();

  /**
   * @returns true if the Node is constrained,
   * false otherwise.
   */
  bool is_constrained_node (const Node * node) const;

  /**
   * Prints (from processor 0) all DoF and Node constraints.  If \p
   * print_nonlocal is true, then each constraint is printed once for
   * each processor that knows about it, which may be useful for \p
   * DistributedMesh debugging.
   */
  void print_dof_constraints(std::ostream & os=libMesh::out,
                             bool print_nonlocal=false) const;

  /**
   * Gets a string reporting all DoF and Node constraints local to
   * this processor.  If \p print_nonlocal is true, then nonlocal
   * constraints which are locally known are included.
   */
  std::string get_local_constraints(bool print_nonlocal=false) const;


  /**
   * Tests the constrained degrees of freedom on the numeric vector \p v, which
   * represents a solution defined on the mesh, returning a pair whose first
   * entry is the maximum absolute error on a constrained DoF and whose second
   * entry is the maximum relative error.  Useful for debugging purposes.
   *
   * If \p v == libmesh_nullptr, the system solution vector is tested.
   */
  std::pair<Real, Real> max_constraint_error(const System & system,
                                             NumericVector<Number> * v = libmesh_nullptr) const;

#endif // LIBMESH_ENABLE_CONSTRAINTS

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
  void constrain_element_matrix (DenseMatrix<Number> & matrix,
                                 std::vector<dof_id_type> & elem_dofs,
                                 bool asymmetric_constraint_rows = true) const;

  /**
   * Constrains the element matrix.  This method allows the
   * element matrix to be non-square, in which case the row_dofs
   * and col_dofs may be of different size and correspond to
   * variables approximated in different spaces.
   */
  void constrain_element_matrix (DenseMatrix<Number> & matrix,
                                 std::vector<dof_id_type> & row_dofs,
                                 std::vector<dof_id_type> & col_dofs,
                                 bool asymmetric_constraint_rows = true) const;

  /**
   * Constrains the element vector.
   */
  void constrain_element_vector (DenseVector<Number> & rhs,
                                 std::vector<dof_id_type> & dofs,
                                 bool asymmetric_constraint_rows = true) const;

  /**
   * Constrains the element matrix and vector.  This method requires
   * the element matrix to be square, in which case the elem_dofs
   * correspond to the global DOF indices of both the rows and
   * columns of the element matrix.  For this case the rows
   * and columns of the matrix necessarily correspond to variables
   * of the same approximation order.
   */
  void constrain_element_matrix_and_vector (DenseMatrix<Number> & matrix,
                                            DenseVector<Number> & rhs,
                                            std::vector<dof_id_type> & elem_dofs,
                                            bool asymmetric_constraint_rows = true) const;

  /**
   * Constrains the element matrix and vector.  This method requires
   * the element matrix to be square, in which case the elem_dofs
   * correspond to the global DOF indices of both the rows and
   * columns of the element matrix.  For this case the rows
   * and columns of the matrix necessarily correspond to variables
   * of the same approximation order.
   *
   * The heterogenous version of this method creates linear systems in
   * which heterogenously constrained degrees of freedom will solve to
   * their correct offset values, as would be appropriate for finding
   * a solution to a linear problem in a single algebraic solve.  The
   * non-heterogenous version of this method creates linear systems in
   * which even heterogenously constrained degrees of freedom are
   * solved without offset values taken into account, as would be
   * appropriate for finding linearized updates to a solution in which
   * heterogenous constraints are already satisfied.
   *
   * By default, the constraints for the primal solution of this
   * system are used.  If a non-negative \p qoi_index is passed in,
   * then the constraints for the corresponding adjoint solution are
   * used instead.
   */
  void heterogenously_constrain_element_matrix_and_vector (DenseMatrix<Number> & matrix,
                                                           DenseVector<Number> & rhs,
                                                           std::vector<dof_id_type> & elem_dofs,
                                                           bool asymmetric_constraint_rows = true,
                                                           int qoi_index = -1) const;

  /**
   * Constrains the element vector.  This method requires
   * the element matrix to be square and not-yet-constrained, in which
   * case the elem_dofs correspond to the global DOF indices of both
   * the rows and columns of the element matrix.
   *
   * The heterogenous version of this method creates linear systems in
   * which heterogenously constrained degrees of freedom will solve to
   * their correct offset values, as would be appropriate for finding
   * a solution to a linear problem in a single algebraic solve.  The
   * non-heterogenous version of this method creates linear systems in
   * which even heterogenously constrained degrees of freedom are
   * solved without offset values taken into account, as would be
   * appropriate for finding linearized updates to a solution in which
   * heterogenous constraints are already satisfied.
   *
   * By default, the constraints for the primal solution of this
   * system are used.  If a non-negative \p qoi_index is passed in,
   * then the constraints for the corresponding adjoint solution are
   * used instead.
   */
  void heterogenously_constrain_element_vector (const DenseMatrix<Number> & matrix,
                                                DenseVector<Number> & rhs,
                                                std::vector<dof_id_type> & elem_dofs,
                                                bool asymmetric_constraint_rows = true,
                                                int qoi_index = -1) const;



  /**
   * Constrains a dyadic element matrix B = v w'.  This method
   * requires the element matrix to be square, in which case the
   * elem_dofs correspond to the global DOF indices of both the rows
   * and columns of the element matrix.  For this case the rows and
   * columns of the matrix necessarily correspond to variables of the
   * same approximation order.
   */
  void constrain_element_dyad_matrix (DenseVector<Number> & v,
                                      DenseVector<Number> & w,
                                      std::vector<dof_id_type> & row_dofs,
                                      bool asymmetric_constraint_rows = true) const;

  /**
   * Does not actually constrain anything, but modifies \p dofs in the
   * same way as any of the constrain functions would do, i.e. adds
   * those dofs in terms of which any of the existing dofs is
   * constrained.
   */
  void constrain_nothing (std::vector<dof_id_type> & dofs) const;

  /**
   * Constrains the numeric vector \p v, which represents a solution defined on
   * the mesh.  This may need to be used after a linear solve, if your linear
   * solver's solutions do not satisfy your DoF constraints to a tight enough
   * tolerance.
   *
   * If \p v == libmesh_nullptr, the system solution vector is constrained
   *
   * If \p homogeneous == true, heterogeneous constraints are enforced
   * as if they were homogeneous.  This might be appropriate for e.g. a
   * vector representing a difference between two
   * heterogeneously-constrained solutions.
   */
  void enforce_constraints_exactly (const System & system,
                                    NumericVector<Number> * v = libmesh_nullptr,
                                    bool homogeneous = false) const;

  /**
   * Heterogenously constrains the numeric vector \p v, which
   * represents an adjoint solution defined on the mesh for quantity
   * fo interest \p q.  For homogeneous constraints, use \p
   * enforce_constraints_exactly instead
   */
  void enforce_adjoint_constraints_exactly (NumericVector<Number> & v,
                                            unsigned int q) const;



#ifdef LIBMESH_ENABLE_PERIODIC

  //--------------------------------------------------------------------
  // PeriodicBoundary-specific methods

  /**
   * Adds a copy of the specified periodic boundary to the system.
   */
  void add_periodic_boundary (const PeriodicBoundaryBase & periodic_boundary);

  /**
   * Add a periodic boundary pair
   *
   * @param boundary - primary boundary
   * @param inverse_boundary - inverse boundary
   */
  void add_periodic_boundary (const PeriodicBoundaryBase & boundary, const PeriodicBoundaryBase & inverse_boundary);

  /**
   * @returns true if the boundary given by \p boundaryid is periodic,
   * false otherwise
   */
  bool is_periodic_boundary (const boundary_id_type boundaryid) const;

  PeriodicBoundaries * get_periodic_boundaries()
  {
    return _periodic_boundaries.get();
  }

#endif // LIBMESH_ENABLE_PERIODIC


#ifdef LIBMESH_ENABLE_DIRICHLET

  //--------------------------------------------------------------------
  // DirichletBoundary-specific methods

  /**
   * Adds a copy of the specified Dirichlet boundary to the system.
   */
  void add_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary);

  /**
   * Adds a copy of the specified Dirichlet boundary to the system,
   * corresponding to the adjoint problem defined by Quantity of
   * Interest \p q.
   */
  void add_adjoint_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary,
                                       unsigned int q);

  /**
   * Removes the specified Dirichlet boundary from the system.
   */
  void remove_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary);

  /**
   * Removes from the system the specified Dirichlet boundary for the
   * adjoint equation defined by Quantity of interest index q
   */
  void remove_adjoint_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary,
                                          unsigned int q);

  const DirichletBoundaries * get_dirichlet_boundaries() const
  {
    return _dirichlet_boundaries.get();
  }

  DirichletBoundaries * get_dirichlet_boundaries()
  {
    return _dirichlet_boundaries.get();
  }

  bool has_adjoint_dirichlet_boundaries(unsigned int q) const;

  const DirichletBoundaries *
  get_adjoint_dirichlet_boundaries(unsigned int q) const;

  DirichletBoundaries *
  get_adjoint_dirichlet_boundaries(unsigned int q);

#endif // LIBMESH_ENABLE_DIRICHLET


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
  void old_dof_indices (const Elem * const elem,
                        std::vector<dof_id_type> & di,
                        const unsigned int vn = libMesh::invalid_uint) const;
  /**
   * @returns the total number of degrees of freedom on old_dof_objects
   *
   */
  dof_id_type n_old_dofs() const { return _n_old_dfs; }

  /**
   * Constrains degrees of freedom on side \p s of element \p elem which
   * correspond to variable number \p var and to p refinement levels
   * above \p p.
   */
  void constrain_p_dofs (unsigned int var,
                         const Elem * elem,
                         unsigned int s,
                         unsigned int p);

#endif // LIBMESH_ENABLE_AMR

  /**
   * Reinitialize the underlying data strucures conformal to the current mesh.
   */
  void reinit (MeshBase & mesh);

  /**
   * Free all memory associated with the object, but keep the mesh pointer.
   */
  void clear ();

  /**
   * Prints summary info about the sparsity bandwidth and constraints.
   */
  void print_info(std::ostream & os=libMesh::out) const;

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
  CouplingMatrix * _dof_coupling;

  /**
   * @returns the number of the system we are responsible for.
   */
  unsigned int sys_number() const;

private:

  /**
   * Helper function that gets the dof indices on the current element
   * for a non-SCALAR type variable.
   *
   * In DEBUG mode, the tot_size parameter will add up the total
   * number of dof indices that should have been added to di.
   */
  void _dof_indices (const Elem * const elem,
                     int p_level,
                     std::vector<dof_id_type> & di,
                     const unsigned int v,
                     const Node * const * nodes,
                     unsigned int       n_nodes
#ifdef DEBUG
                     ,
                     std::size_t & tot_size
#endif
                     ) const;

  /**
   * Builds a sparsity pattern
   */
  UniquePtr<SparsityPattern::Build> build_sparsity(const MeshBase & mesh) const;

  /**
   * Invalidates all active DofObject dofs for this system
   */
  void invalidate_dofs(MeshBase & mesh) const;

  /**
   * An adapter function that returns Node pointers by index
   */
  DofObject * node_ptr(MeshBase & mesh, dof_id_type i) const;

  /**
   * An adapter function that returns Elem pointers by index
   */
  DofObject * elem_ptr(MeshBase & mesh, dof_id_type i) const;

  /**
   * A member function type like node_ptr or elem_ptr
   */
  typedef DofObject * (DofMap::*dofobject_accessor)
    (MeshBase & mesh, dof_id_type i) const;

  /**
   * Helper function for distributing dofs in parallel
   */
  template<typename iterator_type>
  void set_nonlocal_dof_objects(iterator_type objects_begin,
                                iterator_type objects_end,
                                MeshBase & mesh,
                                dofobject_accessor objects);

  /**
   * Distributes the global degrees of freedom, for dofs on
   * this processor.  In this format the local
   * degrees of freedom are in a contiguous block for each
   * variable in the system.
   * Starts at index next_free_dof, and increments it to
   * the post-final index.
   */
  void distribute_local_dofs_var_major (dof_id_type & next_free_dof,
                                        MeshBase & mesh);

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
  void distribute_local_dofs_node_major (dof_id_type & next_free_dof,
                                         MeshBase & mesh);

  /*
   * A utility method for obtaining a set of elements to ghost along
   * with merged coupling matrices.
   */
  static void
  merge_ghost_functor_outputs (GhostingFunctor::map_type & elements_to_ghost,
                               std::set<CouplingMatrix *> & temporary_coupling_matrices,
                               const std::set<GhostingFunctor *>::iterator & gf_begin,
                               const std::set<GhostingFunctor *>::iterator & gf_end,
                               const MeshBase::const_element_iterator & elems_begin,
                               const MeshBase::const_element_iterator & elems_end,
                               processor_id_type p);

  /**
   * Adds entries to the \p _send_list vector corresponding to DoFs
   * on elements neighboring the current processor.
   */
  void add_neighbors_to_send_list(MeshBase & mesh);

#ifdef LIBMESH_ENABLE_CONSTRAINTS

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
  void build_constraint_matrix (DenseMatrix<Number> & C,
                                std::vector<dof_id_type> & elem_dofs,
                                const bool called_recursively=false) const;

  /**
   * Build the constraint matrix C and the forcing vector H
   * associated with the element degree of freedom indices elem_dofs.
   * The optional parameter \p called_recursively should be left at
   * the default value \p false.  This is used to handle the special
   * case of an element's degrees of freedom being constrained in
   * terms of other, local degrees of freedom.  The usual case is for
   * an elements DOFs to be constrained by some other, external DOFs
   * and/or Dirichlet conditions.
   *
   * The forcing vector will depend on which solution's heterogenous
   * constraints are being applied.  For the default \p qoi_index this
   * will be the primal solutoin; for \p qoi_index >= 0 the
   * corresponding adjoint solution's constraints will be used.
   */
  void build_constraint_matrix_and_vector (DenseMatrix<Number> & C,
                                           DenseVector<Number> & H,
                                           std::vector<dof_id_type> & elem_dofs,
                                           int qoi_index = -1,
                                           const bool called_recursively=false) const;

  /**
   * Finds all the DOFS associated with the element DOFs elem_dofs.
   * This will account for off-element couplings via hanging nodes.
   */
  void find_connected_dofs (std::vector<dof_id_type> & elem_dofs) const;

  /**
   * Finds all the DofObjects associated with the set in \p objs.
   * This will account for off-element couplings via hanging nodes.
   */
  void find_connected_dof_objects (std::vector<const DofObject *> & objs) const;

  /**
   * Adds entries to the \p _send_list vector corresponding to DoFs
   * which are dependencies for constraint equations on the current
   * processor.
   */
  void add_constraints_to_send_list();

#endif // LIBMESH_ENABLE_CONSTRAINTS

  /**
   * The finite element type for each variable.
   */
  std::vector<Variable> _variables;

  /**
   * The finite element type for each variable.
   */
  std::vector<VariableGroup> _variable_groups;

  /**
   * The number of the system we manage DOFs for.
   */
  const unsigned int _sys_number;

  /**
   * The mesh that system uses.
   */
  MeshBase & _mesh;

  /**
   * Additional matrices handled by this object.  These pointers do @e
   * not handle the memory, instead, \p System, who
   * told \p DofMap about them, owns them.
   */
  std::vector<SparseMatrix<Number> * > _matrices;

  /**
   * First DOF index on processor \p p.
   */
  std::vector<dof_id_type> _first_df;

  /**
   * Last DOF index (plus 1) on processor \p p.
   */
  std::vector<dof_id_type> _end_df;

  /**
   * First DOF index for SCALAR variable v, or garbage for non-SCALAR
   * variable v
   */
  std::vector<dof_id_type> _first_scalar_df;

  /**
   * A list containing all the global DOF indices that affect the
   * solution on my processor.
   */
  std::vector<dof_id_type> _send_list;

  /**
   * Funtion object to call to add extra entries to the sparsity pattern
   */
  AugmentSparsityPattern * _augment_sparsity_pattern;

  /**
   * A function pointer to a function to call to add extra entries to the sparsity pattern
   */
  void (*_extra_sparsity_function)(SparsityPattern::Graph &,
                                   std::vector<dof_id_type> & n_nz,
                                   std::vector<dof_id_type> & n_oz,
                                   void *);
  /**
   * A pointer associcated with the extra sparsity that can optionally be passed in
   */
  void * _extra_sparsity_context;

  /**
   * Function object to call to add extra entries to the send list
   */
  AugmentSendList * _augment_send_list;

  /**
   * A function pointer to a function to call to add extra entries to the send list
   */
  void (*_extra_send_list_function)(std::vector<dof_id_type> &, void *);

  /**
   * A pointer associcated with the extra send list that can optionally be passed in
   */
  void * _extra_send_list_context;

  /**
   * The default coupling GhostingFunctor, used to implement standard
   * libMesh sparsity pattern construction.
   *
   * We use a UniquePtr here to reduce header dependencies.
   */
  UniquePtr<DefaultCoupling> _default_coupling;

  /**
   * The default algebraic GhostingFunctor, used to implement standard
   * libMesh send_list construction.
   *
   * We use a UniquePtr here to reduce header dependencies.
   */
  UniquePtr<DefaultCoupling> _default_evaluating;

  /**
   * The list of all GhostingFunctor objects to be used when
   * distributing ghosted vectors.
   *
   * The library should automatically copy these functors to the
   * MeshBase, too, so any algebraically ghosted dofs will live on
   * geometrically ghosted elements.
   */
  std::set<GhostingFunctor *> _algebraic_ghosting_functors;

  /**
   * The list of all GhostingFunctor objects to be used when
   * coupling degrees of freedom in matrix sparsity patterns.
   *
   * These objects will *also* be used as algebraic ghosting functors,
   * but not vice-versa.
   *
   * The library should automatically copy these functors to the
   * MeshBase, too, so any dofs coupled to local dofs will live on
   * geometrically ghosted elements.
   */
  std::set<GhostingFunctor *> _coupling_functors;

  /**
   * Default false; set to true if any attached matrix requires a full
   * sparsity pattern.
   */
  bool need_full_sparsity_pattern;

  /**
   * The sparsity pattern of the global matrix, kept around if it
   * might be needed by future additions of the same type of matrix.
   */
  UniquePtr<SparsityPattern::Build> _sp;

  /**
   * The number of on-processor nonzeros in my portion of the
   * global matrix.  If need_full_sparsity_pattern is true, this will
   * just be a pointer into the corresponding sparsity pattern vector.
   * Otherwise we have to new/delete it ourselves.
   */
  std::vector<dof_id_type> * _n_nz;

  /**
   * The number of off-processor nonzeros in my portion of the
   * global matrix; allocated similar to _n_nz.
   */
  std::vector<dof_id_type> * _n_oz;

  /**
   * Total number of degrees of freedom.
   */
  dof_id_type _n_dfs;

  /**
   * The total number of SCALAR dofs associated to
   * all SCALAR variables.
   */
  dof_id_type _n_SCALAR_dofs;

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Total number of degrees of freedom on old dof objects
   */
  dof_id_type _n_old_dfs;

  /**
   * First old DOF index on processor \p p.
   */
  std::vector<dof_id_type> _first_old_df;

  /**
   * Last old DOF index (plus 1) on processor \p p.
   */
  std::vector<dof_id_type> _end_old_df;

  /**
   * First old DOF index for SCALAR variable v, or garbage for
   * non-SCALAR variable v
   */
  std::vector<dof_id_type> _first_old_scalar_df;

#endif

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  /**
   * Data structure containing DOF constraints.  The ith
   * entry is the constraint matrix row for DOF i.
   */
  DofConstraints _dof_constraints, _stashed_dof_constraints;

  DofConstraintValueMap      _primal_constraint_values;

  AdjointDofConstraintValues _adjoint_constraint_values;
#endif

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * Data structure containing DofObject constraints.
   */
  NodeConstraints _node_constraints;
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS


#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * Data structure containing periodic boundaries.  The ith
   * entry is the constraint matrix row for boundaryid i.
   */
  UniquePtr<PeriodicBoundaries> _periodic_boundaries;
#endif

#ifdef LIBMESH_ENABLE_DIRICHLET
  /**
   * Check that all the ids in dirichlet_bcids are actually present in the mesh.
   * If not, this will throw an error.
   */
  void check_dirichlet_bcid_consistency (const MeshBase & mesh,
                                         const DirichletBoundary & boundary) const;
  /**
   * Data structure containing Dirichlet functions.  The ith
   * entry is the constraint matrix row for boundaryid i.
   */
  UniquePtr<DirichletBoundaries> _dirichlet_boundaries;

  /**
   * Data structure containing Dirichlet functions.  The ith
   * entry is the constraint matrix row for boundaryid i.
   */
  std::vector<DirichletBoundaries *> _adjoint_dirichlet_boundaries;
#endif

  friend class SparsityPattern::Build;

  /**
   * Bools to indicate if we override the --implicit_neighbor_dofs
   * commandline options.
   */
  bool _implicit_neighbor_dofs_initialized;
  bool _implicit_neighbor_dofs;
};


// ------------------------------------------------------------
// Dof Map inline member functions
inline
unsigned int DofMap::sys_number() const
{
  return _sys_number;
}



inline
const VariableGroup & DofMap::variable_group (const unsigned int g) const
{
  libmesh_assert_less (g, _variable_groups.size());

  return _variable_groups[g];
}



inline
const Variable & DofMap::variable (const unsigned int c) const
{
  libmesh_assert_less (c, _variables.size());

  return _variables[c];
}



inline
Order DofMap::variable_order (const unsigned int c) const
{
  libmesh_assert_less (c, _variables.size());

  return _variables[c].type().order;
}



inline
Order DofMap::variable_group_order (const unsigned int vg) const
{
  libmesh_assert_less (vg, _variable_groups.size());

  return _variable_groups[vg].type().order;
}



inline
const FEType & DofMap::variable_type (const unsigned int c) const
{
  libmesh_assert_less (c, _variables.size());

  return _variables[c].type();
}



inline
const FEType & DofMap::variable_group_type (const unsigned int vg) const
{
  libmesh_assert_less (vg, _variable_groups.size());

  return _variable_groups[vg].type();
}



inline
bool DofMap::is_constrained_node (const Node *
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

#ifdef LIBMESH_ENABLE_CONSTRAINTS

inline
bool DofMap::is_constrained_dof (const dof_id_type dof) const
{
  if (_dof_constraints.count(dof))
    return true;

  return false;
}


inline
bool DofMap::has_heterogenous_adjoint_constraints (const unsigned int qoi_num) const
{
  AdjointDofConstraintValues::const_iterator it =
    _adjoint_constraint_values.find(qoi_num);
  if (it == _adjoint_constraint_values.end())
    return false;
  if (it->second.empty())
    return false;

  return true;
}


inline
Number DofMap::has_heterogenous_adjoint_constraint (const unsigned int qoi_num,
                                                    const dof_id_type dof) const
{
  AdjointDofConstraintValues::const_iterator it =
    _adjoint_constraint_values.find(qoi_num);
  if (it != _adjoint_constraint_values.end())
    {
      DofConstraintValueMap::const_iterator rhsit =
        it->second.find(dof);
      if (rhsit == it->second.end())
        return 0;
      else
        return rhsit->second;
    }

  return 0;
}



inline
DofConstraintValueMap & DofMap::get_primal_constraint_values()
{
  return _primal_constraint_values;
}



#else

//--------------------------------------------------------------------
// Constraint-specific methods get inlined into nothing if
// constraints are disabled, so there's no reason for users not to
// use them.

inline void DofMap::constrain_element_matrix (DenseMatrix<Number> &,
                                              std::vector<dof_id_type> &,
                                              bool) const {}

inline void DofMap::constrain_element_matrix (DenseMatrix<Number> &,
                                              std::vector<dof_id_type> &,
                                              std::vector<dof_id_type> &,
                                              bool) const {}

inline void DofMap::constrain_element_vector (DenseVector<Number> &,
                                              std::vector<dof_id_type> &,
                                              bool) const {}

inline void DofMap::constrain_element_matrix_and_vector (DenseMatrix<Number> &,
                                                         DenseVector<Number> &,
                                                         std::vector<dof_id_type> &,
                                                         bool) const {}

inline void DofMap::constrain_element_dyad_matrix (DenseVector<Number> &,
                                                   DenseVector<Number> &,
                                                   std::vector<dof_id_type> &,
                                                   bool) const {}

inline void DofMap::enforce_constraints_exactly (const System &,
                                                 NumericVector<Number> *,
                                                 bool = false) const {}

inline void DofMap::enforce_adjoint_constraints_exactly (NumericVector<Number> &,
                                                         unsigned int) const {}

#endif // LIBMESH_ENABLE_CONSTRAINTS

} // namespace libMesh

#endif // LIBMESH_DOF_MAP_H
