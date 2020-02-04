// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SYSTEM_H
#define LIBMESH_SYSTEM_H

// Local Includes
#include "libmesh/elem_range.h"
#include "libmesh/enum_subset_solve_mode.h" // SUBSET_ZERO
#include "libmesh/enum_parallel_type.h" // PARALLEL
#include "libmesh/fem_function_base.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/parallel_object.h"
#include "libmesh/qoi_set.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/tensor_value.h" // For point_hessian
#include "libmesh/variable.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum FEMNormType : int;
}
#else
#include "libmesh/enum_norm_type.h"
#endif

// C++ includes
#include <cstddef>
#include <set>
#include <vector>
#include <memory>

// This define may be useful in --disable-optional builds when it is
// possible that libmesh will not have any solvers available.
#if defined(LIBMESH_HAVE_PETSC) || \
  defined(LIBMESH_HAVE_EIGEN)   || \
  defined(LIBMESH_HAVE_LASPACK) || \
  defined(LIBMESH_TRILINOS_HAVE_AZTECOO)
#define LIBMESH_HAVE_SOLVER 1
#endif

namespace libMesh
{

// Forward Declarations
class System;
class EquationSystems;
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;
class Xdr;
class DofMap;
template <typename Output> class FunctionBase;
class Parameters;
class ParameterVector;
template <typename> class PointTempl;
typedef PointTempl<Real> Point;
class SensitivityData;
template <typename T> class NumericVector;
template <typename T> class SparseMatrix;
template <typename T> class VectorValue;
typedef VectorValue<Number> NumberVectorValue;
typedef NumberVectorValue Gradient;
class SystemSubset;
class FEType;
class SystemNorm;

/**
 * \brief Manages consistently variables, degrees of freedom, and coefficient
 * vectors.
 *
 * This is the base class for classes which contain information related to any
 * physical process that might be simulated. Such information may range from
 * the actual solution values to algorithmic flags that may be used to control
 * the numerical methods employed. In general, use an EquationSystems
 * object to handle one or more of the children of this class.
 *
 * Especially, this class manages the variables of the differential equation,
 * the coefficient vectors and the \p DofMap, and ensures that these are
 * consistent. It provides storage for the solution.  Furthermore, (de-)
 * serialization functionality is provided.
 *
 * \author Benjamin S. Kirk
 * \date 2003-2004
 */
class System : public ReferenceCountedObject<System>,
               public ParallelObject
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  System (EquationSystems & es,
          const std::string & name,
          const unsigned int number);

  /**
   * Abstract base class to be used for system initialization.
   * A user class derived from this class may be used to
   * initialize the system values by attaching an object
   * with the method \p attach_init_object.
   */
  class Initialization
  {
  public:
    /**
     * Destructor.  Virtual because we will have virtual functions.
     */
    virtual ~Initialization () {}

    /**
     * Initialization function.  This function will be called
     * to initialize the system values upon creation and must
     * be provided by the user in a derived class.
     */
    virtual void initialize () = 0;
  };



  /**
   * Abstract base class to be used for system assembly.
   * A user class derived from this class may be used to
   * assemble the system by attaching an object
   * with the method \p attach_assemble_object.
   */
  class Assembly
  {
  public:
    /**
     * Destructor.  Virtual because we will have virtual functions.
     */
    virtual ~Assembly () {}

    /**
     * Assembly function.  This function will be called
     * to assemble the system prior to a solve and must
     * be provided by the user in a derived class.
     */
    virtual void assemble () = 0;
  };



  /**
   * Abstract base class to be used for system constraints.
   * A user class derived from this class may be used to
   * constrain the system by attaching an object
   * with the method \p attach_constraint_object.
   */
  class Constraint
  {
  public:
    /**
     * Destructor.  Virtual because we will have virtual functions.
     */
    virtual ~Constraint () {}

    /**
     * Constraint function.  This function will be called
     * to constrain the system prior to a solve and must
     * be provided by the user in a derived class.
     */
    virtual void constrain () = 0;
  };



  /**
   * Abstract base class to be used for quantities of interest.
   * A user class derived from this class may be used to
   * compute quantities of interest by attaching an object
   * with the method \p attach_QOI_object.
   */
  class QOI
  {
  public:
    /**
     * Destructor.  Virtual because we will have virtual functions.
     */
    virtual ~QOI () {}

    /**
     * Quantity of interest function.  This function will be called
     * to compute quantities of interest and must be provided by the
     * user in a derived class.
     */
    virtual void qoi (const QoISet & qoi_indices) = 0;
  };



  /**
   * Abstract base class to be used for derivatives of quantities
   * of interest. A user class derived from this class may be used
   * to compute quantities of interest by attaching an object
   * with the method \p attach_QOI_derivative_object.
   */
  class QOIDerivative
  {
  public:
    /**
     * Destructor.  Virtual because we will have virtual functions.
     */
    virtual ~QOIDerivative () {}

    /**
     * Quantity of interest derivative function. This function will
     * be called to compute derivatives of quantities of interest and
     * must be provided by the user in a derived class.
     */
    virtual void qoi_derivative (const QoISet & qoi_indices,
                                 bool include_liftfunc,
                                 bool apply_constraints) = 0;
  };



  /**
   * Destructor.
   */
  virtual ~System ();

  /**
   * The type of system.
   */
  typedef System sys_type;

  /**
   * \returns A reference to *this.
   */
  sys_type & system () { return *this; }

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * Initializes degrees of freedom on the current mesh.
   * Sets the
   */
  void init ();

  /**
   * Reinitializes degrees of freedom and other required data on the
   * current mesh.
   *
   * \note The matrix is not initialized at this time since it may not
   * be required for all applications. Should be overridden in derived
   * classes.
   */
  virtual void reinit ();

  /**
   * Reinitializes the constraints for this system.
   */
  virtual void reinit_constraints ();

  /**
   * \returns \p true iff this system has been initialized.
   */
  bool is_initialized();

  /**
   * Update the local values to reflect the solution
   * on neighboring processors.
   */
  virtual void update ();

  /**
   * Prepares \p matrix and \p _dof_map for matrix assembly.
   * Does not actually assemble anything.  For matrix assembly,
   * use the \p assemble() in derived classes.
   * Should be overridden in derived classes.
   */
  virtual void assemble ();

  /**
   * Calls user qoi function.
   * Can be overridden in derived classes.
   */
  virtual void assemble_qoi (const QoISet & qoi_indices = QoISet());

  /**
   * Calls user qoi derivative function.
   * Can be overridden in derived classes.
   */
  virtual void assemble_qoi_derivative (const QoISet & qoi_indices = QoISet(),
                                        bool include_liftfunc = true,
                                        bool apply_constraints = true);

  /**
   * Calls residual parameter derivative function.
   *
   * Library subclasses use finite differences by default.
   *
   * This should assemble the sensitivity rhs vectors to hold
   * -(partial R / partial p_i), making them ready to solve
   * the forward sensitivity equation.
   *
   * This method is only implemented in some derived classes.
   */
  virtual void assemble_residual_derivatives (const ParameterVector & parameters);

  /**
   * After calling this method, any solve will be restricted to the
   * given subdomain.  To disable this mode, call this method with \p
   * subset being a \p nullptr.
   */
  virtual void restrict_solve_to (const SystemSubset * subset,
                                  const SubsetSolveMode subset_solve_mode=SUBSET_ZERO);

  /**
   * Solves the system.  Should be overridden in derived systems.
   */
  virtual void solve () {}

  /**
   * Solves the sensitivity system, for the provided parameters.
   * Must be overridden in derived systems.
   *
   * \returns A pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   *
   * This method is only implemented in some derived classes.
   */
  virtual std::pair<unsigned int, Real>
  sensitivity_solve (const ParameterVector & parameters);

  /**
   * Assembles & solves the linear system(s) (dR/du)*u_w = sum(w_p*-dR/dp), for
   * those parameters p contained within \p parameters weighted by the
   * values w_p found within \p weights.
   *
   * \returns A pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   *
   * This method is only implemented in some derived classes.
   */
  virtual std::pair<unsigned int, Real>
  weighted_sensitivity_solve (const ParameterVector & parameters,
                              const ParameterVector & weights);

  /**
   * Solves the adjoint system, for the specified qoi indices, or for
   * every qoi if \p qoi_indices is nullptr.  Must be overridden in
   * derived systems.
   *
   * \returns A pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   *
   * This method is only implemented in some derived classes.
   */
  virtual std::pair<unsigned int, Real>
  adjoint_solve (const QoISet & qoi_indices = QoISet());

  /**
   * Assembles & solves the linear system(s)
   * (dR/du)^T*z_w = sum(w_p*(d^2q/dudp - d^2R/dudp*z)), for those
   * parameters p contained within \p parameters, weighted by the
   * values w_p found within \p weights.
   *
   * Assumes that adjoint_solve has already calculated z for each qoi
   * in \p qoi_indices.
   *
   * \returns A pair with the total number of linear iterations
   * performed and the (sum of the) final residual norms
   *
   * This method is only implemented in some derived classes.
   */
  virtual std::pair<unsigned int, Real>
  weighted_sensitivity_adjoint_solve (const ParameterVector & parameters,
                                      const ParameterVector & weights,
                                      const QoISet & qoi_indices = QoISet());
  /**
   * Accessor for the adjoint_already_solved boolean
   */
  bool is_adjoint_already_solved() const
  { return adjoint_already_solved;}

  /**
   * Setter for the adjoint_already_solved boolean
   */
  void set_adjoint_already_solved(bool setting)
  { adjoint_already_solved = setting;}


  /**
   * Solves for the derivative of each of the system's quantities of
   * interest q in \p qoi[qoi_indices] with respect to each parameter
   * in \p parameters, placing the result for qoi \p i and parameter
   * \p j into \p sensitivities[i][j].
   *
   * \note \p parameters is a const vector, not a vector-of-const;
   * parameter values in this vector need to be mutable for finite
   * differencing to work.
   *
   * Automatically chooses the forward method for problems with more
   * quantities of interest than parameters, or the adjoint method
   * otherwise.
   *
   * This method is only usable in derived classes which override
   * an implementation.
   */
  virtual void qoi_parameter_sensitivity (const QoISet & qoi_indices,
                                          const ParameterVector & parameters,
                                          SensitivityData & sensitivities);

  /**
   * Solves for parameter sensitivities using the adjoint method.
   *
   * This method is only implemented in some derived classes.
   */
  virtual void adjoint_qoi_parameter_sensitivity (const QoISet & qoi_indices,
                                                  const ParameterVector & parameters,
                                                  SensitivityData & sensitivities);

  /**
   * Solves for parameter sensitivities using the forward method.
   *
   * This method is only implemented in some derived classes.
   */
  virtual void forward_qoi_parameter_sensitivity (const QoISet & qoi_indices,
                                                  const ParameterVector & parameters,
                                                  SensitivityData & sensitivities);

  /**
   * For each of the system's quantities of interest q in
   * \p qoi[qoi_indices], and for a vector of parameters p, the
   * parameter sensitivity Hessian H_ij is defined as
   * H_ij = (d^2 q)/(d p_i d p_j)
   * This Hessian is the output of this method, where for each q_i,
   * H_jk is stored in \p hessian.second_derivative(i,j,k).
   *
   * This method is only implemented in some derived classes.
   */
  virtual void qoi_parameter_hessian(const QoISet & qoi_indices,
                                     const ParameterVector & parameters,
                                     SensitivityData & hessian);

  /**
   * For each of the system's quantities of interest q in
   * \p qoi[qoi_indices], and for a vector of parameters p, the
   * parameter sensitivity Hessian H_ij is defined as
   * H_ij = (d^2 q)/(d p_i d p_j)
   * The Hessian-vector product, for a vector v_k in parameter space, is
   * S_j = H_jk v_k
   * This product is the output of this method, where for each q_i,
   * S_j is stored in \p sensitivities[i][j].
   *
   * This method is only implemented in some derived classes.
   */
  virtual void qoi_parameter_hessian_vector_product(const QoISet & qoi_indices,
                                                    const ParameterVector & parameters,
                                                    const ParameterVector & vector,
                                                    SensitivityData & product);

  /**
   * \returns \p true when the other system contains
   * identical data, up to the given threshold.  Outputs
   * some diagnostic info when \p verbose is set.
   */
  virtual bool compare (const System & other_system,
                        const Real threshold,
                        const bool verbose) const;

  /**
   * \returns The system name.
   */
  const std::string & name () const;

  /**
   * \returns The type of system, helpful in identifying
   * which system type to use when reading equation system
   * data from file.  Should be overridden in derived classes.
   */
  virtual std::string system_type () const { return "Basic"; }

  /**
   * Projects arbitrary functions onto the current solution.
   * The function value \p f and its gradient \p g are
   * user-provided cloneable functors.
   * A gradient \p g is only required/used for projecting onto finite
   * element spaces with continuous derivatives.
   * If non-default \p Parameters are to be used, they can be provided
   * in the \p parameters argument.
   */
  void project_solution (FunctionBase<Number> * f,
                         FunctionBase<Gradient> * g = nullptr) const;

  /**
   * Projects arbitrary functions onto the current solution.
   * The function value \p f and its gradient \p g are
   * user-provided cloneable functors.
   * A gradient \p g is only required/used for projecting onto finite
   * element spaces with continuous derivatives.
   * If non-default \p Parameters are to be used, they can be provided
   * in the \p parameters argument.
   */
  void project_solution (FEMFunctionBase<Number> * f,
                         FEMFunctionBase<Gradient> * g = nullptr) const;

  /**
   * Projects arbitrary functions onto the current solution.
   * The function value \p fptr and its gradient \p gptr are
   * represented by function pointers.
   * A gradient \p gptr is only required/used for projecting onto
   * finite element spaces with continuous derivatives.
   */
  typedef Number (*ValueFunctionPointer)(const Point & p,
                                         const Parameters & Parameters,
                                         const std::string & sys_name,
                                         const std::string & unknown_name);
  typedef Gradient (*GradientFunctionPointer)(const Point & p,
                                              const Parameters & parameters,
                                              const std::string & sys_name,
                                              const std::string & unknown_name);
  void project_solution (ValueFunctionPointer fptr,
                         GradientFunctionPointer gptr,
                         const Parameters & parameters) const;

  /**
   * Projects arbitrary functions onto a vector of degree of freedom
   * values for the current system.
   * The function value \p f and its gradient \p g are
   * user-provided cloneable functors.
   * A gradient \p g is only required/used for projecting onto finite
   * element spaces with continuous derivatives.
   * If non-default \p Parameters are to be used, they can be provided
   * in the \p parameters argument.
   *
   * Constrain the new vector using the requested adjoint rather than
   * primal constraints if is_adjoint is non-negative.
   */
  void project_vector (NumericVector<Number> & new_vector,
                       FunctionBase<Number> * f,
                       FunctionBase<Gradient> * g = nullptr,
                       int is_adjoint = -1) const;

  /**
   * Projects arbitrary functions onto a vector of degree of freedom
   * values for the current system.
   * The function value \p f and its gradient \p g are
   * user-provided cloneable functors.
   * A gradient \p g is only required/used for projecting onto finite
   * element spaces with continuous derivatives.
   * If non-default \p Parameters are to be used, they can be provided
   * in the \p parameters argument.
   *
   * Constrain the new vector using the requested adjoint rather than
   * primal constraints if is_adjoint is non-negative.
   */
  void project_vector (NumericVector<Number> & new_vector,
                       FEMFunctionBase<Number> * f,
                       FEMFunctionBase<Gradient> * g = nullptr,
                       int is_adjoint = -1) const;

  /**
   * Projects arbitrary functions onto a vector of degree of freedom
   * values for the current system.
   * The function value \p fptr and its gradient \p gptr are
   * represented by function pointers.
   * A gradient \p gptr is only required/used for projecting onto
   * finite element spaces with continuous derivatives.
   *
   * Constrain the new vector using the requested adjoint rather than
   * primal constraints if is_adjoint is non-negative.
   */
  void project_vector (ValueFunctionPointer fptr,
                       GradientFunctionPointer gptr,
                       const Parameters & parameters,
                       NumericVector<Number> & new_vector,
                       int is_adjoint = -1) const;

  /**
   * Projects arbitrary boundary functions onto a vector of degree of
   * freedom values for the current system.
   * Only degrees of freedom which affect the function's trace on a
   * boundary in the set \p b are affected.
   * Only degrees of freedom associated with the variables listed in
   * the vector \p variables are projected.
   * The function value \p f and its gradient \p g are
   * user-provided cloneable functors.
   * A gradient \p g is only required/used for projecting onto finite
   * element spaces with continuous derivatives.
   * If non-default \p Parameters are to be used, they can be provided
   * in the \p parameters argument.
   */
  void boundary_project_solution (const std::set<boundary_id_type> & b,
                                  const std::vector<unsigned int> & variables,
                                  FunctionBase<Number> * f,
                                  FunctionBase<Gradient> * g = nullptr);

  /**
   * Projects arbitrary boundary functions onto a vector of degree of
   * freedom values for the current system.
   * Only degrees of freedom which affect the function's trace on a
   * boundary in the set \p b are affected.
   * Only degrees of freedom associated with the variables listed in
   * the vector \p variables are projected.
   * The function value \p fptr and its gradient \p gptr are
   * represented by function pointers.
   * A gradient \p gptr is only required/used for projecting onto
   * finite element spaces with continuous derivatives.
   */
  void boundary_project_solution (const std::set<boundary_id_type> & b,
                                  const std::vector<unsigned int> & variables,
                                  ValueFunctionPointer fptr,
                                  GradientFunctionPointer gptr,
                                  const Parameters & parameters);

  /**
   * Projects arbitrary boundary functions onto a vector of degree of
   * freedom values for the current system.
   * Only degrees of freedom which affect the function's trace on a
   * boundary in the set \p b are affected.
   * Only degrees of freedom associated with the variables listed in
   * the vector \p variables are projected.
   * The function value \p f and its gradient \p g are
   * user-provided cloneable functors.
   * A gradient \p g is only required/used for projecting onto finite
   * element spaces with continuous derivatives.
   * If non-default \p Parameters are to be used, they can be provided
   * in the \p parameters argument.
   *
   * Constrain the new vector using the requested adjoint rather than
   * primal constraints if is_adjoint is non-negative.
   */
  void boundary_project_vector (const std::set<boundary_id_type> & b,
                                const std::vector<unsigned int> & variables,
                                NumericVector<Number> & new_vector,
                                FunctionBase<Number> * f,
                                FunctionBase<Gradient> * g = nullptr,
                                int is_adjoint = -1) const;

  /**
   * Projects arbitrary boundary functions onto a vector of degree of
   * freedom values for the current system.
   * Only degrees of freedom which affect the function's trace on a
   * boundary in the set \p b are affected.
   * Only degrees of freedom associated with the variables listed in
   * the vector \p variables are projected.
   * The function value \p fptr and its gradient \p gptr are
   * represented by function pointers.
   * A gradient \p gptr is only required/used for projecting onto
   * finite element spaces with continuous derivatives.
   *
   * Constrain the new vector using the requested adjoint rather than
   * primal constraints if is_adjoint is non-negative.
   */
  void boundary_project_vector (const std::set<boundary_id_type> & b,
                                const std::vector<unsigned int> & variables,
                                ValueFunctionPointer fptr,
                                GradientFunctionPointer gptr,
                                const Parameters & parameters,
                                NumericVector<Number> & new_vector,
                                int is_adjoint = -1) const;

  /**
   * \returns The system number.
   */
  unsigned int number () const;

  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on all processors.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<Number> & global_soln) const;

  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on processor \p dest_proc.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<Number> & global_soln,
                               const processor_id_type dest_proc) const;

  /**
   * \returns A constant reference to this systems's \p _mesh.
   */
  const MeshBase & get_mesh() const;

  /**
   * \returns A reference to this systems's \p _mesh.
   */
  MeshBase & get_mesh();

  /**
   * \returns A constant reference to this system's \p _dof_map.
   */
  const DofMap & get_dof_map() const;

  /**
   * \returns A writable reference to this system's \p _dof_map.
   */
  DofMap & get_dof_map();

  /**
   * \returns A constant reference to this system's parent EquationSystems object.
   */
  const EquationSystems & get_equation_systems() const { return _equation_systems; }

  /**
   * \returns A reference to this system's parent EquationSystems object.
   */
  EquationSystems & get_equation_systems() { return _equation_systems; }

  /**
   * \returns \p true if the system is active, \p false otherwise.
   * An active system will be solved.
   */
  bool active () const;

  /**
   * Activates the system.  Only active systems are solved.
   */
  void activate ();

  /**
   * Deactivates the system.  Only active systems are solved.
   */
  void deactivate ();

  /**
   * Sets the system to be "basic only": i.e. advanced system
   * components such as ImplicitSystem matrices may not be
   * initialized.  This is useful for efficiency in certain utility
   * programs that never use System::solve().  This method must be
   * called after the System or derived class is created but before it
   * is initialized; e.g. from within EquationSystems::read()
   */
  void set_basic_system_only ();

  /**
   * Vector iterator typedefs.
   */
  typedef std::map<std::string, NumericVector<Number> *>::iterator       vectors_iterator;
  typedef std::map<std::string, NumericVector<Number> *>::const_iterator const_vectors_iterator;

  /**
   * Beginning of vectors container
   */
  vectors_iterator vectors_begin ();

  /**
   * Beginning of vectors container
   */
  const_vectors_iterator vectors_begin () const;

  /**
   * End of vectors container
   */
  vectors_iterator vectors_end ();

  /**
   * End of vectors container
   */
  const_vectors_iterator vectors_end () const;

  /**
   * Adds the additional vector \p vec_name to this system.  All the
   * additional vectors are similarly distributed, like the \p
   * solution, and initialized to zero.
   *
   * By default vectors added by add_vector are projected to changed grids by
   * reinit().  To zero them instead (more efficient), pass "false" as the
   * second argument
   */
  NumericVector<Number> & add_vector (const std::string & vec_name,
                                      const bool projections=true,
                                      const ParallelType type = PARALLEL);

  /**
   * Removes the additional vector \p vec_name from this system
   */
  void remove_vector(const std::string & vec_name);

  /**
   * Tells the System whether or not to project the solution vector onto new
   * grids when the system is reinitialized.  The solution will be projected
   * unless project_solution_on_reinit() = false is called.
   */
  bool & project_solution_on_reinit (void)
  { return _solution_projection; }

  /**
   * \returns \p true if this \p System has a vector associated with the
   * given name, \p false otherwise.
   */
  bool have_vector (const std::string & vec_name) const;

  /**
   * \returns A const pointer to the vector if this \p System has a
   * vector associated with the given name, \p nullptr otherwise.
   */
  const NumericVector<Number> * request_vector (const std::string & vec_name) const;

  /**
   * \returns A pointer to the vector if this \p System has a
   * vector associated with the given name, \p nullptr otherwise.
   */
  NumericVector<Number> * request_vector (const std::string & vec_name);

  /**
   * \returns A const pointer to this system's additional vector
   * number \p vec_num (where the vectors are counted starting with
   * 0), or \p nullptr if the system has no such vector.
   */
  const NumericVector<Number> * request_vector (const unsigned int vec_num) const;

  /**
   * \returns A writable pointer to this system's additional
   * vector number \p vec_num (where the vectors are counted starting
   * with 0), or \p nullptr if the system has no such vector.
   */
  NumericVector<Number> * request_vector (const unsigned int vec_num);

  /**
   * \returns A const reference to this system's additional vector
   * named \p vec_name.  Access is only granted when the vector is already
   * properly initialized.
   */
  const NumericVector<Number> & get_vector (const std::string & vec_name) const;

  /**
   * \returns A writable reference to this system's additional vector
   * named \p vec_name.  Access is only granted when the vector is already
   * properly initialized.
   */
  NumericVector<Number> & get_vector (const std::string & vec_name);

  /**
   * \returns A const reference to this system's additional vector
   * number \p vec_num (where the vectors are counted starting with
   * 0).
   */
  const NumericVector<Number> & get_vector (const unsigned int vec_num) const;

  /**
   * \returns A writable reference to this system's additional
   * vector number \p vec_num (where the vectors are counted starting
   * with 0).
   */
  NumericVector<Number> & get_vector (const unsigned int vec_num);

  /**
   * \returns The name of this system's additional vector number \p
   * vec_num (where the vectors are counted starting with 0).
   */
  const std::string & vector_name (const unsigned int vec_num) const;

  /**
   * \returns The name of a system vector, given a reference to that vector
   */
  const std::string & vector_name (const NumericVector<Number> & vec_reference) const;

  /**
   * Allows one to set the QoI index controlling whether the vector
   * identified by vec_name represents a solution from the adjoint
   * (qoi_num >= 0) or primal (qoi_num == -1) space.  This becomes
   * significant if those spaces have differing heterogeneous
   * Dirichlet constraints.
   *
   * qoi_num == -2 can be used to indicate a vector which should not
   * be affected by constraints during projection operations.
   */
  void set_vector_as_adjoint (const std::string & vec_name, int qoi_num);

  /**
   * \returns The integer describing whether the vector identified by
   * vec_name represents a solution from an adjoint (non-negative) or
   * the primal (-1) space.
   */
  int vector_is_adjoint (const std::string & vec_name) const;

  /**
   * Allows one to set the boolean controlling whether the vector
   * identified by vec_name should be "preserved": projected to new
   * meshes, saved, etc.
   */
  void set_vector_preservation (const std::string & vec_name, bool preserve);

  /**
   * \returns The boolean describing whether the vector identified by
   * vec_name should be "preserved": projected to new meshes, saved,
   * etc.
   */
  bool vector_preservation (const std::string & vec_name) const;

  /**
   * \returns A reference to one of the system's adjoint solution
   * vectors, by default the one corresponding to the first qoi.
   * Creates the vector if it doesn't already exist.
   */
  NumericVector<Number> & add_adjoint_solution(unsigned int i=0);

  /**
   * \returns A reference to one of the system's adjoint solution
   * vectors, by default the one corresponding to the first qoi.
   */
  NumericVector<Number> & get_adjoint_solution(unsigned int i=0);

  /**
   * \returns A reference to one of the system's adjoint solution
   * vectors, by default the one corresponding to the first qoi.
   */
  const NumericVector<Number> & get_adjoint_solution(unsigned int i=0) const;

  /**
   * \returns A reference to one of the system's solution sensitivity
   * vectors, by default the one corresponding to the first parameter.
   * Creates the vector if it doesn't already exist.
   */
  NumericVector<Number> & add_sensitivity_solution(unsigned int i=0);

  /**
   * \returns A reference to one of the system's solution sensitivity
   * vectors, by default the one corresponding to the first parameter.
   */
  NumericVector<Number> & get_sensitivity_solution(unsigned int i=0);

  /**
   * \returns A reference to one of the system's solution sensitivity
   * vectors, by default the one corresponding to the first parameter.
   */
  const NumericVector<Number> & get_sensitivity_solution(unsigned int i=0) const;

  /**
   * \returns A reference to one of the system's weighted sensitivity
   * adjoint solution vectors, by default the one corresponding to the
   * first qoi.
   * Creates the vector if it doesn't already exist.
   */
  NumericVector<Number> & add_weighted_sensitivity_adjoint_solution(unsigned int i=0);

  /**
   * \returns A reference to one of the system's weighted sensitivity
   * adjoint solution vectors, by default the one corresponding to the
   * first qoi.
   */
  NumericVector<Number> & get_weighted_sensitivity_adjoint_solution(unsigned int i=0);

  /**
   * \returns A reference to one of the system's weighted sensitivity
   * adjoint solution vectors, by default the one corresponding to the
   * first qoi.
   */
  const NumericVector<Number> & get_weighted_sensitivity_adjoint_solution(unsigned int i=0) const;

  /**
   * \returns A reference to the solution of the last weighted
   * sensitivity solve
   * Creates the vector if it doesn't already exist.
   */
  NumericVector<Number> & add_weighted_sensitivity_solution();

  /**
   * \returns A reference to the solution of the last weighted
   * sensitivity solve
   */
  NumericVector<Number> & get_weighted_sensitivity_solution();

  /**
   * \returns A reference to the solution of the last weighted
   * sensitivity solve
   */
  const NumericVector<Number> & get_weighted_sensitivity_solution() const;

  /**
   * \returns A reference to one of the system's adjoint rhs
   * vectors, by default the one corresponding to the first qoi.
   * Creates the vector if it doesn't already exist.
   */
  NumericVector<Number> & add_adjoint_rhs(unsigned int i=0);

  /**
   * \returns A reference to one of the system's adjoint rhs
   * vectors, by default the one corresponding to the first qoi.
   * This what the user's QoI derivative code should assemble
   * when setting up an adjoint problem
   */
  NumericVector<Number> & get_adjoint_rhs(unsigned int i=0);

  /**
   * \returns A reference to one of the system's adjoint rhs
   * vectors, by default the one corresponding to the first qoi.
   */
  const NumericVector<Number> & get_adjoint_rhs(unsigned int i=0) const;

  /**
   * \returns A reference to one of the system's sensitivity rhs
   * vectors, by default the one corresponding to the first parameter.
   * Creates the vector if it doesn't already exist.
   */
  NumericVector<Number> & add_sensitivity_rhs(unsigned int i=0);

  /**
   * \returns A reference to one of the system's sensitivity rhs
   * vectors, by default the one corresponding to the first parameter.
   * By default these vectors are built by the library, using finite
   * differences, when \p assemble_residual_derivatives() is called.
   *
   * When assembled, this vector should hold
   * -(partial R / partial p_i)
   */
  NumericVector<Number> & get_sensitivity_rhs(unsigned int i=0);

  /**
   * \returns A reference to one of the system's sensitivity rhs
   * vectors, by default the one corresponding to the first parameter.
   */
  const NumericVector<Number> & get_sensitivity_rhs(unsigned int i=0) const;

  /**
   * \returns The number of vectors (in addition to the solution)
   * handled by this system
   * This is the size of the \p _vectors map
   */
  unsigned int n_vectors () const;

  /**
   * \returns The number of matrices
   * handled by this system.
   *
   * This will return 0 by default but can be overridden.
   */
  virtual unsigned int n_matrices () const;

  /**
   * \returns The number of variables in the system
   */
  unsigned int n_vars() const;

  /**
   * \returns The number of \p VariableGroup variable groups in the system
   */
  unsigned int n_variable_groups() const;

  /**
   * \returns The total number of scalar components in the system's
   * variables.  This will equal \p n_vars() in the case of all
   * scalar-valued variables.
   */
  unsigned int n_components() const;

  /**
   * \returns The number of degrees of freedom in the system
   */
  dof_id_type n_dofs() const;

  /**
   * \returns The number of active degrees of freedom
   * for this System.
   */
  dof_id_type n_active_dofs() const;

  /**
   * \returns The total number of constrained degrees of freedom
   * in the system.
   */
  dof_id_type n_constrained_dofs() const;

  /**
   * \returns The number of constrained degrees of freedom
   * on this processor.
   */
  dof_id_type n_local_constrained_dofs() const;

  /**
   * \returns The number of degrees of freedom local
   * to this processor
   */
  dof_id_type n_local_dofs() const;

  /**
   * Adds the variable \p var to the list of variables
   * for this system. If \p active_subdomains is either \p nullptr
   * (the default) or points to an empty set, then it will be assumed that
   * \p var has no subdomain restrictions
   *
   * \returns The index number for the new variable.
   */
  unsigned int add_variable (const std::string & var,
                             const FEType & type,
                             const std::set<subdomain_id_type> * const active_subdomains = nullptr);

  /**
   * Adds the variable \p var to the list of variables
   * for this system.  Same as before, but assumes \p LAGRANGE
   * as default value for \p FEType.family. If \p active_subdomains is either
   * \p nullptr (the default) or points to an empty set, then it will be assumed
   * that \p var has no subdomain restrictions
   */
  unsigned int add_variable (const std::string & var,
                             const Order order = FIRST,
                             const FEFamily = LAGRANGE,
                             const std::set<subdomain_id_type> * const active_subdomains = nullptr);

  /**
   * Adds the variable \p var to the list of variables
   * for this system. If \p active_subdomains is either \p nullptr
   * (the default) or points to an empty set, then it will be assumed that
   * \p var has no subdomain restrictions
   *
   * \returns The index number for the new variable.
   */
  unsigned int add_variables (const std::vector<std::string> & vars,
                              const FEType & type,
                              const std::set<subdomain_id_type> * const active_subdomains = nullptr);

  /**
   * Adds the variable \p var to the list of variables
   * for this system.  Same as before, but assumes \p LAGRANGE
   * as default value for \p FEType.family. If \p active_subdomains is either
   * \p nullptr (the default) or points to an empty set, then it will be assumed that
   * \p var has no subdomain restrictions
   */
  unsigned int add_variables (const std::vector<std::string> & vars,
                              const Order order = FIRST,
                              const FEFamily = LAGRANGE,
                              const std::set<subdomain_id_type> * const active_subdomains = nullptr);

  /**
   * Return a constant reference to \p Variable \p var.
   */
  const Variable & variable (unsigned int var) const;

  /**
   * Return a constant reference to \p VariableGroup \p vg.
   */
  const VariableGroup & variable_group (unsigned int vg) const;

  /**
   * \returns \p true if a variable named \p var exists in this System
   */
  bool has_variable(const std::string & var) const;

  /**
   * \returns The name of variable \p i.
   */
  const std::string & variable_name(const unsigned int i) const;

  /**
   * \returns The variable number associated with
   * the user-specified variable named \p var.
   */
  unsigned short int variable_number (const std::string & var) const;

  /**
   * Fills \p all_variable_numbers with all the variable numbers for the
   * variables that have been added to this system.
   */
  void get_all_variable_numbers(std::vector<unsigned int> & all_variable_numbers) const;

  /**
   * \returns An index, starting from 0 for the first component of the
   * first variable, and incrementing for each component of each
   * (potentially vector-valued) variable in the system in order.
   * For systems with only scalar-valued variables, this will be the
   * same as \p variable_number(var)
   *
   * Irony: currently our only non-scalar-valued variable type is
   * SCALAR.
   */
  unsigned int variable_scalar_number (const std::string & var,
                                       unsigned int component) const;

  /**
   * \returns An index, starting from 0 for the first component of the
   * first variable, and incrementing for each component of each
   * (potentially vector-valued) variable in the system in order.
   * For systems with only scalar-valued variables, this will be the
   * same as \p var_num
   *
   * Irony: currently our only non-scalar-valued variable type is
   * SCALAR.
   */
  unsigned int variable_scalar_number (unsigned int var_num,
                                       unsigned int component) const;


  /**
   * \returns The finite element type variable number \p i.
   */
  const FEType & variable_type (const unsigned int i) const;

  /**
   * \returns The finite element type for variable \p var.
   */
  const FEType & variable_type (const std::string & var) const;

  /**
   * \returns \p true when \p VariableGroup structures should be
   * automatically identified, \p false otherwise.
   */
  bool identify_variable_groups () const;

  /**
   * Toggle automatic \p VariableGroup identification.
   */
  void identify_variable_groups (const bool);

  /**
   * \returns A norm of variable \p var in the vector \p v, in the specified
   * norm (e.g. L2, L_INF, H1)
   */
  Real calculate_norm(const NumericVector<Number> & v,
                      unsigned int var,
                      FEMNormType norm_type,
                      std::set<unsigned int> * skip_dimensions=nullptr) const;

  /**
   * \returns A norm of the vector \p v, using \p component_norm and \p
   * component_scale to choose and weight the norms of each variable.
   */
  Real calculate_norm(const NumericVector<Number> & v,
                      const SystemNorm & norm,
                      std::set<unsigned int> * skip_dimensions=nullptr) const;

  /**
   * Reads the basic data header for this System.
   */
  void read_header (Xdr & io,
                    const std::string & version,
                    const bool read_header=true,
                    const bool read_additional_data=true,
                    const bool read_legacy_format=false);

  /**
   * Reads additional data, namely vectors, for this System.
   *
   * \deprecated The ability to read XDR data files in the old (aka
   * "legacy") XDR format has been deprecated for many years, this
   * capability may soon disappear altogether.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  void read_legacy_data (Xdr & io,
                         const bool read_additional_data=true);
#endif

  /**
   * Reads additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   */
  template <typename ValType>
  void read_serialized_data (Xdr & io,
                             const bool read_additional_data=true);
  /**
   * Non-templated version for backward compatibility.
   *
   * Reads additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   */
  void read_serialized_data (Xdr & io,
                             const bool read_additional_data=true)
  { read_serialized_data<Number>(io, read_additional_data); }

  /**
   * Read a number of identically distributed vectors.  This method
   * allows for optimization for the multiple vector case by only communicating
   * the metadata once.
   */
  template <typename InValType>
  std::size_t read_serialized_vectors (Xdr & io,
                                       const std::vector<NumericVector<Number> *> & vectors) const;

  /**
   * Non-templated version for backward compatibility.
   *
   * Read a number of identically distributed vectors.  This method
   * allows for optimization for the multiple vector case by only communicating
   * the metadata once.
   */
  std::size_t read_serialized_vectors (Xdr & io,
                                       const std::vector<NumericVector<Number> *> & vectors) const
  { return read_serialized_vectors<Number>(io, vectors); }

  /**
   * Reads additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   * This method will read an individual file for each processor in the simulation
   * where the local solution components for that processor are stored.
   */
  template <typename InValType>
  void read_parallel_data (Xdr & io,
                           const bool read_additional_data);

  /**
   * Non-templated version for backward compatibility.
   *
   * Reads additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   * This method will read an individual file for each processor in the simulation
   * where the local solution components for that processor are stored.
   */
  void read_parallel_data (Xdr & io,
                           const bool read_additional_data)
  { read_parallel_data<Number>(io, read_additional_data); }

  /**
   * Writes the basic data header for this System.
   */
  void write_header (Xdr & io,
                     const std::string & version,
                     const bool write_additional_data) const;

  /**
   * Writes additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   */
  void write_serialized_data (Xdr & io,
                              const bool write_additional_data = true) const;

  /**
   * Serialize & write a number of identically distributed vectors.  This method
   * allows for optimization for the multiple vector case by only communicating
   * the metadata once.
   */
  std::size_t write_serialized_vectors (Xdr & io,
                                        const std::vector<const NumericVector<Number> *> & vectors) const;

  /**
   * Writes additional data, namely vectors, for this System.
   * This method may safely be called on a distributed-memory mesh.
   * This method will create an individual file for each processor in the simulation
   * where the local solution components for that processor will be stored.
   */
  void write_parallel_data (Xdr & io,
                            const bool write_additional_data) const;

  /**
   * \returns A string containing information about the
   * system.
   */
  std::string get_info () const;

  /**
   * Register a user function to use in initializing the system.
   */
  void attach_init_function (void fptr(EquationSystems & es,
                                       const std::string & name));

  /**
   * Register a user class to use to initialize the system.
   *
   * \note This is exclusive with the \p attach_init_function.
   */
  void attach_init_object (Initialization & init);

  /**
   * Register a user function to use in assembling the system
   * matrix and RHS.
   */
  void attach_assemble_function (void fptr(EquationSystems & es,
                                           const std::string & name));

  /**
   * Register a user object to use in assembling the system
   * matrix and RHS.
   */
  void attach_assemble_object (Assembly & assemble);

  /**
   * Register a user function for imposing constraints.
   */
  void attach_constraint_function (void fptr(EquationSystems & es,
                                             const std::string & name));

  /**
   * Register a user object for imposing constraints.
   */
  void attach_constraint_object (Constraint & constrain);

  /**
   * Register a user function for evaluating the quantities of interest,
   * whose values should be placed in \p System::qoi
   */
  void attach_QOI_function (void fptr(EquationSystems & es,
                                      const std::string & name,
                                      const QoISet & qoi_indices));

  /**
   * Register a user object for evaluating the quantities of interest,
   * whose values should be placed in \p System::qoi
   */
  void attach_QOI_object (QOI & qoi);

  /**
   * Register a user function for evaluating derivatives of a quantity
   * of interest with respect to test functions, whose values should
   * be placed in \p System::rhs
   */
  void attach_QOI_derivative (void fptr(EquationSystems & es,
                                        const std::string & name,
                                        const QoISet & qoi_indices,
                                        bool include_liftfunc,
                                        bool apply_constraints));

  /**
   * Register a user object for evaluating derivatives of a quantity
   * of interest with respect to test functions, whose values should
   * be placed in \p System::rhs
   */
  void attach_QOI_derivative_object (QOIDerivative & qoi_derivative);

  /**
   * Calls user's attached initialization function, or is overridden by
   * the user in derived classes.
   */
  virtual void user_initialization ();

  /**
   * Calls user's attached assembly function, or is overridden by
   * the user in derived classes.
   */
  virtual void user_assembly ();

  /**
   * Calls user's attached constraint function, or is overridden by
   * the user in derived classes.
   */
  virtual void user_constrain ();

  /**
   * Calls user's attached quantity of interest function, or is
   * overridden by the user in derived classes.
   */
  virtual void user_QOI (const QoISet & qoi_indices);

  /**
   * Calls user's attached quantity of interest derivative function,
   * or is overridden by the user in derived classes.
   */
  virtual void user_QOI_derivative (const QoISet & qoi_indices = QoISet(),
                                    bool include_liftfunc = true,
                                    bool apply_constraints = true);

  /**
   * Re-update the local values when the mesh has changed.
   * This method takes the data updated by \p update() and
   * makes it up-to-date on the current mesh.
   */
  virtual void re_update ();

  /**
   * Restrict vectors after the mesh has coarsened
   */
  virtual void restrict_vectors ();

  /**
   * Prolong vectors after the mesh has refined
   */
  virtual void prolong_vectors ();

  /**
   * Flag which tells the system to whether or not to
   * call the user assembly function during each call to solve().
   * By default, every call to solve() begins with a call to the
   * user assemble, so this flag is true.  (For explicit systems,
   * "solving" the system occurs during the assembly step, so this
   * flag is always true for explicit systems.)
   *
   * You will only want to set this to false if you need direct
   * control over when the system is assembled, and are willing to
   * track the state of its assembly yourself.  An example of such a
   * case is an implicit system with multiple right hand sides.  In
   * this instance, a single assembly would likely be followed with
   * multiple calls to solve.
   *
   * The frequency system and Newmark system have their own versions
   * of this flag, called _finished_assemble, which might be able to
   * be replaced with this more general concept.
   */
  bool assemble_before_solve;

  /**
   * Avoids use of any cached data that might affect any solve result.  Should
   * be overridden in derived systems.
   */
  virtual void disable_cache ();

  /**
   * A boolean to be set to true by systems using elem_fixed_solution,
   * for optional use by e.g. stabilized methods.
   * False by default.
   *
   * \note For FEMSystem users, if this variable is set to true, it
   * must be before init_data() is called.
   */
  bool use_fixed_solution;

  /**
   * A member int that can be employed to indicate increased or
   * reduced quadrature order.
   *
   * \note For FEMSystem users, by default, when calling the
   * user-defined residual functions, the FEMSystem will first set up
   * an appropriate FEType::default_quadrature_rule() object for
   * performing the integration.  This rule will integrate elements of
   * order up to 2*p+1 exactly (where p is the sum of the base FEType
   * and local p refinement levels), but if additional (or reduced)
   * quadrature accuracy is desired then this extra_quadrature_order
   * (default 0) will be added.
   */
  int extra_quadrature_order;


  //--------------------------------------------------
  // The solution and solution access members

  /**
   * \returns The current solution for the specified global
   * DOF.
   */
  Number current_solution (const dof_id_type global_dof_number) const;

  /**
   * Data structure to hold solution values.
   */
  std::unique_ptr<NumericVector<Number>> solution;

  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.  This vector is necessarily larger
   * than the \p solution vector in the case of a parallel
   * simulation.  The \p update() member is used to synchronize
   * the contents of the \p solution and \p current_local_solution
   * vectors.
   */
  std::unique_ptr<NumericVector<Number>> current_local_solution;

  /**
   * For time-dependent problems, this is the time t at the beginning of
   * the current timestep.
   *
   * \note For DifferentiableSystem users: do \e not access this time
   * during an assembly!  Use the DiffContext::time value instead to
   * get correct results.
   */
  Real time;

  /**
   * Number of currently active quantities of interest.
   */
  unsigned int n_qois() const;

  /**
   * Values of the quantities of interest.  This vector needs
   * to be both resized and filled by the user before any quantity of
   * interest assembly is done and before any sensitivities are
   * calculated.
   */
  std::vector<Number> qoi;

  /**
   * \returns The value of the solution variable \p var at the physical
   * point \p p in the mesh, without knowing a priori which element
   * contains \p p, using the degree of freedom coefficients in \p sol
   * (or in \p current_local_solution if \p sol is left null).
   *
   * \note This function uses \p MeshBase::sub_point_locator(); users
   * may or may not want to call \p MeshBase::clear_point_locator()
   * afterward.  Also, point_locator() is expensive (N log N for
   * initial construction, log N for evaluations).  Avoid using this
   * function in any context where you are already looping over
   * elements.
   *
   * Because the element containing \p p may lie on any processor,
   * this function is parallel-only.
   *
   * By default this method expects the point to reside inside the domain
   * and will abort if no element can be found which contains \p p.  The
   * optional parameter \p insist_on_success can be set to false to allow
   * the method to return 0 when the point is not located.
   */
  Number point_value(unsigned int var, const Point & p,
                     const bool insist_on_success = true,
                     const NumericVector<Number> * sol = nullptr) const;

  /**
   * \returns The value of the solution variable \p var at the physical
   * point \p p contained in local Elem \p e, using the degree of
   * freedom coefficients in \p sol (or in \p current_local_solution
   * if \p sol is left null).
   *
   * This version of point_value can be run in serial, but assumes \p e is in
   * the local mesh partition or is algebraically ghosted.
   */
  Number point_value(unsigned int var, const Point & p, const Elem & e,
                     const NumericVector<Number> * sol = nullptr) const;

  /**
   * Calls the version of point_value() which takes a reference.
   * This function exists only to prevent people from calling the
   * version of point_value() that has a boolean third argument, which
   * would result in unnecessary PointLocator calls.
   */
  Number point_value(unsigned int var, const Point & p, const Elem * e) const;

  /**
   * Calls the parallel version of point_value().
   * This function exists only to prevent people from accidentally
   * calling the version of point_value() that has a boolean third
   * argument, which would result in incorrect output.
   */
  Number point_value(unsigned int var, const Point & p, const NumericVector<Number> * sol) const;

  /**
   * \returns The gradient of the solution variable \p var at the physical
   * point \p p in the mesh, similarly to point_value.
   */
  Gradient point_gradient(unsigned int var, const Point & p,
                          const bool insist_on_success = true,
                          const NumericVector<Number> * sol = nullptr) const;

  /**
   * \returns The gradient of the solution variable \p var at the physical
   * point \p p in local Elem \p e in the mesh, similarly to point_value.
   */
  Gradient point_gradient(unsigned int var, const Point & p, const Elem & e,
                          const NumericVector<Number> * sol = nullptr) const;

  /**
   * Calls the version of point_gradient() which takes a reference.
   * This function exists only to prevent people from calling the
   * version of point_gradient() that has a boolean third argument, which
   * would result in unnecessary PointLocator calls.
   */
  Gradient point_gradient(unsigned int var, const Point & p, const Elem * e) const;

  /**
   * Calls the parallel version of point_gradient().
   * This function exists only to prevent people from accidentally
   * calling the version of point_gradient() that has a boolean third
   * argument, which would result in incorrect output.
   */
  Gradient point_gradient(unsigned int var, const Point & p, const NumericVector<Number> * sol) const;

  /**
   * \returns The second derivative tensor of the solution variable \p var
   * at the physical point \p p in the mesh, similarly to point_value.
   */
  Tensor point_hessian(unsigned int var, const Point & p,
                       const bool insist_on_success = true,
                       const NumericVector<Number> * sol = nullptr) const;

  /**
   * \returns The second derivative tensor of the solution variable \p var
   * at the physical point \p p in local Elem \p e in the mesh, similarly to
   * point_value.
   */
  Tensor point_hessian(unsigned int var, const Point & p, const Elem & e,
                       const NumericVector<Number> * sol = nullptr) const;

  /**
   * Calls the version of point_hessian() which takes a reference.
   * This function exists only to prevent people from calling the
   * version of point_hessian() that has a boolean third argument, which
   * would result in unnecessary PointLocator calls.
   */
  Tensor point_hessian(unsigned int var, const Point & p, const Elem * e) const;

  /**
   * Calls the parallel version of point_hessian().
   * This function exists only to prevent people from accidentally
   * calling the version of point_hessian() that has a boolean third
   * argument, which would result in incorrect output.
   */
  Tensor point_hessian(unsigned int var, const Point & p, const NumericVector<Number> * sol) const;


  /**
   * Fills the std::set with the degrees of freedom on the local
   * processor corresponding the the variable number passed in.
   */
  void local_dof_indices (const unsigned int var,
                          std::set<dof_id_type> & var_indices) const;

  /**
   * Zeroes all dofs in \p v that correspond to variable number \p
   * var_num.
   */
  void zero_variable (NumericVector<Number> & v, unsigned int var_num) const;


  /**
   * \returns A writable reference to a boolean that determines if this system
   * can be written to file or not.  If set to \p true, then
   * \p EquationSystems::write will ignore this system.
   */
  bool & hide_output() { return _hide_output; }

#ifdef LIBMESH_HAVE_METAPHYSICL
  /**
   * This method creates a projection matrix which corresponds to the
   * operation of project_vector between old and new solution spaces.
   *
   * Heterogeneous Dirichlet boundary conditions are *not* taken into
   * account here; if this matrix is used for prolongation (mesh
   * refinement) on a side with a heterogeneous BC, the newly created
   * degrees of freedom on that side will still match the coarse grid
   * approximation of the BC, not the fine grid approximation.
   */
  void projection_matrix (SparseMatrix<Number> & proj_mat) const;
#endif // LIBMESH_HAVE_METAPHYSICL

protected:

  /**
   * Initializes the data for the system.
   *
   * \note This is called before any user-supplied initialization
   * function so that all required storage will be available.
   */
  virtual void init_data ();

  /**
   * Projects the vector defined on the old mesh onto the
   * new mesh.
   *
   * Constrain the new vector using the requested adjoint rather than
   * primal constraints if is_adjoint is non-negative.
   */
  void project_vector (NumericVector<Number> &,
                       int is_adjoint = -1) const;

  /**
   * Projects the vector defined on the old mesh onto the
   * new mesh. The original vector is unchanged and the new vector
   * is passed through the second argument.
   *
   * Constrain the new vector using the requested adjoint rather than
   * primal constraints if is_adjoint is non-negative.
   */
  void project_vector (const NumericVector<Number> &,
                       NumericVector<Number> &,
                       int is_adjoint = -1) const;

private:
  /**
   * This isn't a copyable object, so let's make sure nobody tries.
   *
   * We won't even bother implementing this; we'll just make sure that
   * the compiler doesn't implement a default.
   */
  System (const System &);

  /**
   * This isn't a copyable object, so let's make sure nobody tries.
   *
   * We won't even bother implementing this; we'll just make sure that
   * the compiler doesn't implement a default.
   */
  System & operator=(const System &);

  /**
   * Finds the discrete norm for the entries in the vector
   * corresponding to Dofs associated with var.
   */
  Real discrete_var_norm (const NumericVector<Number> & v,
                          unsigned int var,
                          FEMNormType norm_type) const;

  /**
   * Reads an input vector from the stream \p io and assigns
   * the values to a set of \p DofObjects.  This method uses
   * blocked input and is safe to call on a distributed memory-mesh.
   * Unless otherwise specified, all variables are read.
   *
   * If an entry in \p vecs is a null pointer, the corresponding data
   * is read (incrementing the file read location) but discarded.
   */
  template <typename iterator_type, typename InValType>
  std::size_t read_serialized_blocked_dof_objects (const dof_id_type n_objects,
                                                   const iterator_type begin,
                                                   const iterator_type end,
                                                   const InValType dummy,
                                                   Xdr & io,
                                                   const std::vector<NumericVector<Number> *> & vecs,
                                                   const unsigned int var_to_read=libMesh::invalid_uint) const;

  /**
   * Reads the SCALAR dofs from the stream \p io and assigns the values
   * to the appropriate entries of \p vec.
   *
   * \returns The number of dofs read.
   *
   * Reads data and discards it if \p vec is a null pointer.
   */
  unsigned int read_SCALAR_dofs (const unsigned int var,
                                 Xdr & io,
                                 NumericVector<Number> * vec) const;

  /**
   * Reads a vector for this System.
   * This method may safely be called on a distributed-memory mesh.
   *
   * \returns The length of the vector read.
   *
   * Reads data and discards it if \p vec is a null pointer.
   */
  template <typename InValType>
  numeric_index_type read_serialized_vector (Xdr & io,
                                             NumericVector<Number> * vec);

  /**
   * Non-templated version for backward compatibility.
   *
   * Reads a vector for this System.
   * This method may safely be called on a distributed-memory mesh.
   *
   * \returns The length of the vector read.
   */
  numeric_index_type read_serialized_vector (Xdr & io,
                                             NumericVector<Number> & vec)
  { return read_serialized_vector<Number>(io, &vec); }

  /**
   * Writes an output vector to the stream \p io for a set of \p DofObjects.
   * This method uses blocked output and is safe to call on a distributed memory-mesh.
   *
   * \returns The number of values written
   */
  template <typename iterator_type>
  std::size_t write_serialized_blocked_dof_objects (const std::vector<const NumericVector<Number> *> & vecs,
                                                    const dof_id_type n_objects,
                                                    const iterator_type begin,
                                                    const iterator_type end,
                                                    Xdr & io,
                                                    const unsigned int var_to_write=libMesh::invalid_uint) const;

  /**
   * Writes the SCALAR dofs associated with var to the stream \p io.
   *
   * \returns The number of values written.
   */
  unsigned int write_SCALAR_dofs (const NumericVector<Number> & vec,
                                  const unsigned int var,
                                  Xdr & io) const;

  /**
   * Writes a vector for this System.
   * This method may safely be called on a distributed-memory mesh.
   *
   * \returns The number of values written.
   */
  dof_id_type write_serialized_vector (Xdr & io,
                                       const NumericVector<Number> & vec) const;

  /**
   * Function that initializes the system.
   */
  void (* _init_system_function) (EquationSystems & es,
                                  const std::string & name);

  /**
   * Object that initializes the system.
   */
  Initialization * _init_system_object;

  /**
   * Function that assembles the system.
   */
  void (* _assemble_system_function) (EquationSystems & es,
                                      const std::string & name);

  /**
   * Object that assembles the system.
   */
  Assembly * _assemble_system_object;

  /**
   * Function to impose constraints.
   */
  void (* _constrain_system_function) (EquationSystems & es,
                                       const std::string & name);

  /**
   * Object that constrains the system.
   */
  Constraint * _constrain_system_object;

  /**
   * Function to evaluate quantity of interest
   */
  void (* _qoi_evaluate_function) (EquationSystems & es,
                                   const std::string & name,
                                   const QoISet & qoi_indices);

  /**
   * Object to compute quantities of interest.
   */
  QOI * _qoi_evaluate_object;

  /**
   * Function to evaluate quantity of interest derivative
   */
  void (* _qoi_evaluate_derivative_function) (EquationSystems & es,
                                              const std::string & name,
                                              const QoISet & qoi_indices,
                                              bool include_liftfunc,
                                              bool apply_constraints);

  /**
   * Object to compute derivatives of quantities of interest.
   */
  QOIDerivative * _qoi_evaluate_derivative_object;

  /**
   * Data structure describing the relationship between
   * nodes, variables, etc... and degrees of freedom.
   */
  std::unique_ptr<DofMap> _dof_map;

  /**
   * Constant reference to the \p EquationSystems object
   * used for the simulation.
   */
  EquationSystems & _equation_systems;

  /**
   * Constant reference to the \p mesh data structure used
   * for the simulation.
   */
  MeshBase & _mesh;

  /**
   * A name associated with this system.
   */
  const std::string _sys_name;

  /**
   * The number associated with this system
   */
  const unsigned int _sys_number;

  /**
   * The \p Variable in this \p System.
   */
  std::vector<Variable> _variables;

  /**
   * The \p VariableGroup in this \p System.
   */
  std::vector<VariableGroup> _variable_groups;

  /**
   * The variable numbers corresponding to user-specified
   * names, useful for name-based lookups.
   */
  std::map<std::string, unsigned short int> _variable_numbers;

  /**
   * Flag stating if the system is active or not.
   */
  bool _active;

  /**
   * Some systems need an arbitrary number of vectors.
   * This map allows names to be associated with arbitrary
   * vectors.  All the vectors in this map will be distributed
   * in the same way as the solution vector.
   */
  std::map<std::string, NumericVector<Number> * > _vectors;

  /**
   * Holds true if a vector by that name should be projected
   * onto a changed grid, false if it should be zeroed.
   */
  std::map<std::string, bool> _vector_projections;

  /**
   * Holds non-negative if a vector by that name should be projected
   * using adjoint constraints/BCs, -1 if primal
   */
  std::map<std::string, int> _vector_is_adjoint;

  /**
   * Holds the type of a vector
   */
  std::map<std::string, ParallelType> _vector_types;

  /**
   * Holds true if the solution vector should be projected
   * onto a changed grid, false if it should be zeroed.
   * This is true by default.
   */
  bool _solution_projection;

  /**
   * Holds true if the components of more advanced system types (e.g.
   * system matrices) should not be initialized.
   */
  bool _basic_system_only;

  /**
   * \p true when additional vectors and variables do not require
   * immediate initialization, \p false otherwise.
   */
  bool _is_initialized;

  /**
   * \p true when \p VariableGroup structures should be automatically
   * identified, \p false otherwise.  Defaults to \p true.
   */
  bool _identify_variable_groups;

  /**
   * This flag is used only when *reading* in a system from file.
   * Based on the system header, it keeps track of how many
   * additional vectors were actually written for this file.
   */
  unsigned int _additional_data_written;

  /**
   * This vector is used only when *reading* in a system from file.
   * Based on the system header, it keeps track of any index remapping
   * between variable names in the data file and variable names in the
   * already-constructed system.  I.e. if we have a system with
   * variables "A1", "A2", "B1", and "B2", but we read in a data file with
   * only "A1" and "B1" defined, then we don't want to try and read in
   * A2 or B2, and we don't want to assign A1 and B1 values to
   * different dof indices.
   */
  std::vector<unsigned int> _written_var_indices;

  /**
   * Has the adjoint problem already been solved?  If the user sets
   * \p adjoint_already_solved to \p true, we won't waste time solving
   * it again.
   */
  bool adjoint_already_solved;

  /**
   * Are we allowed to write this system to file?  If \p _hide_output is
   * \p true, then \p EquationSystems::write will ignore this system.
   */
  bool _hide_output;
};



// ------------------------------------------------------------
// System inline methods
inline
const std::string & System::name() const
{
  return _sys_name;
}



inline
unsigned int System::number() const
{
  return _sys_number;
}



inline
const MeshBase & System::get_mesh() const
{
  return _mesh;
}



inline
MeshBase & System::get_mesh()
{
  return _mesh;
}



inline
const DofMap & System::get_dof_map() const
{
  return *_dof_map;
}



inline
DofMap & System::get_dof_map()
{
  return *_dof_map;
}



inline
bool System::active() const
{
  return _active;
}



inline
void System::activate ()
{
  _active = true;
}



inline
void System::deactivate ()
{
  _active = false;
}



inline
bool System::is_initialized ()
{
  return _is_initialized;
}



inline
void System::set_basic_system_only ()
{
  _basic_system_only = true;
}



inline
unsigned int System::n_vars() const
{
  return cast_int<unsigned int>(_variables.size());
}



inline
unsigned int System::n_variable_groups() const
{
  return cast_int<unsigned int>(_variable_groups.size());
}



inline
unsigned int System::n_components() const
{
  if (_variables.empty())
    return 0;

  const Variable & last = _variables.back();
  return last.first_scalar_number() + last.n_components();
}



inline
const Variable & System::variable (const unsigned int i) const
{
  libmesh_assert_less (i, _variables.size());

  return _variables[i];
}



inline
const VariableGroup & System::variable_group (const unsigned int vg) const
{
  libmesh_assert_less (vg, _variable_groups.size());

  return _variable_groups[vg];
}



inline
const std::string & System::variable_name (const unsigned int i) const
{
  libmesh_assert_less (i, _variables.size());

  return _variables[i].name();
}



inline
unsigned int
System::variable_scalar_number (const std::string & var,
                                unsigned int component) const
{
  return variable_scalar_number(this->variable_number(var), component);
}



inline
unsigned int
System::variable_scalar_number (unsigned int var_num,
                                unsigned int component) const
{
  return _variables[var_num].first_scalar_number() + component;
}



inline
const FEType & System::variable_type (const unsigned int i) const
{
  libmesh_assert_less (i, _variables.size());

  return _variables[i].type();
}



inline
const FEType & System::variable_type (const std::string & var) const
{
  return _variables[this->variable_number(var)].type();
}



inline
bool System::identify_variable_groups () const
{
  return _identify_variable_groups;
}



inline
void System::identify_variable_groups (const bool ivg)
{
  _identify_variable_groups = ivg;
}



inline
dof_id_type System::n_active_dofs() const
{
  return this->n_dofs() - this->n_constrained_dofs();
}



inline
bool System::have_vector (const std::string & vec_name) const
{
  return (_vectors.count(vec_name));
}



inline
unsigned int System::n_vectors () const
{
  return cast_int<unsigned int>(_vectors.size());
}

inline
unsigned int System::n_matrices () const
{
  return 0;
}

inline
System::vectors_iterator System::vectors_begin ()
{
  return _vectors.begin();
}

inline
System::const_vectors_iterator System::vectors_begin () const
{
  return _vectors.begin();
}

inline
System::vectors_iterator System::vectors_end ()
{
  return _vectors.end();
}

inline
System::const_vectors_iterator System::vectors_end () const
{
  return _vectors.end();
}

inline
void System::assemble_residual_derivatives (const ParameterVector &)
{
  libmesh_not_implemented();
}

inline
void System::disable_cache () { assemble_before_solve = true; }

inline
unsigned int System::n_qois() const
{
  return cast_int<unsigned int>(this->qoi.size());
}

inline
std::pair<unsigned int, Real>
System::sensitivity_solve (const ParameterVector &)
{
  libmesh_not_implemented();
}

inline
std::pair<unsigned int, Real>
System::weighted_sensitivity_solve (const ParameterVector &,
                                    const ParameterVector &)
{
  libmesh_not_implemented();
}

inline
std::pair<unsigned int, Real>
System::adjoint_solve (const QoISet &)
{
  libmesh_not_implemented();
}

inline
std::pair<unsigned int, Real>
System::weighted_sensitivity_adjoint_solve (const ParameterVector &,
                                            const ParameterVector &,
                                            const QoISet &)
{
  libmesh_not_implemented();
}

inline
void
System::adjoint_qoi_parameter_sensitivity (const QoISet &,
                                           const ParameterVector &,
                                           SensitivityData &)
{
  libmesh_not_implemented();
}

inline
void
System::forward_qoi_parameter_sensitivity (const QoISet &,
                                           const ParameterVector &,
                                           SensitivityData &)
{
  libmesh_not_implemented();
}

inline
void
System::qoi_parameter_hessian(const QoISet &,
                              const ParameterVector &,
                              SensitivityData &)
{
  libmesh_not_implemented();
}

inline
void
System::qoi_parameter_hessian_vector_product(const QoISet &,
                                             const ParameterVector &,
                                             const ParameterVector &,
                                             SensitivityData &)
{
  libmesh_not_implemented();
}


} // namespace libMesh

#endif // LIBMESH_SYSTEM_H
