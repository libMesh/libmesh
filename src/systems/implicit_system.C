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



// C++ includes

// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/implicit_system.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/linear_solver.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameters.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/qoi_set.h"
#include "libmesh/sensitivity_data.h"
#include "libmesh/sparse_matrix.h"

namespace libMesh
{

// ------------------------------------------------------------
// ImplicitSystem implementation
ImplicitSystem::ImplicitSystem (EquationSystems & es,
                                const std::string & name_in,
                                const unsigned int number_in) :

  Parent            (es, name_in, number_in),
  matrix            (libmesh_nullptr),
  zero_out_matrix_and_rhs(true),
  _can_add_matrices (true)
{
}



ImplicitSystem::~ImplicitSystem ()
{
  // Clear data
  this->clear();
}



void ImplicitSystem::clear ()
{
  // clear the parent data
  Parent::clear();

  // clear any user-added matrices
  {
    for (matrices_iterator pos = _matrices.begin();
         pos != _matrices.end(); ++pos)
      {
        pos->second->clear ();
        delete pos->second;
        pos->second = libmesh_nullptr;
      }

    _matrices.clear();
    _can_add_matrices = true;
  }

  // NULL the matrix.
  matrix = libmesh_nullptr;
}



void ImplicitSystem::init_data ()
{
  // initialize parent data
  Parent::init_data();

  // Clear any existing matrices
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    pos->second->clear();

  // Add the system matrix.
  this->add_system_matrix ();

  // Initialize the matrices for the system
  this->init_matrices ();
}



void ImplicitSystem::init_matrices ()
{
  libmesh_assert(matrix);

  // Check for quick return in case the system matrix
  // (and by extension all the matrices) has already
  // been initialized
  if (matrix->initialized())
    return;

  // Get a reference to the DofMap
  DofMap & dof_map = this->get_dof_map();

  // no chance to add other matrices
  _can_add_matrices = false;

  // Tell the matrices about the dof map, and vice versa
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    {
      SparseMatrix<Number> & m = *(pos->second);
      libmesh_assert (!m.initialized());

      // We want to allow repeated init() on systems, but we don't
      // want to attach the same matrix to the DofMap twice
      if (!dof_map.is_attached(m))
        dof_map.attach_matrix (m);
    }

  // Compute the sparsity pattern for the current
  // mesh and DOF distribution.  This also updates
  // additional matrices, \p DofMap now knows them
  dof_map.compute_sparsity (this->get_mesh());

  // Initialize matrices
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    pos->second->init ();

  // Set the additional matrices to 0.
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    pos->second->zero ();
}



void ImplicitSystem::reinit ()
{
  // initialize parent data
  Parent::reinit();

  // Get a reference to the DofMap
  DofMap & dof_map = this->get_dof_map();

  // Clear the matrices
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    {
      pos->second->clear();
      pos->second->attach_dof_map (dof_map);
    }

  // Clear the sparsity pattern
  this->get_dof_map().clear_sparsity();

  // Compute the sparsity pattern for the current
  // mesh and DOF distribution.  This also updates
  // additional matrices, \p DofMap now knows them
  dof_map.compute_sparsity (this->get_mesh());

  // Initialize matrices
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    pos->second->init ();

  // Set the additional matrices to 0.
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    pos->second->zero ();
}



void ImplicitSystem::assemble ()
{
  libmesh_assert(matrix);
  libmesh_assert (matrix->initialized());
  libmesh_assert(rhs);
  libmesh_assert (rhs->initialized());

  if (zero_out_matrix_and_rhs)
    {
      matrix->zero ();
      rhs->zero ();
    }

  // Call the base class assemble function
  Parent::assemble ();
}



SparseMatrix<Number> & ImplicitSystem::add_matrix (const std::string & mat_name)
{
  // only add matrices before initializing...
  if (!_can_add_matrices)
    libmesh_error_msg("ERROR: Too late.  Cannot add matrices to the system after initialization"
                      << "\n any more.  You should have done this earlier.");

  // Return the matrix if it is already there.
  if (this->have_matrix(mat_name))
    return *(_matrices[mat_name]);

  // Otherwise build the matrix and return it.
  SparseMatrix<Number> * buf = SparseMatrix<Number>::build(this->comm()).release();
  _matrices.insert (std::make_pair (mat_name, buf));

  return *buf;
}



const SparseMatrix<Number> * ImplicitSystem::request_matrix (const std::string & mat_name) const
{
  // Make sure the matrix exists
  const_matrices_iterator pos = _matrices.find (mat_name);

  if (pos == _matrices.end())
    return libmesh_nullptr;

  return pos->second;
}



SparseMatrix<Number> * ImplicitSystem::request_matrix (const std::string & mat_name)
{
  // Make sure the matrix exists
  matrices_iterator pos = _matrices.find (mat_name);

  if (pos == _matrices.end())
    return libmesh_nullptr;

  return pos->second;
}



const SparseMatrix<Number> & ImplicitSystem::get_matrix (const std::string & mat_name) const
{
  // Make sure the matrix exists
  const_matrices_iterator pos = _matrices.find (mat_name);

  if (pos == _matrices.end())
    libmesh_error_msg("ERROR: matrix " << mat_name << " does not exist in this system!");

  return *(pos->second);
}



SparseMatrix<Number> & ImplicitSystem::get_matrix (const std::string & mat_name)
{
  // Make sure the matrix exists
  matrices_iterator pos = _matrices.find (mat_name);

  if (pos == _matrices.end())
    libmesh_error_msg("ERROR: matrix " << mat_name << " does not exist in this system!");

  return *(pos->second);
}



void ImplicitSystem::add_system_matrix ()
{
  // Possible that we cleared the _matrices but
  // forgot to NULL-out the matrix?
  if (_matrices.empty()) matrix = libmesh_nullptr;


  // Only need to add the matrix if it isn't there
  // already!
  if (matrix == libmesh_nullptr)
    matrix = &(this->add_matrix ("System Matrix"));

  libmesh_assert(matrix);
}



void ImplicitSystem::disable_cache () {
  this->assemble_before_solve = true;
  this->get_linear_solver()->reuse_preconditioner(false);
}



std::pair<unsigned int, Real>
ImplicitSystem::sensitivity_solve (const ParameterVector & parameters)
{
  // Log how long the linear solve takes.
  LOG_SCOPE("sensitivity_solve()", "ImplicitSystem");

  // The forward system should now already be solved.
  // Now assemble the corresponding sensitivity system.

  if (this->assemble_before_solve)
    {
      // Build the Jacobian
      this->assembly(false, true);
      this->matrix->close();

      // Reset and build the RHS from the residual derivatives
      this->assemble_residual_derivatives(parameters);
    }

  // The sensitivity problem is linear
  LinearSolver<Number> * linear_solver = this->get_linear_solver();

  // Our iteration counts and residuals will be sums of the individual
  // results
  std::pair<unsigned int, Real> solver_params =
    this->get_linear_solve_parameters();
  std::pair<unsigned int, Real> totalrval = std::make_pair(0,0.0);

  // Solve the linear system.
  SparseMatrix<Number> * pc = this->request_matrix("Preconditioner");
  for (std::size_t p=0; p != parameters.size(); ++p)
    {
      std::pair<unsigned int, Real> rval =
        linear_solver->solve (*matrix, pc,
                              this->add_sensitivity_solution(p),
                              this->get_sensitivity_rhs(p),
                              solver_params.second,
                              solver_params.first);

      totalrval.first  += rval.first;
      totalrval.second += rval.second;
    }

  // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  for (std::size_t p=0; p != parameters.size(); ++p)
    this->get_dof_map().enforce_constraints_exactly
      (*this, &this->get_sensitivity_solution(p),
       /* homogeneous = */ true);
#endif

  this->release_linear_solver(linear_solver);

  return totalrval;
}



std::pair<unsigned int, Real>
ImplicitSystem::adjoint_solve (const QoISet & qoi_indices)
{
  // Log how long the linear solve takes.
  LOG_SCOPE("adjoint_solve()", "ImplicitSystem");

  if (this->assemble_before_solve)
    // Assemble the linear system
    this->assembly (/* get_residual = */ false,
                    /* get_jacobian = */ true);

  // The adjoint problem is linear
  LinearSolver<Number> * linear_solver = this->get_linear_solver();

  // Reset and build the RHS from the QOI derivative
  this->assemble_qoi_derivative(qoi_indices,
                                /* include_liftfunc = */ false,
                                /* apply_constraints = */ true);

  // Our iteration counts and residuals will be sums of the individual
  // results
  std::pair<unsigned int, Real> solver_params =
    this->get_linear_solve_parameters();
  std::pair<unsigned int, Real> totalrval = std::make_pair(0,0.0);

  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      {
        const std::pair<unsigned int, Real> rval =
          linear_solver->adjoint_solve (*matrix, this->add_adjoint_solution(i),
                                        this->get_adjoint_rhs(i),
                                        solver_params.second,
                                        solver_params.first);

        totalrval.first  += rval.first;
        totalrval.second += rval.second;
      }

  this->release_linear_solver(linear_solver);

  // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      this->get_dof_map().enforce_adjoint_constraints_exactly
        (this->get_adjoint_solution(i), i);
#endif

  return totalrval;
}



std::pair<unsigned int, Real>
ImplicitSystem::weighted_sensitivity_adjoint_solve (const ParameterVector & parameters_in,
                                                    const ParameterVector & weights,
                                                    const QoISet & qoi_indices)
{
  // Log how long the linear solve takes.
  LOG_SCOPE("weighted_sensitivity_adjoint_solve()", "ImplicitSystem");

  // We currently get partial derivatives via central differencing
  const Real delta_p = TOLERANCE;

  ParameterVector & parameters =
    const_cast<ParameterVector &>(parameters_in);

  // The forward system should now already be solved.
  // The adjoint system should now already be solved.
  // Now we're assembling a weighted sum of adjoint-adjoint systems:
  //
  // dR/du (u, sum_l(w_l*z^l)) = sum_l(w_l*(Q''_ul - R''_ul (u, z)))

  // FIXME: The derivation here does not yet take adjoint boundary
  // conditions into account.
  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      libmesh_assert(!this->get_dof_map().has_adjoint_dirichlet_boundaries(i));

  // We'll assemble the rhs first, because the R'' term will require
  // perturbing the jacobian

  // We'll use temporary rhs vectors, because we haven't (yet) found
  // any good reasons why users might want to save these:

  std::vector<NumericVector<Number> *> temprhs(this->qoi.size());
  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      temprhs[i] = this->rhs->zero_clone().release();

  // We approximate the _l partial derivatives via a central
  // differencing perturbation in the w_l direction:
  //
  // sum_l(w_l*v_l) ~= (v(p + dp*w_l*e_l) - v(p - dp*w_l*e_l))/(2*dp)

  // PETSc doesn't implement SGEMX, so neither does NumericVector,
  // so we want to avoid calculating f -= R'*z.  We'll thus evaluate
  // the above equation by first adding -v(p+dp...), then multiplying
  // the intermediate result vectors by -1, then adding -v(p-dp...),
  // then finally dividing by 2*dp.

  ParameterVector oldparameters, parameterperturbation;
  parameters.deep_copy(oldparameters);
  weights.deep_copy(parameterperturbation);
  parameterperturbation *= delta_p;
  parameters += parameterperturbation;

  this->assembly(false, true);
  this->matrix->close();

  // Take the discrete adjoint, so that we can calculate R_u(u,z) with
  // a matrix-vector product of R_u and z.
  matrix->get_transpose(*matrix);

  this->assemble_qoi_derivative(qoi_indices,
                                /* include_liftfunc = */ false,
                                /* apply_constraints = */ true);
  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      {
        this->get_adjoint_rhs(i).close();
        *(temprhs[i]) -= this->get_adjoint_rhs(i);
        this->matrix->vector_mult_add(*(temprhs[i]), this->get_adjoint_solution(i));
        *(temprhs[i]) *= -1.0;
      }

  oldparameters.value_copy(parameters);
  parameterperturbation *= -1.0;
  parameters += parameterperturbation;

  this->assembly(false, true);
  this->matrix->close();
  matrix->get_transpose(*matrix);

  this->assemble_qoi_derivative(qoi_indices,
                                /* include_liftfunc = */ false,
                                /* apply_constraints = */ true);
  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      {
        this->get_adjoint_rhs(i).close();
        *(temprhs[i]) -= this->get_adjoint_rhs(i);
        this->matrix->vector_mult_add(*(temprhs[i]), this->get_adjoint_solution(i));
        *(temprhs[i]) /= (2.0*delta_p);
      }

  // Finally, assemble the jacobian at the non-perturbed parameter
  // values.  Ignore assemble_before_solve; if we had a good
  // non-perturbed matrix before we've already overwritten it.
  oldparameters.value_copy(parameters);

  // if (this->assemble_before_solve)
  {
    // Build the Jacobian
    this->assembly(false, true);
    this->matrix->close();

    // Take the discrete adjoint
    matrix->get_transpose(*matrix);
  }

  // The weighted adjoint-adjoint problem is linear
  LinearSolver<Number> * linear_solver = this->get_linear_solver();

  // Our iteration counts and residuals will be sums of the individual
  // results
  std::pair<unsigned int, Real> solver_params =
    this->get_linear_solve_parameters();
  std::pair<unsigned int, Real> totalrval = std::make_pair(0,0.0);

  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      {
        const std::pair<unsigned int, Real> rval =
          linear_solver->solve (*matrix, this->add_weighted_sensitivity_adjoint_solution(i),
                                *(temprhs[i]),
                                solver_params.second,
                                solver_params.first);

        totalrval.first  += rval.first;
        totalrval.second += rval.second;
      }

  this->release_linear_solver(linear_solver);

  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      delete temprhs[i];

  // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      this->get_dof_map().enforce_constraints_exactly
        (*this, &this->get_weighted_sensitivity_adjoint_solution(i),
         /* homogeneous = */ true);
#endif

  return totalrval;
}



std::pair<unsigned int, Real>
ImplicitSystem::weighted_sensitivity_solve (const ParameterVector & parameters_in,
                                            const ParameterVector & weights)
{
  // Log how long the linear solve takes.
  LOG_SCOPE("weighted_sensitivity_solve()", "ImplicitSystem");

  // We currently get partial derivatives via central differencing
  const Real delta_p = TOLERANCE;

  ParameterVector & parameters =
    const_cast<ParameterVector &>(parameters_in);

  // The forward system should now already be solved.

  // Now we're assembling a weighted sum of sensitivity systems:
  //
  // dR/du (u, v)(sum(w_l*u'_l)) = -sum_l(w_l*R'_l (u, v)) forall v

  // We'll assemble the rhs first, because the R' term will require
  // perturbing the system, and some applications may not be able to
  // assemble a perturbed residual without simultaneously constructing
  // a perturbed jacobian.

  // We approximate the _l partial derivatives via a central
  // differencing perturbation in the w_l direction:
  //
  // sum_l(w_l*v_l) ~= (v(p + dp*w_l*e_l) - v(p - dp*w_l*e_l))/(2*dp)

  ParameterVector oldparameters, parameterperturbation;
  parameters.deep_copy(oldparameters);
  weights.deep_copy(parameterperturbation);
  parameterperturbation *= delta_p;
  parameters += parameterperturbation;

  this->assembly(true, false, true);
  this->rhs->close();

  UniquePtr<NumericVector<Number> > temprhs = this->rhs->clone();

  oldparameters.value_copy(parameters);
  parameterperturbation *= -1.0;
  parameters += parameterperturbation;

  this->assembly(true, false, true);
  this->rhs->close();

  *temprhs -= *(this->rhs);
  *temprhs /= (2.0*delta_p);

  // Finally, assemble the jacobian at the non-perturbed parameter
  // values
  oldparameters.value_copy(parameters);

  // Build the Jacobian
  this->assembly(false, true);
  this->matrix->close();

  // The weighted sensitivity problem is linear
  LinearSolver<Number> * linear_solver = this->get_linear_solver();

  std::pair<unsigned int, Real> solver_params =
    this->get_linear_solve_parameters();

  const std::pair<unsigned int, Real> rval =
    linear_solver->solve (*matrix, this->add_weighted_sensitivity_solution(),
                          *temprhs,
                          solver_params.second,
                          solver_params.first);

  this->release_linear_solver(linear_solver);

  // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  this->get_dof_map().enforce_constraints_exactly
    (*this, &this->get_weighted_sensitivity_solution(),
     /* homogeneous = */ true);
#endif

  return rval;
}



void ImplicitSystem::assemble_residual_derivatives(const ParameterVector & parameters_in)
{
  ParameterVector & parameters =
    const_cast<ParameterVector &>(parameters_in);

  const unsigned int Np = cast_int<unsigned int>
    (parameters.size());

  for (unsigned int p=0; p != Np; ++p)
    {
      NumericVector<Number> & sensitivity_rhs = this->add_sensitivity_rhs(p);

      // Approximate -(partial R / partial p) by
      // (R(p-dp) - R(p+dp)) / (2*dp)

      Number old_parameter = *parameters[p];

      const Real delta_p =
        TOLERANCE * std::max(std::abs(old_parameter), 1e-3);

      *parameters[p] -= delta_p;

      //      this->assembly(true, false, true);
      this->assembly(true, false, false);
      this->rhs->close();
      sensitivity_rhs = *this->rhs;

      *parameters[p] = old_parameter + delta_p;

      //      this->assembly(true, false, true);
      this->assembly(true, false, false);
      this->rhs->close();

      sensitivity_rhs -= *this->rhs;
      sensitivity_rhs /= (2*delta_p);
      sensitivity_rhs.close();

      *parameters[p] = old_parameter;
    }
}



void ImplicitSystem::adjoint_qoi_parameter_sensitivity (const QoISet & qoi_indices,
                                                        const ParameterVector & parameters_in,
                                                        SensitivityData & sensitivities)
{
  ParameterVector & parameters =
    const_cast<ParameterVector &>(parameters_in);

  const unsigned int Np = cast_int<unsigned int>
    (parameters.size());
  const unsigned int Nq = cast_int<unsigned int>
    (qoi.size());

  // An introduction to the problem:
  //
  // Residual R(u(p),p) = 0
  // partial R / partial u = J = system matrix
  //
  // This implies that:
  // d/dp(R) = 0
  // (partial R / partial p) +
  // (partial R / partial u) * (partial u / partial p) = 0

  // We first do an adjoint solve:
  // J^T * z = (partial q / partial u)
  // if we havent already or dont have an initial condition for the adjoint
  if (!this->is_adjoint_already_solved())
    {
      this->adjoint_solve(qoi_indices);
    }

  this->assemble_residual_derivatives(parameters_in);

  // Get ready to fill in senstivities:
  sensitivities.allocate_data(qoi_indices, *this, parameters);

  // We use the identities:
  // dq/dp = (partial q / partial p) + (partial q / partial u) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + (J^T * z) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + z * J *
  //         (partial u / partial p)

  // Leading to our final formula:
  // dq/dp = (partial q / partial p) - z * (partial R / partial p)

  // In the case of adjoints with heterogenous Dirichlet boundary
  // function phi, where
  // q :=  S(u) - R(u,phi)
  // the final formula works out to:
  // dq/dp = (partial S / partial p) - z * (partial R / partial p)
  // Because we currently have no direct access to
  // (partial S / partial p), we use the identity
  // (partial S / partial p) = (partial q / partial p) +
  //                           phi * (partial R / partial p)
  // to derive an equivalent equation:
  // dq/dp = (partial q / partial p) - (z-phi) * (partial R / partial p)

  // Since z-phi degrees of freedom are zero for constrained indices,
  // we can use the same constrained -(partial R / partial p) that we
  // use for forward sensitivity solves, taking into account the
  // differing sign convention.
  //
  // Since that vector is constrained, its constrained indices are
  // zero, so its product with phi is zero, so we can neglect the
  // evaluation of phi terms.

  for (unsigned int j=0; j != Np; ++j)
    {
      // We currently get partial derivatives via central differencing

      // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)
      // (partial R / partial p) ~= (rhs(p+dp) - rhs(p-dp))/(2*dp)

      Number old_parameter = *parameters[j];

      const Real delta_p =
        TOLERANCE * std::max(std::abs(old_parameter), 1e-3);

      *parameters[j] = old_parameter - delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number> qoi_minus = this->qoi;

      NumericVector<Number> & neg_partialR_partialp = this->get_sensitivity_rhs(j);

      *parameters[j] = old_parameter + delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number> & qoi_plus = this->qoi;

      std::vector<Number> partialq_partialp(Nq, 0);
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          partialq_partialp[i] = (qoi_plus[i] - qoi_minus[i]) / (2.*delta_p);

      // Don't leave the parameter changed
      *parameters[j] = old_parameter;

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          sensitivities[i][j] = partialq_partialp[i] +
            neg_partialR_partialp.dot(this->get_adjoint_solution(i));
    }

  // All parameters have been reset.
  // Reset the original qoi.

  this->assemble_qoi(qoi_indices);
}



void ImplicitSystem::forward_qoi_parameter_sensitivity (const QoISet & qoi_indices,
                                                        const ParameterVector & parameters_in,
                                                        SensitivityData & sensitivities)
{
  ParameterVector & parameters =
    const_cast<ParameterVector &>(parameters_in);

  const unsigned int Np = cast_int<unsigned int>
    (parameters.size());
  const unsigned int Nq = cast_int<unsigned int>
    (qoi.size());

  // An introduction to the problem:
  //
  // Residual R(u(p),p) = 0
  // partial R / partial u = J = system matrix
  //
  // This implies that:
  // d/dp(R) = 0
  // (partial R / partial p) +
  // (partial R / partial u) * (partial u / partial p) = 0

  // We first solve for (partial u / partial p) for each parameter:
  // J * (partial u / partial p) = - (partial R / partial p)

  this->sensitivity_solve(parameters);

  // Get ready to fill in senstivities:
  sensitivities.allocate_data(qoi_indices, *this, parameters);

  // We use the identity:
  // dq/dp = (partial q / partial p) + (partial q / partial u) *
  //         (partial u / partial p)

  // We get (partial q / partial u) from the user
  this->assemble_qoi_derivative(qoi_indices,
                                /* include_liftfunc = */ true,
                                /* apply_constraints = */ false);

  // We don't need these to be closed() in this function, but libMesh
  // standard practice is to have them closed() by the time the
  // function exits
  for (std::size_t i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      this->get_adjoint_rhs(i).close();

  for (unsigned int j=0; j != Np; ++j)
    {
      // We currently get partial derivatives via central differencing

      // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)

      Number old_parameter = *parameters[j];

      const Real delta_p =
        TOLERANCE * std::max(std::abs(old_parameter), 1e-3);

      *parameters[j] = old_parameter - delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number> qoi_minus = this->qoi;

      *parameters[j] = old_parameter + delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number> & qoi_plus = this->qoi;

      std::vector<Number> partialq_partialp(Nq, 0);
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          partialq_partialp[i] = (qoi_plus[i] - qoi_minus[i]) / (2.*delta_p);

      // Don't leave the parameter changed
      *parameters[j] = old_parameter;

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          sensitivities[i][j] = partialq_partialp[i] +
            this->get_adjoint_rhs(i).dot(this->get_sensitivity_solution(j));
    }

  // All parameters have been reset.
  // We didn't cache the original rhs or matrix for memory reasons,
  // but we can restore them to a state consistent solution -
  // principle of least surprise.
  this->assembly(true, true);
  this->rhs->close();
  this->matrix->close();
  this->assemble_qoi(qoi_indices);
}



void ImplicitSystem::qoi_parameter_hessian_vector_product (const QoISet & qoi_indices,
                                                           const ParameterVector & parameters_in,
                                                           const ParameterVector & vector,
                                                           SensitivityData & sensitivities)
{
  // We currently get partial derivatives via finite differencing
  const Real delta_p = TOLERANCE;

  ParameterVector & parameters =
    const_cast<ParameterVector &>(parameters_in);

  // We'll use a single temporary vector for matrix-vector-vector products
  UniquePtr<NumericVector<Number> > tempvec = this->solution->zero_clone();

  const unsigned int Np = cast_int<unsigned int>
    (parameters.size());
  const unsigned int Nq = cast_int<unsigned int>
    (qoi.size());

  // For each quantity of interest q, the parameter sensitivity
  // Hessian is defined as q''_{kl} = {d^2 q}/{d p_k d p_l}.
  // Given a vector of parameter perturbation weights w_l, this
  // function evaluates the hessian-vector product sum_l(q''_{kl}*w_l)
  //
  // We calculate it from values and partial derivatives of the
  // quantity of interest function Q, solution u, adjoint solution z,
  // parameter sensitivity adjoint solutions z^l, and residual R, as:
  //
  // sum_l(q''_{kl}*w_l) =
  // sum_l(w_l * Q''_{kl}) + Q''_{uk}(u)*(sum_l(w_l u'_l)) -
  // R'_k(u, sum_l(w_l*z^l)) - R'_{uk}(u,z)*(sum_l(w_l u'_l) -
  // sum_l(w_l*R''_{kl}(u,z))
  //
  // See the adjoints model document for more details.

  // We first do an adjoint solve to get z for each quantity of
  // interest
  // if we havent already or dont have an initial condition for the adjoint
  if (!this->is_adjoint_already_solved())
    {
      this->adjoint_solve(qoi_indices);
    }

  // Get ready to fill in senstivities:
  sensitivities.allocate_data(qoi_indices, *this, parameters);

  // We can't solve for all the solution sensitivities u'_l or for all
  // of the parameter sensitivity adjoint solutions z^l without
  // requiring O(Nq*Np) linear solves.  So we'll solve directly for their
  // weighted sum - this is just O(Nq) solves.

  // First solve for sum_l(w_l u'_l).
  this->weighted_sensitivity_solve(parameters, vector);

  // Then solve for sum_l(w_l z^l).
  this->weighted_sensitivity_adjoint_solve(parameters, vector, qoi_indices);

  for (unsigned int k=0; k != Np; ++k)
    {
      // We approximate sum_l(w_l * Q''_{kl}) with a central
      // differencing pertubation:
      // sum_l(w_l * Q''_{kl}) ~=
      // (Q(p + dp*w_l*e_l + dp*e_k) - Q(p - dp*w_l*e_l + dp*e_k) -
      // Q(p + dp*w_l*e_l - dp*e_k) + Q(p - dp*w_l*e_l - dp*e_k))/(4*dp^2)

      // The sum(w_l*R''_kl) term requires the same sort of pertubation,
      // and so we subtract it in at the same time:
      // sum_l(w_l * R''_{kl}) ~=
      // (R(p + dp*w_l*e_l + dp*e_k) - R(p - dp*w_l*e_l + dp*e_k) -
      // R(p + dp*w_l*e_l - dp*e_k) + R(p - dp*w_l*e_l - dp*e_k))/(4*dp^2)

      ParameterVector oldparameters, parameterperturbation;
      parameters.deep_copy(oldparameters);
      vector.deep_copy(parameterperturbation);
      parameterperturbation *= delta_p;
      parameters += parameterperturbation;

      Number old_parameter = *parameters[k];

      *parameters[k] = old_parameter + delta_p;
      this->assemble_qoi(qoi_indices);
      this->assembly(true, false, true);
      this->rhs->close();
      std::vector<Number> partial2q_term = this->qoi;
      std::vector<Number> partial2R_term(this->qoi.size());
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          partial2R_term[i] = this->rhs->dot(this->get_adjoint_solution(i));

      *parameters[k] = old_parameter - delta_p;
      this->assemble_qoi(qoi_indices);
      this->assembly(true, false, true);
      this->rhs->close();
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          {
            partial2q_term[i] -= this->qoi[i];
            partial2R_term[i] -= this->rhs->dot(this->get_adjoint_solution(i));
          }

      oldparameters.value_copy(parameters);
      parameterperturbation *= -1.0;
      parameters += parameterperturbation;

      // Re-center old_parameter, which may be affected by vector
      old_parameter = *parameters[k];

      *parameters[k] = old_parameter + delta_p;
      this->assemble_qoi(qoi_indices);
      this->assembly(true, false, true);
      this->rhs->close();
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          {
            partial2q_term[i] -= this->qoi[i];
            partial2R_term[i] -= this->rhs->dot(this->get_adjoint_solution(i));
          }

      *parameters[k] = old_parameter - delta_p;
      this->assemble_qoi(qoi_indices);
      this->assembly(true, false, true);
      this->rhs->close();
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          {
            partial2q_term[i] += this->qoi[i];
            partial2R_term[i] += this->rhs->dot(this->get_adjoint_solution(i));
          }

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          {
            partial2q_term[i] /= (4. * delta_p * delta_p);
            partial2R_term[i] /= (4. * delta_p * delta_p);
          }

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          sensitivities[i][k] = partial2q_term[i] - partial2R_term[i];

      // We get (partial q / partial u), R, and
      // (partial R / partial u) from the user, but centrally
      // difference to get q_uk, R_k, and R_uk terms:
      // (partial R / partial k)
      // R_k*sum(w_l*z^l) = (R(p+dp*e_k)*sum(w_l*z^l) - R(p-dp*e_k)*sum(w_l*z^l))/(2*dp)
      // (partial^2 q / partial u partial k)
      // q_uk = (q_u(p+dp*e_k) - q_u(p-dp*e_k))/(2*dp)
      // (partial^2 R / partial u partial k)
      // R_uk*z*sum(w_l*u'_l) = (R_u(p+dp*e_k)*z*sum(w_l*u'_l) - R_u(p-dp*e_k)*z*sum(w_l*u'_l))/(2*dp)

      // To avoid creating Nq temporary vectors for q_uk or R_uk, we add
      // subterms to the sensitivities output one by one.
      //
      // FIXME: this is probably a bad order of operations for
      // controlling floating point error.

      *parameters[k] = old_parameter + delta_p;
      this->assembly(true, true);
      this->rhs->close();
      this->matrix->close();
      this->assemble_qoi_derivative(qoi_indices,
                                    /* include_liftfunc = */ true,
                                    /* apply_constraints = */ false);

      this->matrix->vector_mult(*tempvec, this->get_weighted_sensitivity_solution());

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          {
            this->get_adjoint_rhs(i).close();
            sensitivities[i][k] += (this->get_adjoint_rhs(i).dot(this->get_weighted_sensitivity_solution()) -
                                    this->rhs->dot(this->get_weighted_sensitivity_adjoint_solution(i)) -
                                    this->get_adjoint_solution(i).dot(*tempvec)) / (2.*delta_p);
          }

      *parameters[k] = old_parameter - delta_p;
      this->assembly(true, true);
      this->rhs->close();
      this->matrix->close();
      this->assemble_qoi_derivative(qoi_indices,
                                    /* include_liftfunc = */ true,
                                    /* apply_constraints = */ false);

      this->matrix->vector_mult(*tempvec, this->get_weighted_sensitivity_solution());

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          {
            this->get_adjoint_rhs(i).close();
            sensitivities[i][k] += (-this->get_adjoint_rhs(i).dot(this->get_weighted_sensitivity_solution()) +
                                    this->rhs->dot(this->get_weighted_sensitivity_adjoint_solution(i)) +
                                    this->get_adjoint_solution(i).dot(*tempvec)) / (2.*delta_p);
          }
    }

  // All parameters have been reset.
  // Don't leave the qoi or system changed - principle of least
  // surprise.
  this->assembly(true, true);
  this->rhs->close();
  this->matrix->close();
  this->assemble_qoi(qoi_indices);
}



void ImplicitSystem::qoi_parameter_hessian (const QoISet & qoi_indices,
                                            const ParameterVector & parameters_in,
                                            SensitivityData & sensitivities)
{
  // We currently get partial derivatives via finite differencing
  const Real delta_p = TOLERANCE;

  ParameterVector & parameters =
    const_cast<ParameterVector &>(parameters_in);

  // We'll use one temporary vector for matrix-vector-vector products
  UniquePtr<NumericVector<Number> > tempvec = this->solution->zero_clone();

  // And another temporary vector to hold a copy of the true solution
  // so we can safely perturb this->solution.
  UniquePtr<NumericVector<Number> > oldsolution = this->solution->clone();

  const unsigned int Np = cast_int<unsigned int>
    (parameters.size());
  const unsigned int Nq = cast_int<unsigned int>
    (qoi.size());

  // For each quantity of interest q, the parameter sensitivity
  // Hessian is defined as q''_{kl} = {d^2 q}/{d p_k d p_l}.
  //
  // We calculate it from values and partial derivatives of the
  // quantity of interest function Q, solution u, adjoint solution z,
  // and residual R, as:
  //
  // q''_{kl} =
  // Q''_{kl} + Q''_{uk}(u)*u'_l + Q''_{ul}(u) * u'_k +
  // Q''_{uu}(u)*u'_k*u'_l -
  // R''_{kl}(u,z) -
  // R''_{uk}(u,z)*u'_l - R''_{ul}(u,z)*u'_k -
  // R''_{uu}(u,z)*u'_k*u'_l
  //
  // See the adjoints model document for more details.

  // We first do an adjoint solve to get z for each quantity of
  // interest
  // if we havent already or dont have an initial condition for the adjoint
  if (!this->is_adjoint_already_solved())
    {
      this->adjoint_solve(qoi_indices);
    }

  // And a sensitivity solve to get u_k for each parameter
  this->sensitivity_solve(parameters);

  // Get ready to fill in second derivatives:
  sensitivities.allocate_hessian_data(qoi_indices, *this, parameters);

  for (unsigned int k=0; k != Np; ++k)
    {
      Number old_parameterk = *parameters[k];

      // The Hessian is symmetric, so we just calculate the lower
      // triangle and the diagonal, and we get the upper triangle from
      // the transpose of the lower

      for (unsigned int l=0; l != k+1; ++l)
        {
          // The second partial derivatives with respect to parameters
          // are all calculated via a central finite difference
          // stencil:
          // F''_{kl} ~= (F(p+dp*e_k+dp*e_l) - F(p+dp*e_k-dp*e_l) -
          //              F(p-dp*e_k+dp*e_l) + F(p-dp*e_k-dp*e_l))/(4*dp^2)
          // We will add Q''_{kl}(u) and subtract R''_{kl}(u,z) at the
          // same time.
          //
          // We have to be careful with the perturbations to handle
          // the k=l case

          Number old_parameterl = *parameters[l];

          *parameters[k] += delta_p;
          *parameters[l] += delta_p;
          this->assemble_qoi(qoi_indices);
          this->assembly(true, false, true);
          this->rhs->close();
          std::vector<Number> partial2q_term = this->qoi;
          std::vector<Number> partial2R_term(this->qoi.size());
          for (unsigned int i=0; i != Nq; ++i)
            if (qoi_indices.has_index(i))
              partial2R_term[i] = this->rhs->dot(this->get_adjoint_solution(i));

          *parameters[l] -= 2.*delta_p;
          this->assemble_qoi(qoi_indices);
          this->assembly(true, false, true);
          this->rhs->close();
          for (unsigned int i=0; i != Nq; ++i)
            if (qoi_indices.has_index(i))
              {
                partial2q_term[i] -= this->qoi[i];
                partial2R_term[i] -= this->rhs->dot(this->get_adjoint_solution(i));
              }

          *parameters[k] -= 2.*delta_p;
          this->assemble_qoi(qoi_indices);
          this->assembly(true, false, true);
          this->rhs->close();
          for (unsigned int i=0; i != Nq; ++i)
            if (qoi_indices.has_index(i))
              {
                partial2q_term[i] += this->qoi[i];
                partial2R_term[i] += this->rhs->dot(this->get_adjoint_solution(i));
              }

          *parameters[l] += 2.*delta_p;
          this->assemble_qoi(qoi_indices);
          this->assembly(true, false, true);
          this->rhs->close();
          for (unsigned int i=0; i != Nq; ++i)
            if (qoi_indices.has_index(i))
              {
                partial2q_term[i] -= this->qoi[i];
                partial2R_term[i] -= this->rhs->dot(this->get_adjoint_solution(i));
                partial2q_term[i] /= (4. * delta_p * delta_p);
                partial2R_term[i] /= (4. * delta_p * delta_p);
              }

          for (unsigned int i=0; i != Nq; ++i)
            if (qoi_indices.has_index(i))
              {
                Number current_terms = partial2q_term[i] - partial2R_term[i];
                sensitivities.second_derivative(i,k,l) += current_terms;
                if (k != l)
                  sensitivities.second_derivative(i,l,k) += current_terms;
              }

          // Don't leave the parameters perturbed
          *parameters[l] = old_parameterl;
          *parameters[k] = old_parameterk;
        }

      // We get (partial q / partial u) and
      // (partial R / partial u) from the user, but centrally
      // difference to get q_uk and R_uk terms:
      // (partial^2 q / partial u partial k)
      // q_uk*u'_l = (q_u(p+dp*e_k)*u'_l - q_u(p-dp*e_k)*u'_l)/(2*dp)
      // R_uk*z*u'_l = (R_u(p+dp*e_k)*z*u'_l - R_u(p-dp*e_k)*z*u'_l)/(2*dp)
      //
      // To avoid creating Nq temporary vectors, we add these
      // subterms to the sensitivities output one by one.
      //
      // FIXME: this is probably a bad order of operations for
      // controlling floating point error.

      *parameters[k] = old_parameterk + delta_p;
      this->assembly(false, true);
      this->matrix->close();
      this->assemble_qoi_derivative(qoi_indices,
                                    /* include_liftfunc = */ true,
                                    /* apply_constraints = */ false);

      for (unsigned int l=0; l != Np; ++l)
        {
          this->matrix->vector_mult(*tempvec, this->get_sensitivity_solution(l));
          for (unsigned int i=0; i != Nq; ++i)
            if (qoi_indices.has_index(i))
              {
                this->get_adjoint_rhs(i).close();
                Number current_terms =
                  (this->get_adjoint_rhs(i).dot(this->get_sensitivity_solution(l)) -
                   tempvec->dot(this->get_adjoint_solution(i))) / (2.*delta_p);
                sensitivities.second_derivative(i,k,l) += current_terms;

                // We use the _uk terms twice; symmetry lets us reuse
                // these calculations for the _ul terms.

                sensitivities.second_derivative(i,l,k) += current_terms;
              }
        }

      *parameters[k] = old_parameterk - delta_p;
      this->assembly(false, true);
      this->matrix->close();
      this->assemble_qoi_derivative(qoi_indices,
                                    /* include_liftfunc = */ true,
                                    /* apply_constraints = */ false);

      for (unsigned int l=0; l != Np; ++l)
        {
          this->matrix->vector_mult(*tempvec, this->get_sensitivity_solution(l));
          for (unsigned int i=0; i != Nq; ++i)
            if (qoi_indices.has_index(i))
              {
                this->get_adjoint_rhs(i).close();
                Number current_terms =
                  (-this->get_adjoint_rhs(i).dot(this->get_sensitivity_solution(l)) +
                   tempvec->dot(this->get_adjoint_solution(i))) / (2.*delta_p);
                sensitivities.second_derivative(i,k,l) += current_terms;

                // We use the _uk terms twice; symmetry lets us reuse
                // these calculations for the _ul terms.

                sensitivities.second_derivative(i,l,k) += current_terms;
              }
        }

      // Don't leave the parameter perturbed
      *parameters[k] = old_parameterk;

      // Our last remaining terms are -R_uu(u,z)*u_k*u_l and
      // Q_uu(u)*u_k*u_l
      //
      // We take directional central finite differences of R_u and Q_u
      // to approximate these terms, e.g.:
      //
      // Q_uu(u)*u_k ~= (Q_u(u+dp*u_k) - Q_u(u-dp*u_k))/(2*dp)

      *this->solution = this->get_sensitivity_solution(k);
      *this->solution *= delta_p;
      *this->solution += *oldsolution;
      this->assembly(false, true);
      this->matrix->close();
      this->assemble_qoi_derivative(qoi_indices,
                                    /* include_liftfunc = */ true,
                                    /* apply_constraints = */ false);

      // The Hessian is symmetric, so we just calculate the lower
      // triangle and the diagonal, and we get the upper triangle from
      // the transpose of the lower
      //
      // Note that, because we took the directional finite difference
      // with respect to k and not l, we've added an O(delta_p^2)
      // error to any permutational symmetry in the Hessian...
      for (unsigned int l=0; l != k+1; ++l)
        {
          this->matrix->vector_mult(*tempvec, this->get_sensitivity_solution(l));
          for (unsigned int i=0; i != Nq; ++i)
            if (qoi_indices.has_index(i))
              {
                this->get_adjoint_rhs(i).close();
                Number current_terms =
                  (this->get_adjoint_rhs(i).dot(this->get_sensitivity_solution(l)) -
                   tempvec->dot(this->get_adjoint_solution(i))) / (2.*delta_p);
                sensitivities.second_derivative(i,k,l) += current_terms;
                if (k != l)
                  sensitivities.second_derivative(i,l,k) += current_terms;
              }
        }

      *this->solution = this->get_sensitivity_solution(k);
      *this->solution *= -delta_p;
      *this->solution += *oldsolution;
      this->assembly(false, true);
      this->matrix->close();
      this->assemble_qoi_derivative(qoi_indices,
                                    /* include_liftfunc = */ true,
                                    /* apply_constraints = */ false);

      for (unsigned int l=0; l != k+1; ++l)
        {
          this->matrix->vector_mult(*tempvec, this->get_sensitivity_solution(l));
          for (unsigned int i=0; i != Nq; ++i)
            if (qoi_indices.has_index(i))
              {
                this->get_adjoint_rhs(i).close();
                Number current_terms =
                  (-this->get_adjoint_rhs(i).dot(this->get_sensitivity_solution(l)) +
                   tempvec->dot(this->get_adjoint_solution(i))) / (2.*delta_p);
                sensitivities.second_derivative(i,k,l) += current_terms;
                if (k != l)
                  sensitivities.second_derivative(i,l,k) += current_terms;
              }
        }

      // Don't leave the solution perturbed
      *this->solution = *oldsolution;
    }

  // All parameters have been reset.
  // Don't leave the qoi or system changed - principle of least
  // surprise.
  this->assembly(true, true);
  this->rhs->close();
  this->matrix->close();
  this->assemble_qoi(qoi_indices);
}



LinearSolver<Number> * ImplicitSystem::get_linear_solver() const
{
  // This function allocates memory and hands it back to the user as a
  // naked pointer.  This makes it too easy to leak memory, and
  // therefore this function is deprecated.  After a period of
  // deprecation, this function will eventually be marked with a
  // libmesh_error_msg().
  libmesh_deprecated();
  // libmesh_error_msg("This function should be overridden by derived classes. "
  //                   "It does not contain a valid LinearSolver to hand back to "
  //                   "the user, so it creates one, opening up the possibility "
  //                   "of a memory leak.");

  LinearSolver<Number> * new_solver =
    LinearSolver<Number>::build(this->comm()).release();

  if (libMesh::on_command_line("--solver_system_names"))
    new_solver->init((this->name()+"_").c_str());
  else
    new_solver->init();

  return new_solver;
}



std::pair<unsigned int, Real> ImplicitSystem::get_linear_solve_parameters() const
{
  return std::make_pair(this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations"),
                        this->get_equation_systems().parameters.get<Real>("linear solver tolerance"));
}



void ImplicitSystem::release_linear_solver(LinearSolver<Number> * s) const
{
  // This is the counterpart of the get_linear_solver() function, which is now deprecated.
  libmesh_deprecated();

  delete s;
}

} // namespace libMesh
