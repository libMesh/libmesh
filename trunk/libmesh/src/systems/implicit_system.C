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



// C++ includes

// Local includes
#include "dof_map.h"
#include "equation_systems.h"
#include "implicit_system.h"
#include "libmesh_logging.h"
#include "linear_solver.h"
#include "mesh.h"
#include "numeric_vector.h"
#include "parameters.h"
#include "parameter_vector.h"
#include "qoi_set.h"
#include "sensitivity_data.h"
#include "sparse_matrix.h"

// ------------------------------------------------------------
// ImplicitSystem implementation
ImplicitSystem::ImplicitSystem (EquationSystems& es,
				const std::string& name,
				const unsigned int number) :
  
  Parent            (es, name, number),
  matrix            (NULL),
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
	pos->second = NULL;
      }

    _matrices.clear();
    _can_add_matrices = true;
  }

  // NULL the matrix.
  matrix = NULL;
}



void ImplicitSystem::init_data ()
{
  // initialize parent data
  Parent::init_data();

  // Add the system matrix.
  this->add_system_matrix ();
  
  // Initialize the matrices for the system
  this->init_matrices ();
}



void ImplicitSystem::init_matrices ()
{
  libmesh_assert (matrix != NULL);

  // Check for quick return in case the system matrix
  // (and by extension all the matrices) has already
  // been initialized
  if (matrix->initialized())
    return;

  // Get a reference to the DofMap
  DofMap& dof_map = this->get_dof_map();
  
  // no chance to add other matrices
  _can_add_matrices = false;
  
  // Tell the matrices about the dof map, and vice versa
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    {
      libmesh_assert (!pos->second->initialized());
      dof_map.attach_matrix (*(pos->second));
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

  // Clear the matrices
  for (matrices_iterator pos = _matrices.begin();
       pos != _matrices.end(); ++pos)
    pos->second->clear();

  // Re-initialize the matrices
  this->init_matrices ();
}



void ImplicitSystem::assemble ()
{
  libmesh_assert (matrix != NULL);
  libmesh_assert (matrix->initialized());
  libmesh_assert (rhs    != NULL);
  libmesh_assert (rhs->initialized());

  // The user assembly gets to expect to accumulate on an initially
  // empty system
  matrix->zero ();
  rhs->zero ();

  // Call the base class assemble function
  Parent::assemble ();
}



SparseMatrix<Number> & ImplicitSystem::add_matrix (const std::string& mat_name)
{
  // only add matrices before initializing...
  if (!_can_add_matrices)
    {
      std::cerr << "ERROR: Too late.  Cannot add matrices to the system after initialization"
		<< std::endl
		<< " any more.  You should have done this earlier."
		<< std::endl;
      libmesh_error();
    }

  // Return the matrix if it is already there.
  if (this->have_matrix(mat_name))
    return *(_matrices[mat_name]);

  // Otherwise build the matrix and return it.
  SparseMatrix<Number>* buf = SparseMatrix<Number>::build().release();
  _matrices.insert (std::make_pair (mat_name, buf));

  return *buf;
}



const SparseMatrix<Number> * ImplicitSystem::request_matrix (const std::string& mat_name) const
{
  // Make sure the matrix exists
  const_matrices_iterator pos = _matrices.find (mat_name);
  
  if (pos == _matrices.end())
    return NULL;
  
  return pos->second;
}



SparseMatrix<Number> * ImplicitSystem::request_matrix (const std::string& mat_name)
{
  // Make sure the matrix exists
  matrices_iterator pos = _matrices.find (mat_name);
  
  if (pos == _matrices.end())
    return NULL;
  
  return pos->second;
}



const SparseMatrix<Number> & ImplicitSystem::get_matrix (const std::string& mat_name) const
{
  // Make sure the matrix exists
  const_matrices_iterator pos = _matrices.find (mat_name);
  
  if (pos == _matrices.end())
    {
      std::cerr << "ERROR: matrix "
		<< mat_name
		<< " does not exist in this system!"
		<< std::endl;      
      libmesh_error();
    }
  
  return *(pos->second);
}



SparseMatrix<Number> & ImplicitSystem::get_matrix (const std::string& mat_name)
{
  // Make sure the matrix exists
  matrices_iterator pos = _matrices.find (mat_name);
  
  if (pos == _matrices.end())
    {
      std::cerr << "ERROR: matrix "
		<< mat_name
		<< " does not exist in this system!"
		<< std::endl;      
      libmesh_error();
    }
  
  return *(pos->second);
}



void ImplicitSystem::add_system_matrix ()
{
  // Possible that we cleared the _matrices but
  // forgot to NULL-out the matrix?
  if (_matrices.empty()) matrix = NULL;


  // Only need to add the matrix if it isn't there
  // already!
  if (matrix == NULL)
    matrix = &(this->add_matrix ("System Matrix"));

  libmesh_assert (matrix != NULL);
}



std::pair<unsigned int, Real>
ImplicitSystem::sensitivity_solve (const ParameterVector& parameters)
{
  // Log how long the linear solve takes.
  START_LOG("sensitivity_solve()", "ImplicitSystem");

  // The forward system should now already be solved.
  // Now assemble the corresponding sensitivity system.

  if (this->assemble_before_solve)
    {
      // Build the Jacobian
      this->assembly(false, true);

      // Reset and build the RHS from the residual derivatives
      this->assemble_residual_derivatives(parameters);
    }

  // The sensitivity problem is linear
  LinearSolver<Number> *linear_solver = this->get_linear_solver();

  // Our iteration counts and residuals will be sums of the individual
  // results
  std::pair<unsigned int, Real> solver_params =
    this->get_linear_solve_parameters();
  std::pair<unsigned int, Real> totalrval = std::make_pair(0,0.0);

  // Solve the linear system.
  SparseMatrix<Number> *pc = this->request_matrix("Preconditioner");
  for (unsigned int p=0; p != parameters.size(); ++p)
    {
      std::pair<unsigned int, Real> rval =
        linear_solver->solve (*matrix, pc,
                              this->get_sensitivity_solution(p),
                              this->get_sensitivity_rhs(p),
                              solver_params.second,
                              solver_params.first);

      totalrval.first  += rval.first;
      totalrval.second += rval.second;
    }

  // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_AMR
  for (unsigned int p=0; p != parameters.size(); ++p)
    this->get_dof_map().enforce_constraints_exactly
      (*this, &this->get_sensitivity_solution(p));
#endif

  this->release_linear_solver(linear_solver);

  // Stop logging the nonlinear solve
  STOP_LOG("sensitivity_solve()", "ImplicitSystem");

  return totalrval;
}



std::pair<unsigned int, Real>
ImplicitSystem::adjoint_solve (const QoISet& qoi_indices)
{
  // Log how long the linear solve takes.
  START_LOG("adjoint_solve()", "ImplicitSystem");

  // The forward system should now already be solved.
  // Now assemble it's adjoint.

  if (this->assemble_before_solve)
    {
      // Build the Jacobian
      this->assembly(false, true);

      // Take the discrete adjoint
      matrix->get_transpose(*matrix);

      // Reset and build the RHS from the QOI derivative
      this->assemble_qoi_derivative(qoi_indices);
    }

  // The adjoint problem is linear
  LinearSolver<Number> *linear_solver = this->get_linear_solver();

  // Our iteration counts and residuals will be sums of the individual
  // results
  std::pair<unsigned int, Real> solver_params =
    this->get_linear_solve_parameters();
  std::pair<unsigned int, Real> totalrval = std::make_pair(0,0.0);

  for (unsigned int i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      {
        const std::pair<unsigned int, Real> rval =
          linear_solver->solve (*matrix, this->add_adjoint_solution(i),
                                 this->get_adjoint_rhs(i),
                                 solver_params.second,
                                 solver_params.first);

        totalrval.first  += rval.first;
        totalrval.second += rval.second;
      }

  this->release_linear_solver(linear_solver);

  // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_AMR
  for (unsigned int i=0; i != this->qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      this->get_dof_map().enforce_constraints_exactly
        (*this, &this->get_adjoint_solution(i));
#endif

  // Stop logging the nonlinear solve
  STOP_LOG("adjoint_solve()", "ImplicitSystem");

  return totalrval;
}



void ImplicitSystem::assemble_residual_derivatives(const ParameterVector& parameters)
{
  const unsigned int Np = parameters.size();
  Real deltap = TOLERANCE;

  for (unsigned int p=0; p != Np; ++p)
    {
      NumericVector<Number> &sensitivity_rhs = this->add_sensitivity_rhs(p);

      // Approximate -(partial R / partial p) by
      // (R(p-dp) - R(p+dp)) / (2*dp)

      Number old_parameter = *parameters[p];
      *parameters[p] -= deltap;

      this->assembly(true, false);
      sensitivity_rhs = *this->rhs;

      *parameters[p] = old_parameter + deltap;

      this->assembly(true, false);

      sensitivity_rhs -= *this->rhs;
      sensitivity_rhs /= (2*deltap);
      sensitivity_rhs.close();

      *parameters[p] = old_parameter;
    }
}



void ImplicitSystem::adjoint_qoi_parameter_sensitivity
  (const QoISet&          qoi_indices,
   const ParameterVector& parameters,
   SensitivityData&       sensitivities)
{
  const unsigned int Np = parameters.size();
  const unsigned int Nq = qoi.size();

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

  this->adjoint_solve(qoi_indices);

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

  for (unsigned int j=0; j != Np; ++j)
    {
      // We currently get partial derivatives via central differencing
      Number delta_p = 1e-6;

      // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)
      // (partial R / partial p) ~= (rhs(p+dp) - rhs(p-dp))/(2*dp)

      Number old_parameter = *parameters[j];
      // Number old_qoi = this->qoi;

      *parameters[j] = old_parameter - delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number> qoi_minus = this->qoi;

      this->assembly(true, false);
      AutoPtr<NumericVector<Number> > partialR_partialp = this->rhs->clone();
      *partialR_partialp *= -1;

      *parameters[j] = old_parameter + delta_p;
      this->assemble_qoi(qoi_indices);
      std::vector<Number>& qoi_plus = this->qoi;

      std::vector<Number> partialq_partialp(Nq, 0);
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          partialq_partialp[i] = (qoi_plus[i] - qoi_minus[i]) / (2.*delta_p);

      this->assembly(true, false);
      *partialR_partialp += *this->rhs;
      *partialR_partialp /= (2.*delta_p);

      // Don't leave the parameter changed
      *parameters[j] = old_parameter;

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          sensitivities[i][j] = partialq_partialp[i] -
            partialR_partialp->dot(this->get_adjoint_solution(i));
    }

  // All parameters have been reset.
  // We didn't cache the original rhs or matrix for memory reasons,
  // but we can restore them to a state consistent solution -
  // principle of least surprise.
  this->assembly(true, true);
  this->assemble_qoi(qoi_indices);
}



void ImplicitSystem::forward_qoi_parameter_sensitivity
  (const QoISet&          qoi_indices,
   const ParameterVector& parameters,
   SensitivityData&       sensitivities)
{
  const unsigned int Np = parameters.size();
  const unsigned int Nq = qoi.size();

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
  this->assemble_qoi_derivative(qoi_indices);

  for (unsigned int j=0; j != Np; ++j)
    {
      // We currently get partial derivatives via central differencing
      Number delta_p = 1e-6;

      // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)

      Number old_parameter = *parameters[j];

      *parameters[j] = old_parameter - delta_p;
      this->assemble_qoi();
      std::vector<Number> qoi_minus = this->qoi;

      *parameters[j] = old_parameter + delta_p;
      this->assemble_qoi();
      std::vector<Number>& qoi_plus = this->qoi;

      std::vector<Number> partialq_partialp(Nq, 0);
      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          partialq_partialp[i] = (qoi_plus[i] - qoi_minus[i]) / (2.*delta_p);

      // Don't leave the parameter changed
      *parameters[j] = old_parameter;

      for (unsigned int i=0; i != Nq; ++i)
        if (qoi_indices.has_index(i))
          sensitivities[i][j] = partialq_partialp[i] +
            this->get_adjoint_rhs(i).dot(this->get_sensitivity_solution(i));
    }

  // All parameters have been reset.
  // We didn't cache the original rhs or matrix for memory reasons,
  // but we can restore them to a state consistent solution -
  // principle of least surprise.
  this->assembly(true, true);
  this->assemble_qoi(qoi_indices);
}



LinearSolver<Number>* ImplicitSystem::get_linear_solver() const
{
  return LinearSolver<Number>::build().release();
}



std::pair<unsigned int, Real> ImplicitSystem::get_linear_solve_parameters() const
{
  return std::make_pair(this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations"),
                        this->get_equation_systems().parameters.get<Real>("linear solver tolerance"));
}



void ImplicitSystem::release_linear_solver(LinearSolver<Number>* s) const
{
  delete s;
}

