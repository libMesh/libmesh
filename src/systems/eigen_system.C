// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/libmesh_config.h"

// Currently, the EigenSystem should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC)

// C++ includes

// Local includes
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/numeric_vector.h"


namespace libMesh
{


// ------------------------------------------------------------
// EigenSystem implementation
EigenSystem::EigenSystem (EquationSystems& es,
                          const std::string& name_in,
                          const unsigned int number_in
                          ) :
  Parent           (es, name_in, number_in),
  matrix_A         (NULL),
  matrix_B         (NULL),
  eigen_solver     (EigenSolver<Number>::build(es.comm())),
  _n_converged_eigenpairs (0),
  _n_iterations           (0),
  _is_generalized_eigenproblem (false),
  _eigen_problem_type (NHEP),
  _eigenproblem_sensitivity_assemble_system_function(NULL),
  _eigenproblem_sensitivity_assemble_system_object(NULL)
{
}



EigenSystem::~EigenSystem ()
{
  // clear data
  this->clear();
}



void EigenSystem::clear ()
{
  // Clear the parent data
  Parent::clear();

  // delete the matricies
  delete matrix_A;
  delete matrix_B;

  // NULL-out the matricies.
  matrix_A = NULL;
  matrix_B = NULL;

  // clear the solver
  eigen_solver->clear();

}


void EigenSystem::set_eigenproblem_type (EigenProblemType ept)
{
  _eigen_problem_type = ept;

  eigen_solver->set_eigenproblem_type(ept);

  // libMesh::out << "The Problem type is set to be: " << std::endl;

  switch (_eigen_problem_type)
    {
    case HEP:
      // libMesh::out << "Hermitian" << std::endl;
      break;

    case NHEP:
      // libMesh::out << "Non-Hermitian" << std::endl;
      break;

    case GHEP:
      // libMesh::out << "Gerneralized Hermitian" << std::endl;
      break;

    case GNHEP:
      // libMesh::out << "Generalized Non-Hermitian" << std::endl;
      break;

    default:
      // libMesh::out << "not properly specified" << std::endl;
      libmesh_error_msg("Unrecognized _eigen_problem_type = " << _eigen_problem_type);
    }
}




void EigenSystem::init_data ()
{
  // initialize parent data
  Parent::init_data();

  // define the type of eigenproblem
  if (_eigen_problem_type == GNHEP || _eigen_problem_type == GHEP)
    _is_generalized_eigenproblem = true;

  // build the system matrix
  matrix_A = SparseMatrix<Number>::build(this->comm()).release();

  this->init_matrices();
}



void EigenSystem::init_matrices ()
{
  DofMap& dof_map = this->get_dof_map();

  dof_map.attach_matrix(*matrix_A);

  // build matrix_B only in case of a
  // generalized problem
  if (_is_generalized_eigenproblem)
    {
      matrix_B = SparseMatrix<Number>::build(this->comm()).release();
      dof_map.attach_matrix(*matrix_B);
    }

  dof_map.compute_sparsity(this->get_mesh());

  // initialize and zero system matrix
  matrix_A->init();
  matrix_A->zero();

  // eventually initialize and zero system matrix_B
  if (_is_generalized_eigenproblem)
    {
      matrix_B->init();
      matrix_B->zero();
    }
}



void EigenSystem::reinit ()
{
  // initialize parent data
  Parent::reinit();

  // Clear the matrices
  matrix_A->clear();

  if (_is_generalized_eigenproblem)
    matrix_B->clear();

  DofMap& dof_map = this->get_dof_map();

  // Clear the sparsity pattern
  dof_map.clear_sparsity();

  // Compute the sparsity pattern for the current
  // mesh and DOF distribution.  This also updates
  // both matrices, \p DofMap now knows them
  dof_map.compute_sparsity(this->get_mesh());

  matrix_A->init();
  matrix_A->zero();

  if (_is_generalized_eigenproblem)
    {
      matrix_B->init();
      matrix_B->zero();
    }
}



void EigenSystem::solve ()
{

  // A reference to the EquationSystems
  EquationSystems& es = this->get_equation_systems();

  // check that necessary parameters have been set
  libmesh_assert (es.parameters.have_parameter<unsigned int>("eigenpairs"));
  libmesh_assert (es.parameters.have_parameter<unsigned int>("basis vectors"));

  if (this->assemble_before_solve)
    // Assemble the linear system
    this->assemble ();

  // Get the tolerance for the solver and the maximum
  // number of iterations. Here, we simply adopt the linear solver
  // specific parameters.
  const Real tol            =
    es.parameters.get<Real>("linear solver tolerance");

  const unsigned int maxits =
    es.parameters.get<unsigned int>("linear solver maximum iterations");

  const unsigned int nev    =
    es.parameters.get<unsigned int>("eigenpairs");

  const unsigned int ncv    =
    es.parameters.get<unsigned int>("basis vectors");

  std::pair<unsigned int, unsigned int> solve_data;

  // call the solver depending on the type of eigenproblem
  if (_is_generalized_eigenproblem)
    {
      //in case of a generalized eigenproblem
      solve_data = eigen_solver->solve_generalized (*matrix_A,*matrix_B, nev, ncv, tol, maxits);

    }

  else
    {
      libmesh_assert (!matrix_B);

      //in case of a standard eigenproblem
      solve_data = eigen_solver->solve_standard (*matrix_A, nev, ncv, tol, maxits);

    }

  this->_n_converged_eigenpairs = solve_data.first;
  this->_n_iterations           = solve_data.second;

}

  
std::pair<unsigned int, Real>
EigenSystem::sensitivity_solve (const ParameterVector& parameters,
                                std::vector<Number>& sens)
  {
    // make sure that eigensolution is already available
    libmesh_assert(_n_converged_eigenpairs);
    
    // the sensitivity is calculated based on the inner product of the left and
    // right eigen vectors.
    //
    //    y^T [A] x - lambda y^T [B] x = 0
    //    where y and x are the left and right eigenvectors
    //    d lambda/dp = (y^T (d[A]/dp - lambda d[B]/dp) x) / (y^T [B] x)
    //
    //    the denominator remain constant for all sensitivity calculations.
    //
    std::vector<Number> denom(_n_converged_eigenpairs, 0.);
    sens.resize(_n_converged_eigenpairs*parameters.size(), 0.);
    std::pair<Real, Real> eig_val;
    
    AutoPtr< NumericVector<Number> > x_right = NumericVector<Number>::build(this->comm()),
    x_left = NumericVector<Number>::build(this->comm()),
    tmp = NumericVector<Number>::build(this->comm());
    x_right->init(this->n_dofs(), this->n_local_dofs(), false, solution->type());
    x_left->init(this->n_dofs(), this->n_local_dofs(), false, solution->type());
    tmp->init(this->n_dofs(), this->n_local_dofs(), false, solution->type());
    
    for (unsigned int i=0; i<_n_converged_eigenpairs; i++)
    {
      switch (_eigen_problem_type) {
        case HEP:
          // right and left eigenvectors are same
          // imaginary part of eigenvector for real matrices is zero
          this->get_eigenpair(i, x_right.get(), NULL);
          denom[i] = x_right->dot(*x_right);               // x^H x
          break;
          
        case GHEP:
          // imaginary part of eigenvector for real matrices is zero
          this->get_eigenpair(i, x_right.get(), NULL);
          matrix_B->vector_mult(*tmp, *x_right);
          denom[i] = x_right->dot(*tmp);                  // x^H B x
          break;
          
        default:
          // to be implemented for the non-Hermitian problems
          libmesh_error();
          break;
      }
    }
    
    unsigned int num;
    for (unsigned int p=0; p<parameters.size(); p++)
    {
      // calculate sensitivity of matrix quantities
      this->solution->zero();
      this->assemble_eigensystem_sensitivity(parameters, p);
      
      // now calculate sensitivity of each eigenvalue for the parameter
      for (unsigned int i=0; i<_n_converged_eigenpairs; i++)
      {
        num = p*_n_converged_eigenpairs+i;
        switch (_eigen_problem_type)
        {
          case HEP:
            eig_val = this->get_eigenpair(i, x_right.get());
            matrix_A->vector_mult(*tmp, *x_right);
            sens[num] = x_right->dot(*tmp);                     // x^H A' x
            sens[num]-= eig_val.first * x_right->dot(*x_right); // - lambda x^H x
            sens[num] /= denom[i];                              // x^H x
            break;

          case GHEP:
            eig_val = this->get_eigenpair(i, x_right.get());
            matrix_A->vector_mult(*tmp, *x_right);
            sens[num] = x_right->dot(*tmp);                 // x^H A' x
            matrix_B->vector_mult(*tmp, *x_right);
            sens[num]-= eig_val.first * x_right->dot(*tmp); // - lambda x^H B' x
            sens[num] /= denom[i];                          // x^H B x
            break;
            
          default:
            // to be implemented for the non-Hermitian problems
            libmesh_error();
            break;
        }
      }
    }
    
    return std::pair<unsigned int, Real> (0, 0.);
  }

  

void EigenSystem::assemble ()
{

  // Assemble the linear system
  Parent::assemble ();

}


std::pair<Real, Real> EigenSystem::get_eigenpair (unsigned int i,
                                                  NumericVector<Number>* vec_re,
                                                  NumericVector<Number>* vec_im)
{
  NumericVector<Number>* sol_re = NULL;
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  libmesh_assert (vec_im == NULL);
#else
  if (_eigen_problem_type == HEP ||
      _eigen_problem_type == GHEP)
    libmesh_assert (vec_im == NULL);
#endif
  
  sol_re = vec_re;
  
  if (sol_re == NULL)
    sol_re = this->solution.get();
  
  // call the eigen_solver get_eigenpair method
  return eigen_solver->get_eigenpair (i, *sol_re, vec_im);
}

  
  
bool EigenSystem::user_eigensystem_sensitivity_assemble(const ParameterVector& parameters,
                                                        const unsigned int p)
{
  bool rval = false;
  
  // Call the user-provided assembly function,
  // if it was provided
  if (_eigenproblem_sensitivity_assemble_system_function != NULL)
    rval = this->_eigenproblem_sensitivity_assemble_system_function (this->get_equation_systems(),
                                                                     this->name(),
                                                                     parameters,
                                                                     p,
                                                                     matrix_A,
                                                                     matrix_B);
  
  // ...or the user-provided assembly object.
  else if (_eigenproblem_sensitivity_assemble_system_object != NULL)
    rval = this->_eigenproblem_sensitivity_assemble_system_object->sensitivity_assemble
    (parameters, p, matrix_A, matrix_B);
  

  return rval;
}

  

  
void EigenSystem::assemble_eigensystem_sensitivity(const ParameterVector& parameters,
                                                   const unsigned int p)
{
  // if this has not been calculated by the use provided routines, then
  // use central differencing to calculate the sensitivity matrices
  if (!this->user_eigensystem_sensitivity_assemble(parameters, p))
  {
    libMesh::err << "Error: EigenSystem is currently not setup to calculate matrix "
    << "sensitivity using finite differencing!";
    libmesh_error();
  }
}

  
void EigenSystem::attach_eigenproblem_sensitivity_assemble_function
  (bool fptr(EquationSystems& es,
             const std::string& name,
             const ParameterVector& parameters,
             const unsigned int i,
             SparseMatrix<Number>* sensitivity_A,
             SparseMatrix<Number>* sensitivity_B))
{
  libmesh_assert(fptr);
  
  if (_eigenproblem_sensitivity_assemble_system_object != NULL)
  {
    libmesh_here();
    libMesh::out << "WARNING:  Cannot specify both assembly sensitivity function and object!"
    << std::endl;
    
    _eigenproblem_sensitivity_assemble_system_object = NULL;
  }
  
  _eigenproblem_sensitivity_assemble_system_function = fptr;
}
  

void EigenSystem::attach_eigenproblem_sensitivity_assemble_object (EigenproblemSensitivityAssembly& assemble)
{
  if (_eigenproblem_sensitivity_assemble_system_function != NULL)
  {
    libmesh_here();
    libMesh::out << "WARNING:  Cannot specify both assembly object and function!"
    << std::endl;
    
    _eigenproblem_sensitivity_assemble_system_function = NULL;
  }
  
  _eigenproblem_sensitivity_assemble_system_object = &assemble;
}

  
} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC
